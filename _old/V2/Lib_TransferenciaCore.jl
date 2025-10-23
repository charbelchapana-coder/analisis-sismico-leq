module Lib_TransferenciaCore

using LinearAlgebra, FFTW

# Estructuras necesarias (restauradas desde Lib_LEQ.jl)
struct parametros_degradacion
    # Parámetros para curva de degradación G/Go
    a_deg::Float64
    γref_deg::Float64
    
    # Parámetros para curva de amortiguamiento
    a_amor::Float64
    b_amor::Float64
    c_amor::Float64
    
    # Tipo de comportamiento: "lineal" o "equivalente"
    tipo::String
end

# Definición ampliada de estructuras para capas de suelo y roca
mutable struct capa_suelo
    # Propiedades básicas
    ρ::Float64                           # Densidad [kg/m³]
    Vs_inicial::Float64                  # Velocidad de corte inicial [m/s]
    ζ_inicial::Float64                   # Amortiguamiento inicial
    h::Float64                           # Espesor [m]
    
    # Propiedades dinámicas actuales (para análisis lineal equivalente)
    Vs_actual::Float64                   # Velocidad de corte actual [m/s]
    ζ_actual::Float64                    # Amortiguamiento actual
    G_actual::Float64                    # Módulo de corte actual [Pa]
    
    # Parámetros de degradación (opcional para comportamiento no lineal)
    params_degradacion::Union{parametros_degradacion, Nothing}
    
    # Identificación
    id::String                           # Identificador de la capa
    
    # Constructor interno para asegurar consistencia
    function capa_suelo(ρ, Vs_inicial, ζ_inicial, h, params_degradacion=nothing, id="")
        G_inicial = ρ * Vs_inicial^2
        new(ρ, Vs_inicial, ζ_inicial, h, Vs_inicial, ζ_inicial, G_inicial, params_degradacion, id)
    end
end

# Constructor para compatibilidad con código existente
function capa_suelo(ρ::Float64, Vs::Float64, ζ::Float64, h::Float64)
    return capa_suelo(ρ, Vs, ζ, h, nothing, "")
end

export capa_suelo, parametros_degradacion, matriz_T, funcion_transferencia, deformaciones_corte_shake, 
       curva_degradacion, curva_amortiguamiento, actualizar_propiedades_dinamicas!,
       calc_waves_pystrata, calc_strain_tf_pystrata

# Función para calcular la matriz de transferencia de una capa de suelo a otra
function matriz_T(suelo_1::capa_suelo, suelo_2::capa_suelo, ω::Float64)
    # Parámetros de la primera capa (usar propiedades actuales)
    ρ1, Vs1, ζ1, h1 = suelo_1.ρ, suelo_1.Vs_actual, suelo_1.ζ_actual, suelo_1.h
    # Parámetros de la segunda capa (usar propiedades actuales)
    ρ2, Vs2, ζ2, h2 = suelo_2.ρ, suelo_2.Vs_actual, suelo_2.ζ_actual, suelo_2.h

    # Módulo de corte dinámico complejo
    G1 = ρ1 * Vs1^2 * (1 + 2im * ζ1)
    G2 = ρ2 * Vs2^2 * (1 + 2im * ζ2)

    # Número de onda complejo
    k1 = ρ1 * ω^2 / G1
    k2 = ρ2 * ω^2 / G2

    # Impedancia compleja
    αm = (k1*G1)/(k2*G2)

    # Funciones auxiliares
    epos = 0.5*exp(im*k1*h1)
    eneg = 0.5*exp(-im*k1*h1)

    # Matriz de transferencia 
    M = [(1+αm)*epos (1-αm)*eneg;
          (1-αm)*epos (1+αm)*eneg]
    return M
end

"""
    calc_waves_pystrata(capas, omega)

Calcula las amplitudes de ondas ascendentes (A) y descendentes (B) siguiendo la metodología de PyStrata.
Implementación basada en propagation.py líneas 408-444.
"""
function calc_waves_pystrata(capas, omega)
    n_capas = length(capas)
    n_freq = length(omega)
    
    # Inicializar matrices de ondas
    waves_a = ones(ComplexF64, n_capas, n_freq)  # Ondas ascendentes
    waves_b = ones(ComplexF64, n_capas, n_freq)  # Ondas descendentes
    wave_nums = zeros(ComplexF64, n_capas, n_freq)  # Números de onda
    
    # Calcular números de onda complejos para cada capa
    for i = 1:n_capas
        capa = capas[i]
        Vs_comp = capa.Vs_actual * (1 + 2im * capa.ζ_actual)
        for j = 1:n_freq
            if abs(omega[j]) > 1e-12
                wave_nums[i, j] = omega[j] / Vs_comp
            else
                wave_nums[i, j] = 0.0
            end
        end
    end
    
    # Propagación de ondas desde la superficie hacia abajo
    # En la superficie: amplitudes = 1 (reflexión total)
    for i = 1:(n_capas-1)
        # Módulos de corte complejos
        G1 = capas[i].ρ * capas[i].Vs_actual^2 * (1 + 2im * capas[i].ζ_actual)
        G2 = capas[i+1].ρ * capas[i+1].Vs_actual^2 * (1 + 2im * capas[i+1].ζ_actual)
        
        for j = 1:n_freq
            # Impedancia compleja (evitar división por cero)
            if abs(wave_nums[i+1, j] * G2) > 1e-12
                cimped = (wave_nums[i, j] * G1) / (wave_nums[i+1, j] * G2)
            else
                cimped = 1.0
            end
            
            # Término complejo usando espesor total de la capa
            cterm = 1im * wave_nums[i, j] * capas[i].h
            
            # Propagación de ondas (ecuaciones de PyStrata líneas 434-441)
            if isfinite(cimped)
                waves_a[i+1, j] = 0.5 * waves_a[i, j] * (1 + cimped) * exp(cterm) + 
                                  0.5 * waves_b[i, j] * (1 - cimped) * exp(-cterm)
                waves_b[i+1, j] = 0.5 * waves_a[i, j] * (1 - cimped) * exp(cterm) + 
                                  0.5 * waves_b[i, j] * (1 + cimped) * exp(-cterm)
            else
                # Para frecuencias muy bajas o impedancias infinitas
                waves_a[i+1, j] = 1.0
                waves_b[i+1, j] = 1.0
            end
        end
    end
    
    # Manejar frecuencias cercanas a cero
    for i = 1:n_capas
        for j = 1:n_freq
            if abs(omega[j]) < 1e-12
                waves_a[i, j] = 1.0
                waves_b[i, j] = 1.0
            end
        end
    end
    
    return waves_a, waves_b, wave_nums
end

"""
    calc_strain_tf_pystrata(waves_a, waves_b, wave_nums, omega, capa_input_idx, capa_output_idx, depth_within_output)

Calcula la función de transferencia de deformaciones siguiendo la metodología exacta de PyStrata.
Implementación basada en propagation.py líneas 476-503.

Parámetros:
- waves_a, waves_b: Amplitudes de ondas ascendentes y descendentes
- wave_nums: Números de onda complejos
- omega: Frecuencias angulares
- capa_input_idx: Índice de la capa de entrada (típicamente la última: roca)
- capa_output_idx: Índice de la capa de salida
- depth_within_output: Profundidad dentro de la capa de salida (típicamente h/2)
"""
function calc_strain_tf_pystrata(waves_a, waves_b, wave_nums, omega, capa_input_idx, capa_output_idx, depth_within_output)
    n_freq = length(omega)
    strain_tf = zeros(ComplexF64, n_freq)
    
    for j = 1:n_freq
        if abs(omega[j]) > 1e-12
            # Término complejo para la profundidad dentro de la capa
            cterm = 1im * wave_nums[capa_output_idx, j] * depth_within_output
            
            # Numerador: i k*_m [A_m exp(i k*_m h_m/2) - B_m exp(-i k*_m h_m/2)]
            # Siguiendo exactamente PyStrata líneas 495-501
            numer = (1im * wave_nums[capa_output_idx, j] * 
                    (waves_a[capa_output_idx, j] * exp(cterm) - 
                     waves_b[capa_output_idx, j] * exp(-cterm)))
            
            # Denominador: -angFreq^2 * (2 * A_n) 
            # Donde A_n es la amplitud de onda ascendente en la capa de entrada
            denom = -(omega[j]^2) * (2 * waves_a[capa_input_idx, j])
            
            # Función de transferencia con escalamiento de gravedad (PyStrata línea 502)
            GRAVITY = 9.81  # m/s² (constante de PyStrata)
            if abs(denom) > 1e-12
                strain_tf[j] = GRAVITY * numer / denom
            else
                strain_tf[j] = 0.0
            end
        else
            strain_tf[j] = 0.0
        end
    end
    
    return strain_tf
end

function funcion_transferencia(suelos::Vector{capa_suelo}, ω::Vector{Float64})
    H = zeros(ComplexF64, 2, length(ω)) # Inicializar vector para todos los resultados
    Amp = zeros(ComplexF64, length(ω))

    for i in 1:length(ω)
        # Inicialización de la matriz de transferencia total
        M_total = Matrix{ComplexF64}(I, 2, 2) # Matriz identidad compleja
        N = length(suelos)

        # Cálculo de la matriz de transferencia total a través de las capas de suelo
        for j in 1:(length(suelos)-1)
            M = matriz_T(suelos[j], suelos[j+1], ω[i])
            M_total = M * M_total
        end
        H[1,i] = M_total[1,1] + M_total[1,2]
        H[2,i] = M_total[2,1] + M_total[2,2]
        denom = 2*H[1,i]
        if abs(denom) < 1e-12
            Amp[i] = 1.0  # Cambio: evitar Inf, usar 1.0
        else
            Amp[i] = (H[1,i] + H[2,i]) / denom
        end
    end
    return Amp
end

"""
    deformaciones_corte_shake(capas, omega, aceleracion_base_fft)

Calcula las deformaciones de corte máximas usando el método estándar SHAKE/STRATA.
NUEVA IMPLEMENTACIÓN: Siguiendo exactamente la metodología de PyStrata para corregir 
la inversión física de deformaciones (suelo blando debe tener mayor deformación que suelo rígido).
"""
function deformaciones_corte_shake(capas, omega, aceleracion_base_fft)
    n_capas = length(capas) - 1  # Excluir capa de roca
    gamma_max = zeros(n_capas)
    gamma_interfaces = zeros(max(n_capas - 1, 0))

    # Evitar división por cero
    omega_safe = [max(abs(w), 1e-6) for w in omega]

    # 1. Calcular amplitudes de ondas A y B usando metodología PyStrata
    waves_a, waves_b, wave_nums = calc_waves_pystrata(capas, omega)
    
    # 2. Calcular deformaciones usando función de transferencia de PyStrata
    capa_input_idx = length(capas)  # Roca basal (última capa)
    
    for i in 1:n_capas
        # Calcular deformación en el centro de la capa (h/2)
        depth_within = capas[i].h / 2.0
        
        # Función de transferencia de deformaciones desde roca hasta centro de capa i
        strain_tf = calc_strain_tf_pystrata(waves_a, waves_b, wave_nums, omega, 
                                          capa_input_idx, i, depth_within)
        
        # Calcular deformación máxima en dominio de tiempo
        # Aplicar función de transferencia a la aceleración de entrada
        strain_fas = strain_tf .* aceleracion_base_fft
        
        # Convertir a tiempo usando simetría Hermitiana
        strain_full = [strain_fas; conj(reverse(strain_fas[2:end-1]))]
        strain_time = real(ifft(strain_full))
        
        # Deformación máxima con factor de tensión efectiva estándar
        # Factor 0.65 es el estándar en análisis lineal equivalente (Idriss & Sun, 1992)
        strain_ratio = 0.65
        gamma_max[i] = maximum(abs.(strain_time)) * strain_ratio * 100.0  # Convertir a %
    end

    # Deformación en superficie (condición de frontera libre: γ = 0)
    gamma_superficie = 0.0

    # Deformación en roca basal usando velocidad de partícula
    vel_base_fft = zeros(ComplexF64, length(omega))
    for i in 2:length(omega)
        vel_base_fft[i] = 1im * omega[i] * (-aceleracion_base_fft[i] / (omega_safe[i]^2))
    end
    vel_base_full = [vel_base_fft; conj(reverse(vel_base_fft[2:end-1]))]
    vel_base_time = real(ifft(vel_base_full))
    v_base_max = maximum(abs.(vel_base_time))
    Vs_roca = capas[end].Vs_actual
    gamma_roca = (v_base_max / Vs_roca) * 0.65 * 100.0  # Factor estándar para roca

    # Deformaciones en interfaces entre capas
    for i in 1:(n_capas - 1)
        gamma_interfaces[i] = (gamma_max[i] + gamma_max[i+1]) / 2.0
    end

    return gamma_max, gamma_superficie, gamma_roca, gamma_interfaces
end

"""
    curva_degradacion(γ, params::parametros_degradacion)

Calcula el factor de degradación del módulo de corte G/Go.
"""
function curva_degradacion(γ, params::parametros_degradacion)
    if params.tipo == "lineal"
        return 1.0  # Sin degradación para comportamiento lineal
    else
        # Curva hiperbólica típica
        return 1.0 / (1 + params.a_deg * (γ / params.γref_deg))
    end
end

"""
    curva_amortiguamiento(γ, params::parametros_degradacion)

Calcula el amortiguamiento en función de la deformación.
"""
function curva_amortiguamiento(γ, params::parametros_degradacion)
    if params.tipo == "lineal"
        return params.a_amor  # Valor constante para comportamiento lineal
    else
        # Curva típica de amortiguamiento
        return params.a_amor + params.b_amor * (γ / params.γref_deg) / (1 + params.c_amor * (γ / params.γref_deg))
    end
end

"""
    actualizar_propiedades_dinamicas!(capa::capa_suelo, γ::Float64)

Actualiza las propiedades dinámicas de una capa basándose en la deformación.
"""
function actualizar_propiedades_dinamicas!(capa::capa_suelo, γ::Float64)
    if capa.params_degradacion === nothing
        # Sin parámetros de degradación, mantener valores iniciales
        capa.Vs_actual = capa.Vs_inicial
        capa.ζ_actual = capa.ζ_inicial
        capa.G_actual = capa.ρ * capa.Vs_actual^2
    else
        # Aplicar degradación
        factor_G = curva_degradacion(γ, capa.params_degradacion)
        ζ_nuevo = curva_amortiguamiento(γ, capa.params_degradacion)
        
        # Actualizar propiedades
        capa.Vs_actual = capa.Vs_inicial * sqrt(factor_G)
        capa.ζ_actual = ζ_nuevo
        capa.G_actual = capa.ρ * capa.Vs_actual^2
    end
end

end  # module