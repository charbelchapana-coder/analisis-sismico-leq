module Lib_LElastico

using LinearAlgebra, FFTW

export capa_suelo, curva_degradacion, curva_amortiguamiento, matriz_T, funcion_transferencia, espectro_chopra, deformaciones_corte_shake

# Definición de estructuras para capas de suelo y roca
struct capa_suelo
    ρ::Float64
    Vs::Float64
    ζ::Float64
    h::Float64
end

# Definición de función de degradacion de rigidez (G/Go)
function curva_degradacion(γ, a, γref)
    if γ > 1
        Ge = 1/(1+(1/γref)^a)
    else
        Ge = 1/(1+(γ/γref)^a)
    end
    return Ge
end

# Definición de función de amortiguamiento (D)
function curva_amortiguamiento(γ, a, b, c)
    if γ > 1
        De = a + b + c
    else
        De = a*γ^2 + b*γ + c
    end
    return De
end

# Función para calcular la matriz de transferencia de una capa de suelo a otra
function matriz_T(suelo_1::capa_suelo, suelo_2::capa_suelo, ω::Float64)
    # Parámetros de la primera capa
    ρ1, Vs1, ζ1, h1 = suelo_1.ρ, suelo_1.Vs, suelo_1.ζ, suelo_1.h
    # Parámetros de la segunda capa
    ρ2, Vs2, ζ2, h2 = suelo_2.ρ, suelo_2.Vs, suelo_2.ζ, suelo_2.h

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

function funcion_transferencia(suelos::Vector{capa_suelo}, ω::Vector{Float64})
    H = zeros(ComplexF64, 2, length(ω)) # Inicializar vector para todos los resultados
    Amp = zeros(ComplexF64, length(ω))

    for i in 1:length(ω)
        # Inicialización de la matriz de transferencia total
        M_total = Matrix{ComplexF64}(I, 2, 2) # Corregido: matriz identidad compleja
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
            Amp[i] = 0.0
        else
            Amp[i] = (H[1,i] + H[2,i]) / denom
        end
    end
    return Amp
end

"""
    deformaciones_corte_shake(capas, omega, aceleracion_base_fft)

Calcula las deformaciones de corte máximas en el centro de cada capa usando 
la función de transferencia y el gradiente de desplazamientos.

# Argumentos
- `capas`: Vector de objetos `capa_suelo` 
- `omega`: Vector de frecuencias angulares (rad/s)
- `aceleracion_base_fft`: FFT de la aceleración en la base (roca)

# Retorna
- `gamma_max`: Vector con las deformaciones de corte máximas para cada capa (%)
"""
function deformaciones_corte_shake(capas, omega, aceleracion_base_fft)
    n_capas = length(capas) - 1  # Excluir capa de roca
    gamma_max = zeros(n_capas)
    
    # Desplazamientos en frecuencia: u(ω) = -a(ω)/ω²
    despl_base_fft = -aceleracion_base_fft ./ (omega.^2)
    despl_base_fft[1] = 0.0  # Evitar división por cero
    
    # Calcular deformaciones en el centro de cada capa
    for i in 1:n_capas
        h = capas[i].h
        Vs = capas[i].Vs
        
        # Calcular función de transferencia hasta el centro de esta capa
        # Crear una capa ficticia con altura hasta el centro
        capas_hasta_centro = copy(capas[1:i])
        capas_hasta_centro[i] = capa_suelo(capas[i].ρ, capas[i].Vs, capas[i].ζ, h/2)
        
        H_centro_capa = funcion_transferencia(capas_hasta_centro, omega)
        
        # Desplazamiento en el centro de la capa
        despl_centro_capa_fft = H_centro_capa .* despl_base_fft
        
        # Número de onda de la capa
        k_capa = omega ./ Vs
        
        # Deformación de corte en el centro de la capa
        gamma_freq_capa = 1im .* k_capa .* despl_centro_capa_fft
        
        # Transformada inversa
        gamma_full_capa = [gamma_freq_capa; conj(reverse(gamma_freq_capa[2:end-1]))]
        gamma_t_capa = real(ifft(gamma_full_capa))
        
        # Tomar valor absoluto y encontrar máximo
        gamma_max[i] = maximum(abs.(gamma_t_capa)) * 100
    end
    
    return gamma_max
end

function espectro_chopra(aceleracion, dt, periodos, xi)
    Sa = zeros(length(periodos))
    m = 1.0
    for k in eachindex(periodos)
        T = periodos[k]
        if T == 0
            Sa[k] = maximum(abs.(aceleracion))
            continue
        end
        w_n = 2 * π / T
        K = w_n^2 * m
        w_d = w_n * sqrt(1 - xi^2)
        e_xi_wn_dt = exp(-xi * w_n * dt)
        sin_wd_dt = sin(w_d * dt)
        cos_wd_dt = cos(w_d * dt)

        A = e_xi_wn_dt * (xi / sqrt(1 - xi^2) * sin_wd_dt + cos_wd_dt)
        B = e_xi_wn_dt * (1 / w_d * sin_wd_dt)
        C = (1 / K) * (2 * xi / (w_n * dt) + e_xi_wn_dt * (((1 - 2 * xi^2) / (w_d * dt) - xi / sqrt(1 - xi^2)) * sin_wd_dt - (1 + 2 * xi / (w_n * dt)) * cos_wd_dt))
        D = (1 / K) * (1 - 2 * xi / (w_n * dt) - e_xi_wn_dt * (((1 - 2 * xi^2) / (w_d * dt)) * sin_wd_dt - (2 * xi / (w_n * dt)) * cos_wd_dt))

        Aprime = -w_n / sqrt(1 - xi^2) * e_xi_wn_dt * sin_wd_dt
        Bprime = e_xi_wn_dt * (cos_wd_dt - xi / sqrt(1 - xi^2) * sin_wd_dt)
        Cprime = (1 / K) * (-1 / dt + e_xi_wn_dt * ((w_n / sqrt(1 - xi^2) + xi / (dt * sqrt(1 - xi^2))) * sin_wd_dt + (1 / dt) * cos_wd_dt))
        Dprime = (1 / (K * dt)) * (1 - e_xi_wn_dt * ((xi / sqrt(1 - xi^2)) * sin_wd_dt + cos_wd_dt))

        u = zeros(length(aceleracion))
        u_dot = zeros(length(aceleracion))
        p = -m .* aceleracion

        for i in 1:(length(aceleracion) - 1)
            u[i+1] = A * u[i] + B * u_dot[i] + C * p[i] + D * p[i+1]
            u_dot[i+1] = Aprime * u[i] + Bprime * u_dot[i] + Cprime * p[i] + Dprime * p[i+1]
        end

        acel_total = -(2 * xi * w_n .* u_dot + w_n^2 .* u)
        Sa[k] = maximum(abs.(acel_total))
    end
    return Sa
end

end # module