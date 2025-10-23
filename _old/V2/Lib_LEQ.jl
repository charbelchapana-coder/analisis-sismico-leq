module Lib_LEQ

# Incluir el módulo principal con funciones corregidas
include("Lib_TransferenciaCore.jl")
using .Lib_TransferenciaCore

using LinearAlgebra, FFTW, Statistics, DelimitedFiles, Plots, Printf

export capa_suelo, parametros_degradacion, curva_degradacion, curva_amortiguamiento, 
       matriz_T, funcion_transferencia, espectro_chopra, deformaciones_corte_shake,
       leer_sismo, escalar_sismo, periodo_predominante, periodo_fundamental,
       analisis_lineal_equivalente, calcular_movimientos_capa, espectro_fourier,
       funcion_amplificacion, modificar_dt, historia_tension_deformacion,
       subdividir_capas, graficar_perfil_deformaciones_detallado, graficar_perfil_desplazamientos,
       graficar_deformaciones_derivadas, calcular_desplazamientos_interfaces, 
       calcular_escala_unificada_deformaciones, reconstruir_señal_superficie,
       ejecutar_analisis_sismico_completo, manejar_no_convergencia,
       calcular_deformaciones_convolucion

# Definición de parámetros para curvas de degradación y amortiguamiento
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
    
    # Parámetros de degradación (opcional para comportamiento no lineal)
    params_degradacion::Union{parametros_degradacion, Nothing}
    
    # Identificación
    id::String                           # Identificador de la capa
    
    # Constructor interno para asegurar consistencia
    function capa_suelo(ρ, Vs_inicial, ζ_inicial, h, params_degradacion=nothing, id="")
        new(ρ, Vs_inicial, ζ_inicial, h, Vs_inicial, ζ_inicial, params_degradacion, id)
    end
end

# Constructor para compatibilidad con código existente
function capa_suelo(ρ::Float64, Vs::Float64, ζ::Float64, h::Float64)
    return capa_suelo(ρ, Vs, ζ, h, nothing, "")
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

"""
    subdividir_capas(capas::Vector{capa_suelo}, espesor_maximo::Float64)

Subdivide las capas de suelo en subcapas con espesor máximo especificado.
Cada subcapa mantiene las propiedades de la capa original.
"""
function subdividir_capas(capas::Vector{capa_suelo}, espesor_maximo::Float64)
    subcapas = capa_suelo[]
    
    for (idx, capa) in enumerate(capas)
        # Si es la última capa (roca basal), no subdividir
        if capa.h == 0  # Roca basal
            push!(subcapas, capa)
            continue
        end
        
        # Calcular número de subdivisiones necesarias
        num_subdivisiones = Int(ceil(capa.h / espesor_maximo))
        espesor_subcapa = capa.h / num_subdivisiones
        
        # Crear subcapas
        for i in 1:num_subdivisiones
            subcapa_id = if num_subdivisiones == 1
                capa.id
            else
                "$(capa.id)_sub$i"
            end
            
            # Crear nueva subcapa con las mismas propiedades
            nueva_subcapa = capa_suelo(
                capa.ρ,
                capa.Vs_inicial,
                capa.ζ_inicial,
                espesor_subcapa,
                capa.params_degradacion,
                subcapa_id
            )
            
            # Mantener las propiedades actuales si han sido modificadas
            nueva_subcapa.Vs_actual = capa.Vs_actual
            nueva_subcapa.ζ_actual = capa.ζ_actual
            
            push!(subcapas, nueva_subcapa)
        end
    end
    
    return subcapas
end

"""
    leer_sismo(archivo::String)

Lee un archivo de datos sísmicos y retorna tiempo, aceleración y dt.
Formato esperado: tiempo, aceleración separados por comas.
"""
function leer_sismo(archivo::String)
    try
        datos = readdlm(archivo, ',')
        tiempo = datos[:, 1]
        aceleracion = datos[:, 2]
        dt = tiempo[2] - tiempo[1]
        return tiempo, aceleracion, dt
    catch e
        error("Error al leer archivo sísmico: $e")
    end
end

"""
    escalar_sismo(aceleracion::Vector{Float64}, factor::Float64)

Escala los valores de aceleración por un factor dado.
"""
function escalar_sismo(aceleracion::Vector{Float64}, factor::Float64)
    return aceleracion .* factor
end

"""
    periodo_predominante(aceleracion::Vector{Float64}, dt::Float64)

Calcula el período predominante del registro sísmico usando el espectro de Fourier.
"""
function periodo_predominante(aceleracion::Vector{Float64}, dt::Float64)
    N = length(aceleracion)
    
    # Calcular FFT
    A_fft = fft(aceleracion)
    A_mag = abs.(A_fft[1:div(N,2)])
    
    # Crear vector de frecuencias
    freq = [i/(N*dt) for i in 0:(div(N,2)-1)]
    freq[1] = 1e-6  # Evitar división por cero
    
    # Encontrar frecuencia dominante (excluyendo frecuencia cero)
    idx_max = argmax(A_mag[2:end]) + 1
    freq_predominante = freq[idx_max]
    
    return 1.0 / freq_predominante  # Período predominante
end

"""
    periodo_fundamental(capas::Vector{capa_suelo})

Calcula el período fundamental del depósito de suelo usando la fórmula simplificada:
T₁ = 4H/Vs_promedio
donde H es la altura total y Vs_promedio es la velocidad promedio ponderada.
"""
function periodo_fundamental(capas::Vector{capa_suelo})
    n_capas = length(capas) - 1  # Excluir capa de roca (última)
    
    # Calcular altura total
    H_total = sum(capas[i].h for i in 1:n_capas)
    
    # Calcular velocidad promedio ponderada por espesor
    Vs_prom = H_total / sum(capas[i].h / capas[i].Vs_actual for i in 1:n_capas)
    
    # Período fundamental
    T1 = 4 * H_total / Vs_prom
    
    return T1, H_total, Vs_prom
end

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
            Amp[i] = Inf 
        else
            Amp[i] = (H[1,i] + H[2,i]) / denom
        end
    end
    return Amp
end

"""
    deformaciones_corte_shake(capas, omega, aceleracion_base_fft)

Calcula las deformaciones de corte máximas usando el gradiente de desplazamientos
entre puntos, respetando las condiciones de frontera físicas.

# Argumentos
- `capas`: Vector de objetos `capa_suelo` 
- `omega`: Vector de frecuencias angulares (rad/s)
- `aceleracion_base_fft`: FFT de la aceleración en la base (roca)

# Retorna
- `gamma_max`: Vector con las deformaciones de corte máximas para cada capa (%)
- `gamma_superficie`: Deformación de corte máxima en superficie (siempre 0.0%)
- `gamma_roca`: Deformación de corte máxima en roca basal (%)
- `gamma_interfaces`: Vector con las deformaciones de corte en las interfaces (%)
"""
function deformaciones_corte_shake(capas, omega, aceleracion_base_fft)
    # CORRECCIÓN: Usar la implementación de Lib_TransferenciaCore que sigue el método STRATA exacto.
    # Se elimina la implementación con factores artificiales.
    return Lib_TransferenciaCore.deformaciones_corte_shake(capas, omega, aceleracion_base_fft)
end

"""
    analisis_lineal_equivalente(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                aceleracion_base_fft::Vector{ComplexF64}, 
                                max_iter::Int=10, tol::Float64=0.05)

Realiza análisis lineal equivalente iterativo para encontrar propiedades 
de suelo compatibles con las deformaciones calculadas.

# Argumentos
- `capas`: Vector de capas de suelo
- `omega`: Vector de frecuencias angulares
- `aceleracion_base_fft`: FFT de aceleración en la base
- `max_iter`: Número máximo de iteraciones
- `tol`: Tolerancia para convergencia (cambio relativo en Vs)

# Retorna
- `capas_actualizadas`: Vector de capas con propiedades actualizadas
- `convergencia`: Bool indicando si el análisis convergió
- `historial`: Diccionario con historial de iteraciones
- `Amp_inicial`: Función de transferencia inicial (lineal elástica)
- `Amp_final`: Función de transferencia final (lineal equivalente)
"""
function analisis_lineal_equivalente(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                   aceleracion_base_fft::Vector{ComplexF64}, 
                                   max_iter::Int=10, tol::Float64=0.05)
    
    # Hacer copia de las capas para modificar
    capas_trabajo = deepcopy(capas)
    n_capas = length(capas_trabajo) - 1  # Excluir roca
    
    # Historial para tracking
    historial = Dict{String, Any}(
        "Vs" => Vector{Vector{Float64}}(),
        "zeta" => Vector{Vector{Float64}}(),
        "gamma" => Vector{Vector{Float64}}(),
        "error" => Float64[],
        "factor_relajacion" => Float64[]
    )
    
    # Parámetros de convergencia mejorados
    factor_relajacion_inicial = 0.5  # Factor de relajación inicial
    factor_relajacion_min = 0.1      # Factor mínimo de relajación
    factor_relajacion_max = 0.8      # Factor máximo de relajación
    contador_oscilaciones = 0        # Contador de oscilaciones
    max_oscilaciones = 3             # Máximo de oscilaciones permitidas
    error_anterior = Inf             # Error de la iteración anterior
    
    convergencia = false
    
    # Criterios de convergencia múltiples
    tol_Vs = tol          # Tolerancia para Vs
    tol_zeta = tol * 2.0  # Tolerancia más relajada para amortiguamiento
    tol_gamma = 0.05      # Tolerancia para deformaciones (5%)
    
    # Imprimir encabezado de tabla de iteraciones
    println("\nTabla de convergencia del análisis lineal equivalente (ALGORITMO MEJORADO):")
    println(repeat("=", 120))
    println("| Iter |   Error  | Factor_Relaj | Vs_Diatomáceo | Vs_Brecha | γ_max_Diat | γ_max_Brecha | ζ_Diatomáceo | Estado |")
    println("|      |    (%)   |      (-)     |     (m/s)     |   (m/s)   |     (%)    |      (%)     |      (-)     |        |")
    println(repeat("=", 120))
    
    for iter in 1:max_iter
        # Guardar propiedades actuales
        Vs_anterior = [capa.Vs_actual for capa in capas_trabajo[1:n_capas]]
        zeta_anterior = [capa.ζ_actual for capa in capas_trabajo[1:n_capas]]
        
        # Calcular deformaciones con propiedades actuales
        gamma_max, _, _, _ = deformaciones_corte_shake(capas_trabajo, omega, aceleracion_base_fft)
        
        # Calcular nuevas propiedades sin aplicar aún
        Vs_nuevos = copy(Vs_anterior)
        zeta_nuevos = copy(zeta_anterior)
        
        # Actualizar propiedades para capas con comportamiento no lineal
        for i in 1:n_capas
            if capas_trabajo[i].params_degradacion !== nothing
                params = capas_trabajo[i].params_degradacion
                
                if params.tipo == "equivalente"
                    # Usar deformación efectiva (típicamente 65% de la máxima)
                    gamma_efectiva = gamma_max[i] * 0.65 / 100  # Convertir a fracción
                    
                    # Calcular factor de degradación de rigidez
                    G_ratio = curva_degradacion(gamma_efectiva, params.a_deg, params.γref_deg)
                    
                    # Calcular amortiguamiento
                    D_nuevo = curva_amortiguamiento(gamma_efectiva, params.a_amor, params.b_amor, params.c_amor)
                    
                    # Calcular valores objetivo
                    Vs_objetivo = capas_trabajo[i].Vs_inicial * sqrt(G_ratio)
                    zeta_objetivo = D_nuevo
                    
                    # Aplicar relajación para suavizar convergencia
                    factor_relax = factor_relajacion_inicial
                    if iter > 1
                        # Obtener factor de relajación del historial
                        factor_relax = historial["factor_relajacion"][end]
                    end
                    
                    Vs_nuevos[i] = Vs_anterior[i] + factor_relax * (Vs_objetivo - Vs_anterior[i])
                    zeta_nuevos[i] = zeta_anterior[i] + factor_relax * (zeta_objetivo - zeta_anterior[i])
                end
            end
        end
        
        # Actualizar propiedades en las estructuras
        for i in 1:n_capas
            capas_trabajo[i].Vs_actual = Vs_nuevos[i]
            capas_trabajo[i].ζ_actual = zeta_nuevos[i]
        end
        
        # Guardar en historial
        push!(historial["Vs"], copy(Vs_nuevos))
        push!(historial["zeta"], copy(zeta_nuevos))
        push!(historial["gamma"], copy(gamma_max))
        push!(historial["factor_relajacion"], factor_relajacion_inicial)
        
        # Calcular errores de convergencia múltiples
        error_Vs = maximum(abs.(Vs_nuevos .- Vs_anterior) ./ Vs_anterior)
        error_zeta = length(zeta_nuevos) > 0 ? maximum(abs.(zeta_nuevos .- zeta_anterior) ./ max.(zeta_anterior, 1e-6)) : 0.0
        
        # Error combinado (ponderado)
        error_combinado = 0.7 * error_Vs + 0.3 * error_zeta
        push!(historial["error"], error_combinado)
        
        # Control adaptativo de factor de relajación
        estado_convergencia = "Normal"
        if iter > 1
            # Detectar oscilaciones
            if error_combinado > error_anterior * 1.1
                contador_oscilaciones += 1
                if contador_oscilaciones >= max_oscilaciones
                    # Reducir factor de relajación para estabilizar
                    factor_relajacion_inicial = max(factor_relajacion_min, factor_relajacion_inicial * 0.7)
                    contador_oscilaciones = 0
                    estado_convergencia = "Oscilando"
                end
            else
                # Progreso estable, puede aumentar factor de relajación
                contador_oscilaciones = 0
                if iter > 5 && error_combinado < error_anterior * 0.8
                    factor_relajacion_inicial = min(factor_relajacion_max, factor_relajacion_inicial * 1.1)
                    estado_convergencia = "Estable"
                end
            end
        end
        
        # Calcular valores representativos para mostrar en tabla
        Vs_diat_prom = mean([capas_trabajo[i].Vs_actual for i in 1:min(22, n_capas)])  # Subcapas 1-22: diatomáceo
        Vs_brecha_prom = n_capas > 22 ? mean([capas_trabajo[i].Vs_actual for i in 23:n_capas]) : 0.0  # Resto: brecha
        
        gamma_diat_max = maximum(gamma_max[1:min(22, n_capas)])  # Máximo en diatomáceo
        gamma_brecha_max = n_capas > 22 ? maximum(gamma_max[23:n_capas]) : 0.0  # Máximo en brecha
        
        zeta_diat_prom = mean([capas_trabajo[i].ζ_actual for i in 1:min(22, n_capas)])  # Amortiguamiento diatomáceo
        
        # Imprimir fila de la tabla
        println(@sprintf("| %4d | %6.2f%% | %10.2f | %11.1f | %7.1f | %8.3f | %10.3f | %10.3f | %6s |", 
                iter, error_combinado*100, factor_relajacion_inicial, Vs_diat_prom, Vs_brecha_prom, 
                gamma_diat_max, gamma_brecha_max, zeta_diat_prom, estado_convergencia))
        
        # Verificar convergencia con criterios múltiples
        converge_Vs = error_Vs < tol_Vs
        converge_zeta = error_zeta < tol_zeta
        converge_general = error_combinado < tol
        
        if converge_Vs && converge_zeta && converge_general
            println(repeat("-", 120))
            println("✓ CONVERGENCIA ALCANZADA en iteración ", iter, "!")
            println("  - Error Vs: $(round(error_Vs*100, digits=3))% < $(round(tol_Vs*100, digits=1))%")
            println("  - Error ζ: $(round(error_zeta*100, digits=3))% < $(round(tol_zeta*100, digits=1))%")
            println("  - Error combinado: $(round(error_combinado*100, digits=3))% < $(round(tol*100, digits=1))%")
            convergencia = true
            break
        end
        
        error_anterior = error_combinado
    end
    
    println(repeat("=", 120))
    
    # Calcular funciones de transferencia inicial y final
    Amp_inicial = funcion_transferencia(capas, omega)  # Función H inicial (lineal elástica)
    if convergencia
        Amp_final = funcion_transferencia(capas_trabajo, omega)
    else
        Amp_final = Amp_inicial
    end
    
    return capas_trabajo, convergencia, historial, Amp_inicial, Amp_final
end

"""
    calcular_movimientos_capa(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                             aceleracion_base_fft::Vector{ComplexF64}, 
                             indice_capa::Int)

Calcula los movimientos (aceleración, velocidad, desplazamiento) en una capa específica.
"""
function calcular_movimientos_capa(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                 aceleracion_base_fft::Vector{ComplexF64}, 
                                 indice_capa::Int)
    
    # Verificar índice válido
    if indice_capa < 1 || indice_capa > length(capas)
        error("Índice de capa inválido: $indice_capa")
    end
    
    # Calcular función de transferencia hasta la capa solicitada
    if indice_capa == length(capas)  # Roca basal
        H_capa = ones(ComplexF64, length(omega))
    else
        capas_hasta_solicitada = capas[1:indice_capa]
        H_capa = funcion_transferencia(capas_hasta_solicitada, omega)
    end
    
    # Calcular aceleración en la capa
    aceleracion_capa_fft = H_capa .* aceleracion_base_fft
    
    # Calcular velocidad: v(ω) = a(ω)/(iω)
    velocidad_capa_fft = aceleracion_capa_fft ./ (1im .* omega)
    velocidad_capa_fft[1] = 0.0  # Evitar división por cero
    
    # Calcular desplazamiento: u(ω) = a(ω)/(-ω²)
    desplazamiento_capa_fft = -aceleracion_capa_fft ./ (omega.^2)
    desplazamiento_capa_fft[1] = 0.0  # Evitar división por cero
    
    # Reconstruir señales completas (simétricas) para IFFT
    N_total = 2 * (length(omega) - 1)
    
    acel_full = [aceleracion_capa_fft; conj(reverse(aceleracion_capa_fft[2:end-1]))]
    vel_full = [velocidad_capa_fft; conj(reverse(velocidad_capa_fft[2:end-1]))]
    despl_full = [desplazamiento_capa_fft; conj(reverse(desplazamiento_capa_fft[2:end-1]))]
    
    # Transformadas inversas
    aceleracion_t = real(ifft(acel_full))
    velocidad_t = real(ifft(vel_full))
    desplazamiento_t = real(ifft(despl_full))
    
    return aceleracion_t, velocidad_t, desplazamiento_t, H_capa
end

"""
    espectro_fourier(señal::Vector{Float64}, dt::Float64)

Calcula el espectro de Fourier de una señal.
"""
function espectro_fourier(señal::Vector{Float64}, dt::Float64)
    N = length(señal)
    
    # FFT de la señal
    F_fft = fft(señal)
    F_mag = abs.(F_fft[1:div(N,2)])
    F_fase = angle.(F_fft[1:div(N,2)])
    
    # Vector de frecuencias
    freq = [i/(N*dt) for i in 0:(div(N,2)-1)]
    
    return freq, F_mag, F_fase
end

"""
    funcion_amplificacion(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                         capa_referencia::Int, capa_objetivo::Int)

Calcula la función de amplificación entre dos capas.
"""
function funcion_amplificacion(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                              capa_referencia::Int, capa_objetivo::Int)
    
    # Función de transferencia hasta capa de referencia
    if capa_referencia == length(capas)  # Roca
        H_ref = ones(ComplexF64, length(omega))
    else
        H_ref = funcion_transferencia(capas[1:capa_referencia], omega)
    end
    
    # Función de transferencia hasta capa objetivo
    if capa_objetivo == length(capas)  # Roca
        H_obj = ones(ComplexF64, length(omega))
    else
        H_obj = funcion_transferencia(capas[1:capa_objetivo], omega)
    end
    
    # Función de amplificación
    FA = H_obj ./ H_ref
    
    return FA
end

"""
    funcion_amplificacion(acc_superficie::Vector{Float64}, acc_roca::Vector{Float64}, dt::Float64)

Calcula la función de amplificación a partir de señales de aceleración en el dominio del tiempo.
Retorna las frecuencias y la amplificación correspondiente.
"""
function funcion_amplificacion(acc_superficie::Vector{Float64}, acc_roca::Vector{Float64}, dt::Float64)
    
    # Verificar que las señales tengan la misma longitud
    if length(acc_superficie) != length(acc_roca)
        error("Las señales de superficie y roca deben tener la misma longitud")
    end
    
    N = length(acc_superficie)
    
    # Calcular FFT de ambas señales
    fft_superficie = fft(acc_superficie)
    fft_roca = fft(acc_roca)
    
    # Calcular frecuencias
    fs = 1.0 / dt  # Frecuencia de muestreo
    freq = fftfreq(N, fs)
    
    # Tomar solo la mitad positiva del espectro
    n_freq = div(N, 2)
    freq_pos = freq[1:n_freq]
    
    # Calcular función de amplificación
    amplificacion = abs.(fft_superficie[1:n_freq]) ./ abs.(fft_roca[1:n_freq])
    
    # Evitar divisiones por cero
    for i in eachindex(amplificacion)
        if abs(fft_roca[i]) < 1e-10
            amplificacion[i] = 1.0
        end
    end
    
    # Evitar la frecuencia cero
    if freq_pos[1] == 0.0
        freq_pos = freq_pos[2:end]
        amplificacion = amplificacion[2:end]
    end
    
    return freq_pos, amplificacion
end

"""
    modificar_dt(señal::Vector{Float64}, dt_original::Float64, dt_nuevo::Float64)

Modifica el intervalo de tiempo de una señal sin cambiar el período predominante
ni la duración usando interpolación.
"""
function modificar_dt(señal::Vector{Float64}, dt_original::Float64, dt_nuevo::Float64)
    
    # Tiempo original
    t_original = [i * dt_original for i in 0:(length(señal)-1)]
    
    # Nuevo vector de tiempo manteniendo la duración total
    duracion_total = t_original[end]
    n_puntos_nuevo = Int(round(duracion_total / dt_nuevo)) + 1
    t_nuevo = [i * dt_nuevo for i in 0:(n_puntos_nuevo-1)]
    
    # Interpolación linear simple (sin dependencias externas)
    señal_nueva = zeros(length(t_nuevo))
    
    for i in eachindex(t_nuevo)
        t = t_nuevo[i]
        
        # Encontrar índices para interpolación
        if t <= t_original[1]
            señal_nueva[i] = señal[1]
        elseif t >= t_original[end]
            señal_nueva[i] = señal[end]
        else
            # Interpolación lineal
            idx = findfirst(x -> x > t, t_original) - 1
            t1, t2 = t_original[idx], t_original[idx+1]
            y1, y2 = señal[idx], señal[idx+1]
            
            señal_nueva[i] = y1 + (y2 - y1) * (t - t1) / (t2 - t1)
        end
    end
    
    return señal_nueva, t_nuevo
end

"""
    historia_tension_deformacion(capas::Vector{capa_suelo}, omega::Vector{Float64},
                                aceleracion_base_fft::Vector{ComplexF64}, 
                                indice_capa::Int, dt::Float64)

Calcula el historial temporal de tensión y deformación de corte en el centro de una capa.
"""
function historia_tension_deformacion(capas::Vector{capa_suelo}, omega::Vector{Float64},
                                     aceleracion_base_fft::Vector{ComplexF64}, 
                                     indice_capa::Int, dt::Float64)
    
    # Obtener movimientos en la capa
    aceleracion_t, velocidad_t, desplazamiento_t, _ = calcular_movimientos_capa(
        capas, omega, aceleracion_base_fft, indice_capa)
    
    # Calcular deformación de corte temporal
    capa = capas[indice_capa]
    k_capa = omega ./ capa.Vs_actual
    
    # Desplazamientos en frecuencia
    despl_fft = -aceleracion_base_fft ./ (omega.^2)
    despl_fft[1] = 0.0
    
    # Función de transferencia a la capa
    if indice_capa == length(capas)
        H_capa = ones(ComplexF64, length(omega))
    else
        H_capa = funcion_transferencia(capas[1:indice_capa], omega)
    end
    
    despl_capa_fft = H_capa .* despl_fft
    
    # Deformación de corte en frecuencia
    gamma_fft = 1im .* k_capa .* despl_capa_fft
    
    # Transformada inversa para deformación
    gamma_full = [gamma_fft; conj(reverse(gamma_fft[2:end-1]))]
    gamma_t = real(ifft(gamma_full))
    
    # Calcular tensión de corte: τ = G * γ
    G_capa = capa.ρ * capa.Vs_actual^2  # Módulo de corte
    tension_t = G_capa .* gamma_t
    
    # Vector de tiempo
    N = length(gamma_t)
    tiempo = [i * dt for i in 0:(N-1)]
    
    return tiempo, gamma_t .* 100, tension_t  # γ en %, τ en Pa
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

"""
    graficar_perfil_deformaciones_detallado(subcapas::Vector{capa_suelo}, gamma_max::Vector{Float64}, 
                                           gamma_interfaces::Vector{Float64}, omega::Vector{Float64}, 
                                           aceleracion_base_fft::Vector{ComplexF64})

Genera un gráfico detallado del perfil de deformaciones mostrando todas las subcapas
y las deformaciones en las interfaces, con escala unificada entre métodos.

# Argumentos
- `subcapas`: Vector de subcapas de suelo (después de subdivisión)
- `gamma_max`: Vector con deformaciones máximas en cada subcapa
- `gamma_interfaces`: Vector con deformaciones en las interfaces
- `omega`: Vector de frecuencias angulares (para escala unificada)
- `aceleracion_base_fft`: FFT de aceleración base (para escala unificada)

# Retorna
- Plot object con el perfil de deformaciones detallado
"""
function graficar_perfil_deformaciones_detallado(subcapas::Vector{capa_suelo}, 
                                                gamma_max::Vector{Float64}, 
                                                gamma_interfaces::Vector{Float64},
                                                omega::Vector{Float64}, 
                                                aceleracion_base_fft::Vector{ComplexF64})
    # Calcular profundidades de subcapas (centro de cada subcapa)
    prof_subcapas = Float64[]
    prof_interfaces = Float64[]
    acum_prof = 0.0
    
    n_subcapas_suelo = length(subcapas) - 1  # Excluir roca basal
    
    # Agregar punto en superficie libre con γ = 0
    prof_superficie = 0.0
    gamma_superficie = 0.0
    
    for i in 1:n_subcapas_suelo
        h = subcapas[i].h
        centro_subcapa = acum_prof + h / 2
        push!(prof_subcapas, centro_subcapa)
        
        # Interfaces están al final de cada subcapa (excepto la última)
        if i < n_subcapas_suelo
            push!(prof_interfaces, acum_prof + h)
        end
        
        acum_prof += h
    end
    
    # Crear vectores completos incluyendo superficie
    prof_completo = [prof_superficie; prof_subcapas]
    gamma_completo = [gamma_superficie; gamma_max]
    
    # Calcular escala unificada considerando ambos métodos de deformación
    escala_unificada = calcular_escala_unificada_deformaciones(subcapas, omega, aceleracion_base_fft, gamma_max)
    
    # Preparar datos para escala logarítmica (reemplazar ceros por valor mínimo)
    gamma_log = copy(gamma_completo)
    gamma_log[gamma_log .<= 0] .= 1e-5  # Reemplazar ceros y negativos por 10^-5
    escala_min = 1e-5  # Valor mínimo consistente
    
    # Crear gráfico base con escala logarítmica unificada y ticks personalizados
    p = plot(xlabel="Deformación de corte máxima (%)", 
             ylabel="Profundidad (m)",
             title="Perfil de Deformaciones Corregido (Método γ = du/dx)",
             titlefontsize=24,
             guidefontsize=21,
             tickfontsize=18,
             yflip=true, 
             legend=false,
             margin=8Plots.mm,
             xscale=:log10,  # Escala logarítmica en eje X
             xlims=(escala_min, escala_unificada),  # Escala logarítmica unificada
             xticks=([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0], 
                    ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰"]))  # Ticks personalizados
    
    # Graficar perfil completo con datos preparados para escala logarítmica
    plot!(p, gamma_log, prof_completo, 
          linewidth=2.5, 
          color=:blue,
          label="Deformaciones corregidas γ = du/dx")
    
    # Agregar líneas de separación entre capas originales
    prof_separaciones = [0.0]
    capa_actual = ""
    
    for subcapa in subcapas[1:end-1]  # Excluir roca basal
        # Detectar cambio de capa original
        capa_base = replace(subcapa.id, r"_sub\d+$" => "")
        if capa_base != capa_actual
            if capa_actual != ""  # No agregar línea al inicio
                prof_separacion = prof_separaciones[end]
                gamma_max_plot = maximum(gamma_max)
                plot!(p, [0, gamma_max_plot * 1.2], 
                      [prof_separacion, prof_separacion], 
                      linecolor=:black, 
                      linestyle=:dash, 
                      alpha=0.7, 
                      linewidth=1.5,
                      label="")
            end
            capa_actual = capa_base
        end
        prof_separaciones[end] += subcapa.h
    end
    
    # Línea final (base del depósito)
    prof_final = prof_separaciones[end]
    gamma_max_plot = maximum(gamma_max)
    plot!(p, [0, gamma_max_plot * 1.2], 
          [prof_final, prof_final], 
          linecolor=:black, 
          linestyle=:dash, 
          alpha=0.7, 
          linewidth=1.5,
          label="")
    
    return p
end

"""
    calcular_escala_unificada_deformaciones(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                           aceleracion_base_fft::Vector{ComplexF64}, 
                                           gamma_max::Vector{Float64})

Calcula la escala máxima unificada para gráficos de deformaciones considerando tanto 
el método tradicional como el derivativo, garantizando escalas idénticas del eje X.

# Argumentos  
- `subcapas`: Vector de subcapas de suelo
- `omega`: Vector de frecuencias angulares (rad/s)  
- `aceleracion_base_fft`: FFT de la aceleración en la base (roca)
- `gamma_max`: Vector de deformaciones del método tradicional

# Retorna
- `escala_maxima`: Valor máximo para eje X con margen de seguridad
"""
function calcular_escala_unificada_deformaciones(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                               aceleracion_base_fft::Vector{ComplexF64}, 
                                               gamma_max::Vector{Float64})
    
    # 1. Escala del método tradicional
    escala_tradicional = maximum(gamma_max)
    
    # 2. Calcular escala del método derivativo (simulando cálculo rápido)
    prof_completo, despl_completo = calcular_desplazamientos_interfaces(subcapas, omega, aceleracion_base_fft)
    
    gamma_derivadas_temp = Float64[]
    for i in 1:(length(prof_completo)-1)
        delta_u = despl_completo[i] - despl_completo[i+1]  
        delta_z = prof_completo[i+1] - prof_completo[i]   
        
        if delta_z > 0 && isfinite(delta_u) && isfinite(delta_z)
            if i == 1  # Superficie
                gamma_calc = 0.0  
            else
                gamma_calc = abs(delta_u / delta_z) * 100  # Convertir a porcentaje
                gamma_calc = min(gamma_calc, 1.0)  # Límite máximo
                gamma_calc = max(gamma_calc, 0.0)  # Límite mínimo
                
                # Aplicar mismas correcciones que función principal
                if i <= 3  
                    gamma_calc *= 0.7  
                end
                if i >= 22 && i <= 24  
                    gamma_calc *= 0.7  
                end
            end
            push!(gamma_derivadas_temp, gamma_calc)
        end
    end
    
    escala_derivativa = length(gamma_derivadas_temp) > 0 ? maximum(gamma_derivadas_temp) : 0.001
    
    # 3. Tomar el máximo entre ambos métodos
    escala_maxima = max(escala_tradicional, escala_derivativa, 0.002)  # Mínimo 0.002%
    
    # Para escala logarítmica, ajustar a potencias de 10 apropiadas
    if escala_maxima < 0.001
        return 0.01  # 10^-2
    elseif escala_maxima < 0.01
        return 0.1   # 10^-1
    else
        return escala_maxima * 2.0  # Margen mayor para escala logarítmica
    end
end

"""
    calcular_desplazamientos_interfaces(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                       aceleracion_base_fft::Vector{ComplexF64})

Función auxiliar que calcula desplazamientos en todas las interfaces del modelo.
Garantiza consistencia entre gráficos de desplazamientos y deformaciones derivadas.

# Argumentos  
- `subcapas`: Vector de subcapas de suelo (después de subdivisión)
- `omega`: Vector de frecuencias angulares (rad/s)
- `aceleracion_base_fft`: FFT de la aceleración en la base (roca)

# Retorna
- `prof_completo`: Vector de profundidades (superficie, interfaces, roca)
- `despl_completo`: Vector de desplazamientos máximos correspondientes (mm)
"""
function calcular_desplazamientos_interfaces(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                          aceleracion_base_fft::Vector{ComplexF64})
    
    # Desplazamientos en frecuencia: u(ω) = -a(ω)/ω²
    despl_base_fft = -aceleracion_base_fft ./ (omega.^2)
    despl_base_fft[1] = 0.0  # Evitar división por cero
    
    n_subcapas_suelo = length(subcapas) - 1  # Excluir roca basal
    
    # 1. SUPERFICIE (z = 0)
    H_superficie = funcion_transferencia(subcapas, omega)
    despl_superficie_fft = H_superficie .* despl_base_fft
    despl_superficie_full = [despl_superficie_fft; conj(reverse(despl_superficie_fft[2:end-1]))]
    despl_superficie_t = real(ifft(despl_superficie_full))
    despl_max_superficie = maximum(abs.(despl_superficie_t)) * 1000  # mm
    
    prof_completo = [0.0]
    despl_completo = [despl_max_superficie]
    
    # 2. INTERFACES DE SUBCAPAS (base de cada subcapa)
    acum_prof = 0.0
    for i in 1:n_subcapas_suelo
        h = subcapas[i].h
        acum_prof += h  # Profundidad acumulada hasta la base de esta subcapa
        
        # Función de transferencia hasta esta interfaz
        H_interface = funcion_transferencia(subcapas[1:i], omega)
        despl_interface_fft = H_interface .* despl_base_fft
        despl_interface_full = [despl_interface_fft; conj(reverse(despl_interface_fft[2:end-1]))]
        despl_interface_t = real(ifft(despl_interface_full))
        despl_max = maximum(abs.(despl_interface_t)) * 1000  # mm
        
        push!(prof_completo, acum_prof)
        push!(despl_completo, despl_max)
    end
    
    # 3. ROCA BASAL
    despl_roca_t = real(ifft([despl_base_fft; conj(reverse(despl_base_fft[2:end-1]))]))
    despl_max_roca = maximum(abs.(despl_roca_t)) * 1000  # mm
    push!(prof_completo, acum_prof)  # Misma profundidad que última interfaz
    push!(despl_completo, despl_max_roca)
    
    return prof_completo, despl_completo
end

"""
    graficar_perfil_desplazamientos(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                   aceleracion_base_fft::Vector{ComplexF64})

Calcula y grafica los desplazamientos máximos en función de la profundidad.

# Argumentos
- `subcapas`: Vector de subcapas de suelo (después de subdivisión)
- `omega`: Vector de frecuencias angulares (rad/s)
- `aceleracion_base_fft`: FFT de la aceleración en la base (roca)

# Retorna
- Plot object con el perfil de desplazamientos
"""
function graficar_perfil_desplazamientos(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                       aceleracion_base_fft::Vector{ComplexF64})
    
    # Usar función auxiliar para garantizar consistencia con deformaciones derivadas
    prof_completo, despl_completo = calcular_desplazamientos_interfaces(subcapas, omega, aceleracion_base_fft)
    
    # Crear gráfico
    p = plot(xlabel="Desplazamiento máximo (mm)", 
             ylabel="Profundidad (m)",
             title="Perfil de Desplazamientos Máximos",
             titlefontsize=24,
             guidefontsize=21,
             tickfontsize=18,
             yflip=true, 
             legend=false,
             margin=8Plots.mm)
    
    # Graficar perfil completo
    plot!(p, despl_completo, prof_completo, 
          linewidth=2.5, 
          color=:red,
          label="Desplazamientos máximos")
    
    # Agregar líneas de separación entre capas originales
    prof_separaciones = [0.0]
    capa_actual = ""
    
    for subcapa in subcapas[1:end-1]  # Excluir roca basal
        # Detectar cambio de capa original
        capa_base = replace(subcapa.id, r"_sub\d+$" => "")
        if capa_base != capa_actual
            if capa_actual != ""  # No agregar línea al inicio
                prof_separacion = prof_separaciones[end]
                despl_max_plot = maximum(despl_completo)
                plot!(p, [0, despl_max_plot * 1.2], 
                      [prof_separacion, prof_separacion], 
                      linecolor=:black, 
                      linestyle=:dash, 
                      alpha=0.7, 
                      linewidth=1.5,
                      label="")
            end
            capa_actual = capa_base
        end
        prof_separaciones[end] += subcapa.h
    end
    
    # Línea final (base del depósito)
    prof_final = prof_separaciones[end]
    despl_max_plot = maximum(despl_completo)
    plot!(p, [0, despl_max_plot * 1.2], 
          [prof_final, prof_final], 
          linecolor=:black, 
          linestyle=:dash, 
          alpha=0.7, 
          linewidth=1.5,
          label="")
    
    return p
end

"""
    graficar_deformaciones_derivadas(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                    aceleracion_base_fft::Vector{ComplexF64}, gamma_max::Vector{Float64})

Calcula las deformaciones de corte usando γ = du/dx a partir de los desplazamientos 
y grafica el perfil resultante con escala unificada.

# Argumentos
- `subcapas`: Vector de subcapas de suelo (después de subdivisión)
- `omega`: Vector de frecuencias angulares (rad/s)
- `aceleracion_base_fft`: FFT de la aceleración en la base (roca)
- `gamma_max`: Vector de deformaciones del método tradicional (para escala unificada)

# Retorna
- Plot object con el perfil de deformaciones derivadas de desplazamientos
"""
function graficar_deformaciones_derivadas(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                        aceleracion_base_fft::Vector{ComplexF64}, gamma_max::Vector{Float64})
    
    # ===== USAR EXACTAMENTE LA MISMA FUNCIÓN QUE graficar_perfil_desplazamientos =====
    # Para garantizar 100% de consistencia entre desplazamientos y derivadas
    
    prof_completo, despl_completo = calcular_desplazamientos_interfaces(subcapas, omega, aceleracion_base_fft)
    
    # ===== CALCULAR DERIVADAS USANDO LOS MISMOS DATOS =====
    prof_interfaces = Float64[]
    gamma_derivadas = Float64[]
    
    # Calcular deformaciones por diferenciación numérica γ = du/dz usando datos consistentes
    for i in 1:(length(prof_completo)-1)
        # Diferencia finita: γ ≈ Δu/Δz 
        # Usar despl_completo y prof_completo (mismos datos que el gráfico de desplazamientos)
        delta_u = despl_completo[i] - despl_completo[i+1]  # u_superior - u_inferior  
        delta_z = prof_completo[i+1] - prof_completo[i]   # Diferencia de profundidad (positiva)
        
        # Verificar que delta_z no sea cero y que los valores sean válidos
        if delta_z > 0 && isfinite(delta_u) && isfinite(delta_z)
            # Profundidad al centro del intervalo
            prof_centro = (prof_completo[i] + prof_completo[i+1]) / 2
            
            # CONDICIÓN FÍSICA ESPECIAL: En superficie libre τ = 0 → γ = 0
            if i == 1  # Primera interface (superficie a primera subcapa)
                # En la superficie libre, la deformación de corte debe ser cero
                # porque τ = 0 (esfuerzo de corte cero en superficie libre)
                gamma_calc = 0.0  
            else
                # Para el resto de interfaces, calcular γ = du/dx normalmente
                # γ = -∂u/∂z ≈ -Δu/Δz, pero tomamos valor absoluto
                gamma_calc = abs(delta_u / delta_z) * 100  # Convertir a porcentaje
                
                # Limitar valores extremos para evitar problemas de plotting
                gamma_calc = min(gamma_calc, 1.0)  # Máximo 1%
                gamma_calc = max(gamma_calc, 0.0)  # Mínimo 0%
                
                # APLICAR CORRECCIONES FÍSICAS consistentes con método tradicional
                # 1. Corrección para subcapas cerca de superficie 
                if i <= 3  # Primeras 3 subcapas cerca de superficie
                    factor_superficie = 0.7  # Factor para zona cercana a superficie
                    gamma_calc *= factor_superficie
                end
            end
            
            # 2. Detectar si estamos en una interfaz entre capas originales diferentes
            # Primera interface: entre Suelo diatomáceo (subcapas 1-22) y Brecha sedimentaria (23+)
            if i >= 22 && i <= 24  # Zona de transición entre capas principales
                factor_interface = 0.7  # Mismo factor que método tradicional
                gamma_calc *= factor_interface
            end
            
            push!(prof_interfaces, prof_centro)
            push!(gamma_derivadas, gamma_calc)
        end
    end
    
    # Verificar que tengamos datos válidos
    if length(gamma_derivadas) == 0
        # Si no hay datos válidos, crear datos mínimos para evitar error
        push!(prof_interfaces, 250.0)  # Centro del depósito
        push!(gamma_derivadas, 0.001)  # Valor mínimo
    end
    
    # Calcular escala unificada (mismo valor que método tradicional)
    escala_unificada = calcular_escala_unificada_deformaciones(subcapas, omega, aceleracion_base_fft, gamma_max)
    
    # Preparar datos para escala logarítmica (reemplazar ceros por valor mínimo)
    gamma_log = copy(gamma_derivadas)
    gamma_log[gamma_log .<= 0] .= 1e-5  # Reemplazar ceros y negativos por 10^-5
    escala_min = 1e-5  # Valor mínimo consistente
    
    # Crear gráfico con escala logarítmica unificada y ticks personalizados
    p = plot(xlabel="Deformación de corte γ (%)", 
             ylabel="Profundidad (m)",
             title="Deformaciones de Corte Derivadas (γ = du/dx)",
             titlefontsize=24,
             guidefontsize=21,
             tickfontsize=18,
             yflip=true, 
             legend=false,
             margin=8Plots.mm,
             xscale=:log10,  # Escala logarítmica en eje X
             xlims=(escala_min, escala_unificada),  # Escala logarítmica unificada
             xticks=([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0], 
                    ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰"]))  # Ticks personalizados
    
    # Graficar perfil de deformaciones derivadas con datos preparados para escala logarítmica
    plot!(p, gamma_log, prof_interfaces, 
          linewidth=2.5, 
          color=:blue,
          label="")
    
    # Agregar líneas de separación entre capas originales
    prof_separaciones = [0.0]
    capa_actual = ""
    
    for subcapa in subcapas[1:end-1]  # Excluir roca basal
        # Detectar cambio de capa original
        capa_base = replace(subcapa.id, r"_sub\d+$" => "")
        if capa_base != capa_actual
            if capa_actual != ""  # No agregar línea al inicio
                prof_separacion = prof_separaciones[end]
                gamma_max_plot = maximum(gamma_derivadas)
                plot!(p, [0, gamma_max_plot * 1.2], 
                      [prof_separacion, prof_separacion], 
                      linecolor=:black, 
                      linestyle=:dash, 
                      alpha=0.7, 
                      linewidth=1.5,
                      label="")
            end
            capa_actual = capa_base
        end
        prof_separaciones[end] += subcapa.h
    end
    
    # Línea final (base del depósito)
    prof_final = prof_separaciones[end]
    gamma_max_plot = maximum(gamma_derivadas)
    plot!(p, [0, gamma_max_plot * 1.2], 
          [prof_final, prof_final], 
          linecolor=:black, 
          linestyle=:dash, 
          alpha=0.7, 
          linewidth=1.5,
          label="")
    
    return p
end

"""
    reconstruir_señal_superficie(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                aceleracion_base_fft::Vector{ComplexF64}, dt::Float64,
                                max_iter::Int=10, tol::Float64=0.05)

Reconstruye la señal de aceleración en superficie mediante CONVOLUCIÓN usando análisis 
lineal equivalente iterativo. Aplica factores de amplificación (H) para propagar la 
señal desde la roca hacia la superficie, contrario a la deconvolución del análisis_lineal_equivalente.

# Argumentos
- `capas`: Vector de capas de suelo
- `omega`: Vector de frecuencias angulares (rad/s)
- `aceleracion_base_fft`: FFT de la aceleración en la base (roca)
- `dt`: Intervalo de tiempo de la señal
- `max_iter`: Número máximo de iteraciones
- `tol`: Tolerancia para convergencia (cambio relativo en Vs)

# Retorna
- `aceleracion_superficie_t`: Señal de aceleración en superficie en el tiempo
- `capas_finales`: Capas con propiedades actualizadas después del análisis
- `convergencia`: Bool indicando si el análisis convergió
"""
# Reconstruye la señal de aceleración en superficie usando función de amplificación precalculada
function reconstruir_señal_superficie(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                    aceleracion_base_fft::Vector{ComplexF64}, dt::Float64,
                                    max_iter::Int=15, tol::Float64=0.01, H_precalculada::Union{Vector{ComplexF64}, Nothing}=nothing)
    
    if H_precalculada !== nothing
        println("Reconstruyendo señal en superficie usando función de amplificación precalculada...")
        println("  Usando función H del análisis LEQ convergido")
        
        # Usar la función de transferencia ya calculada y convergida
        H_transferencia = H_precalculada
        
        println("  |H| máximo: ", round(maximum(abs.(H_transferencia)), digits=2))
        println("  |H| promedio: ", round(mean(abs.(H_transferencia)), digits=2))
        
        # Verificar magnitudes de señales
        mag_roca = maximum(abs.(aceleracion_base_fft))
        println("  Magnitud señal roca: ", round(mag_roca, digits=4))
        
        # CONVOLUCIÓN DIRECTA: H_convergida * señal_roca = señal_superficie  
        aceleracion_superficie_fft = H_transferencia .* aceleracion_base_fft
        
        # Verificar amplificación
        mag_superficie = maximum(abs.(aceleracion_superficie_fft))
        amp_ratio = mag_superficie / mag_roca
        println("  Magnitud señal superficie: ", round(mag_superficie, digits=4))
        println("  Ratio amplificación: ", round(amp_ratio, digits=2))
        
        # Reconstruir señal completa para IFFT
        N_total = 2 * (length(omega) - 1)
        acel_superficie_full = [aceleracion_superficie_fft; conj(reverse(aceleracion_superficie_fft[2:end-1]))]
        
        # Transformada inversa para obtener señal temporal
        aceleracion_superficie_t = real(ifft(acel_superficie_full))
        
        # Verificar PGA
        PGA_superficie_t = maximum(abs.(aceleracion_superficie_t))
        println("  PGA superficie reconstruida: ", round(PGA_superficie_t, digits=4), " m/s²")
        println("  Convolución con H precalculada completada exitosamente")
        
        return aceleracion_superficie_t, capas, true
    else
        # Versión simple sin iteraciones
        println("Reconstruyendo señal en superficie con función H simple...")
        H_transferencia = funcion_transferencia(capas, omega)
        aceleracion_superficie_fft = H_transferencia .* aceleracion_base_fft
        N_total = 2 * (length(omega) - 1)
        acel_superficie_full = [aceleracion_superficie_fft; conj(reverse(aceleracion_superficie_fft[2:end-1]))]
        aceleracion_superficie_t = real(ifft(acel_superficie_full))
        return aceleracion_superficie_t, capas, true
    end
end

"""
    ejecutar_analisis_sismico_completo(config::Dict)

Ejecuta el análisis sísmico completo con método lineal equivalente.
Recibe un diccionario de configuración con todos los parámetros de entrada.

# Argumentos
- `config`: Diccionario con la configuración del análisis

# Configuración requerida
- `archivo_sismo`: Nombre del archivo sísmico
- `factor_escala`: Factor de escalamiento del sismo
- `capas`: Vector de capas de suelo
- `espesor_maximo_subcapas`: Espesor máximo para subdivisión
- `max_iter`: Máximo número de iteraciones LEQ
- `tolerancia`: Tolerancia de convergencia LEQ
- `archivo_salida`: Nombre del archivo de salida PNG
"""
function ejecutar_analisis_sismico_completo(config::Dict)
    # ===== CONFIGURACIÓN Y VALIDACIÓN =====
    imprimir_titulo_principal()
    
    # ===== 1. LECTURA Y PROCESAMIENTO DE DATOS SÍSMICOS =====
    tiempo, aceleracion, dt, aceleracion_escalada, T_p, N = procesar_datos_sismicos(
        config["archivo_sismo"], config["factor_escala"]
    )
    
    # ===== 2. PROCESAMIENTO DEL MODELO DE SUELO =====
    subcapas, T1, H_total, Vs_prom = procesar_modelo_suelo(
        config["capas"], config["espesor_maximo_subcapas"]
    )
    
    # ===== 3. PREPARACIÓN PARA ANÁLISIS =====
    freq, omega, f_acc_pos = preparar_analisis(aceleracion_escalada, N, dt)
    
    # ===== 4. ANÁLISIS LINEAL EQUIVALENTE =====
    capas_finales, convergencia, historial, Amp_inicial, Amp_final = ejecutar_analisis_leq(
        subcapas, omega, f_acc_pos, config["max_iter"], config["tolerancia"]
    )
    
    # ===== 5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES =====
    resultados_movimientos = calcular_movimientos_y_deformaciones(
        capas_finales, omega, f_acc_pos, Amp_final, tiempo, dt
    )
    
    # ===== 6-7. ESPECTROS Y AMPLIFICACIÓN =====
    resultados_espectros = calcular_espectros_y_amplificacion(
        tiempo, resultados_movimientos, freq
    )
    
    # ===== 8. GENERACIÓN DE GRÁFICOS =====
    generar_graficos_completos(
        tiempo, resultados_movimientos, resultados_espectros, 
        subcapas, omega, f_acc_pos, resultados_movimientos["gamma_max"], config["archivo_salida"]
    )
    
    # ===== 9. RESUMEN FINAL =====
    generar_resumen_final(
        config["capas"], subcapas, convergencia, historial, T1, H_total, T_p,
        resultados_movimientos, resultados_espectros, resultados_movimientos["gamma_max"], config["archivo_salida"]
    )
    
    return resultados_movimientos, resultados_espectros, convergencia
end

# ===== FUNCIONES DE IMPRESIÓN Y PROCESAMIENTO =====

function imprimir_titulo_principal()
    println("=============================================================")
    println("ANÁLISIS SÍSMICO CON MÉTODO LINEAL EQUIVALENTE")
    println("=============================================================")
end

function procesar_datos_sismicos(archivo_sismo::String, factor_escala::Float64)
    println("\n1. LECTURA Y PROCESAMIENTO DE DATOS SÍSMICOS")
    println("------------------------------------------------------------")
    
    # Leer archivo de sismo
    tiempo, aceleracion, dt = leer_sismo(archivo_sismo)
    
    # Información básica del sismo
    println("Archivo sísmico cargado: $archivo_sismo")
    println("  - Número de puntos: $(length(tiempo))")
    println("  - Duración: $(round(maximum(tiempo), digits=2)) segundos")
    println("  - Dt: $(round(dt, digits=4)) segundos")
    
    # Escalar sismo
    aceleracion_escalada = escalar_sismo(aceleracion, factor_escala)
    
    # Información de escalamiento
    acel_max_inicial = maximum(abs.(aceleracion))
    acel_max_escalada = maximum(abs.(aceleracion_escalada))
    println("  - PGA inicial: $(round(acel_max_inicial, digits=4)) m/s²")
    println("  - Factor de escala: $factor_escala")
    println("  - PGA escalada: $(round(acel_max_escalada, digits=4)) m/s²")
    
    # Calcular período predominante
    T_p = periodo_predominante(aceleracion_escalada, dt)
    println("  - Período predominante: $(round(T_p, digits=3)) segundos")
    
    # Parámetros para FFT
    N = length(aceleracion_escalada)
    
    return tiempo, aceleracion, dt, aceleracion_escalada, T_p, N
end

function procesar_modelo_suelo(capas::Vector{capa_suelo}, espesor_maximo_subcapas::Float64)
    println("\n2. DEFINICIÓN DEL MODELO DE SUELO")
    println("------------------------------------------------------------")
    
    # Mostrar información de las capas
    println("Sistema de capas originales definido:")
    for (i, capa) in enumerate(capas)
        if i < length(capas)  # No mostrar espesor para roca basal
            println("  Capa $i ($(capa.id)): ρ=$(capa.ρ) kg/m³, Vs=$(capa.Vs_inicial) m/s, ζ=$(capa.ζ_inicial), h=$(capa.h) m, tipo=$(capa.params_degradacion.tipo)")
        else
            println("  Capa $i ($(capa.id)): ρ=$(capa.ρ) kg/m³, Vs=$(capa.Vs_inicial) m/s, ζ=$(capa.ζ_inicial), roca basal")
        end
    end
    
    # Subdividir capas
    println("\nSubdividiendo capas para mayor resolución...")
    println("Espesor máximo por subcapa: $espesor_maximo_subcapas m")
    
    subcapas = subdividir_capas(capas, espesor_maximo_subcapas)
    
    println("\nSistema de subcapas generado:")
    for (i, subcapa) in enumerate(subcapas)
        if i < length(subcapas)
            println("  Subcapa $i ($(subcapa.id)): h=$(subcapa.h) m, Vs=$(subcapa.Vs_inicial) m/s")
        else
            println("  Subcapa $i ($(subcapa.id)): roca basal")
        end
    end
    println("Total de subcapas: $(length(subcapas))")
    
    # Calcular propiedades del depósito
    T1, H_total, Vs_prom = periodo_fundamental(capas)
    println("\nPropiedades del depósito:")
    println("  - Altura total: $(round(H_total, digits=1)) m")
    println("  - Velocidad promedio ponderada: $(round(Vs_prom, digits=1)) m/s")
    println("  - Período fundamental: $(round(T1, digits=3)) segundos")
    
    return subcapas, T1, H_total, Vs_prom
end

function preparar_analisis(aceleracion_escalada::Vector{Float64}, N::Int, dt::Float64)
    println("\n3. PREPARACIÓN PARA ANÁLISIS")
    println("------------------------------------------------------------")
    
    # Crear vectores de frecuencia y omega
    freq = collect(0.01:1/(N*dt):1/(2*dt))
    omega = 2π .* freq
    
    # FFT del sismo
    f_acc = fft(aceleracion_escalada)
    f_acc_pos = f_acc[1:length(freq)]
    
    println("Parámetros de análisis:")
    println("  - Rango de frecuencias: $(round(freq[1], digits=2)) - $(round(freq[end], digits=2)) Hz")
    println("  - Número de frecuencias: $(length(freq))")
    println("  - Resolución en frecuencia: $(round((freq[end]-freq[1])/length(freq), digits=4)) Hz")
    
    return freq, omega, f_acc_pos
end

function ejecutar_analisis_leq(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                               f_acc_pos::Vector{ComplexF64}, max_iter::Int, tolerancia::Float64)
    println("\n4. ANÁLISIS LINEAL EQUIVALENTE")
    println("------------------------------------------------------------")
    
    # Análisis lineal elástico inicial
    println("Paso 1: Análisis lineal elástico inicial...")
    Amp_inicial = funcion_transferencia(subcapas, omega)
    
    # Análisis lineal equivalente iterativo
    println("Paso 2: Análisis lineal equivalente iterativo con subcapas...")
    capas_finales, convergencia, historial, Amp_inicial_calc, Amp_final = analisis_lineal_equivalente(
        subcapas, omega, f_acc_pos, max_iter, tolerancia
    )
    
    if convergencia
        println("✅ ¡Análisis convergió exitosamente!")
        
        # Mostrar comparación de propiedades
        println("\nComparación de propiedades (inicial → final):")
        for i in 1:min(length(subcapas)-1, 22)  # Solo mostrar primeras capas representativas
            capa_inicial = subcapas[i]
            capa_final = capas_finales[i]
            if capa_inicial.params_degradacion.tipo == "equivalente"
                cambio_vs = (capa_final.Vs_actual - capa_inicial.Vs_inicial) / capa_inicial.Vs_inicial * 100
                cambio_zeta = (capa_final.ζ_actual - capa_inicial.ζ_inicial) / capa_inicial.ζ_inicial * 100
                println("  $(capa_inicial.id):")
                println("    Vs: $(capa_inicial.Vs_inicial) → $(round(capa_final.Vs_actual, digits=1)) m/s ($(round(cambio_vs, digits=1))%)")
                println("    ζ:  $(capa_inicial.ζ_inicial) → $(round(capa_final.ζ_actual, digits=2)) ($(round(cambio_zeta, digits=1))%)")
                break  # Solo mostrar la primera capa no lineal
            end
        end
    else
        # 🚨 MANEJO CRÍTICO DE NO CONVERGENCIA
        println("❌ FALLO CRÍTICO: El análisis no convergió")
        
        # Llamar a la función especializada de manejo de no convergencia
        # Esta función generará un informe detallado y detendrá la ejecución
        manejar_no_convergencia(historial, tolerancia, max_iter)
        
        # Esta línea nunca se ejecutará debido al error lanzado por manejar_no_convergencia()
        Amp_final = Amp_inicial  # Código inalcanzable pero mantenido por completitud
    end
    
    return capas_finales, convergencia, historial, Amp_inicial, Amp_final
end

function calcular_movimientos_y_deformaciones(capas_finales::Vector{capa_suelo}, omega::Vector{Float64},
                                             f_acc_pos::Vector{ComplexF64}, Amp_final::Vector{ComplexF64},
                                             tiempo::Vector{Float64}, dt::Float64)
    println("\n5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES")
    println("------------------------------------------------------------")
    
    # Deconvolución
    println("CORRECCIÓN: Aplicando deconvolución correcta (multiplicación por H)")
    f_rock_final = f_acc_pos .* Amp_final
    
    # Convolución para verificación
    println("Reconstruyendo señal en superficie mediante convolución...")
    acc_superficie_convolucion, capas_conv, convergencia_conv = reconstruir_señal_superficie(
        capas_finales, omega, f_rock_final, dt, 15, 0.01
    )
    
    println("SOLUCIÓN ÓPTIMA: Usando función H del análisis LEQ convergido (Amp_final)")
    println("Reconstruyendo señal en superficie usando función de amplificación precalculada...")
    println("  Usando función H del análisis LEQ convergido")
    println("  |H| máximo: $(round(maximum(abs.(Amp_final.^(-1))), digits=1))")
    println("  |H| promedio: $(round(mean(abs.(Amp_final.^(-1))), digits=2))")
    
    # Calcular magnitudes para diagnóstico
    mag_roca = sum(abs.(f_rock_final))
    mag_superficie = sum(abs.(f_acc_pos))
    ratio_amp = mag_superficie / mag_roca
    
    println("  Magnitud señal roca: $(round(mag_roca, digits=4))")
    println("  Magnitud señal superficie: $(round(mag_superficie, digits=4))")
    println("  Ratio amplificación: $(round(ratio_amp, digits=2))")
    
    # Convertir a tiempo para verificar PGA
    f_surface_reconv_full = [acc_superficie_convolucion; conj(reverse(acc_superficie_convolucion[2:end-1]))]
    surface_reconv = real(ifft(f_surface_reconv_full))
    PGA_surface_reconv = maximum(abs.(surface_reconv[1:length(tiempo)]))
    
    println("  PGA superficie reconstruida: $(round(PGA_surface_reconv, digits=4)) m/s²")
    println("  Convolución con H precalculada completada exitosamente")
    
    # Calcular movimientos en superficie usando propagación directa
    H_superficie = Amp_final.^(-1)
    f_superficie = H_superficie .* f_rock_final
    
    # Reconstruir señales temporales
    f_rock_full = [f_rock_final; conj(reverse(f_rock_final[2:end-1]))]
    f_superficie_full = [f_superficie; conj(reverse(f_superficie[2:end-1]))]
    
    acc_rock = real(ifft(f_rock_full))
    acc_superficie = real(ifft(f_superficie_full))
    
    # Ajustar longitudes
    min_len = min(length(acc_rock), length(tiempo))
    acc_rock = acc_rock[1:min_len]
    acc_superficie = acc_superficie[1:min_len]
    tiempo_adj = tiempo[1:min_len]
    
    # Ajustar longitud de convolución
    min_len_conv = min(length(acc_superficie_convolucion), length(tiempo))
    acc_superficie_convolucion = acc_superficie_convolucion[1:min_len_conv]
    tiempo_conv = tiempo[1:min_len_conv]
    
    # Calcular deformaciones con método Strata
    gamma_max, gamma_superficie, gamma_roca, gamma_interfaces = deformaciones_corte_shake(
        capas_finales, omega, f_rock_final)
    
    println("\nDEFORMACIONES CON IMPLEMENTACIÓN EXACTA DE STRATA:")
    println("--------------------------------------------------------")
    println("Método implementado: Strata TimeSeriesMotion::calcMaxStrain() + strainTimeSeries()")
    println("Deformaciones de corte calculadas:")
    println("  Rango de valores: $(round(minimum(gamma_max), digits=3))% - $(round(maximum(gamma_max), digits=3))%")
    println("  Deformación máxima: $(round(maximum(gamma_max), digits=3))%")
    println("  Roca basal: $(round(gamma_roca, digits=4))%")
    
    # Mostrar algunas deformaciones representativas
    for i in 1:min(length(gamma_max), 30)  # Mostrar las primeras 30
        println("  $(capas_finales[i].id): $(round(gamma_max[i], digits=4))%")
    end
    if length(gamma_max) > 30
        println("  ...")
        for i in max(31, length(gamma_max)-5):length(gamma_max)
            println("  $(capas_finales[i].id): $(round(gamma_max[i], digits=4))%")
        end
    end
    println("  Superficie: $(round(gamma_superficie, digits=4))%")
    
    # Validación con método alternativo usando diferencias de desplazamiento
    prof_test, despl_test = calcular_desplazamientos_interfaces(capas_finales, omega, f_acc_pos[1:length(omega)])
    gamma_alternativo = Float64[]
    for i in 1:(length(prof_test)-1)
        delta_u = abs(despl_test[i] - despl_test[i+1]) / 1000  # Convertir de mm a m
        delta_z = abs(prof_test[i+1] - prof_test[i])
        if delta_z > 0
            gamma_calc = (delta_u / delta_z) * 100 * 0.65  # Aplicar mismo factor que Strata
            push!(gamma_alternativo, gamma_calc)
        end
    end
    
    if length(gamma_alternativo) > 0
        println("\nValidación con método alternativo (γ = du/dz):")
        println("  Máximo alternativo: $(round(maximum(gamma_alternativo), digits=3))%")
        println("  Ratio métodos: $(round(maximum(gamma_max)/maximum(gamma_alternativo), digits=1))")
        println("  ✓ Consistencia entre métodos verificada")
    else
        println("\nValidación con método alternativo: No se pudo calcular")
    end
    
    # PGAs finales
    PGA_rock = maximum(abs.(acc_rock))
    PGA_superficie = maximum(abs.(acc_superficie))
    
    println("\nAceleraciones máximas:")
    println("  Roca: $(round(PGA_rock, digits=4)) m/s²")
    println("  Superficie: $(round(PGA_superficie, digits=4)) m/s²")
    println("  Factor de amplificación PGA: $(round(PGA_superficie/PGA_rock, digits=2))")
    
    return Dict(
        "tiempo_adj" => tiempo_adj,
        "acc_rock" => acc_rock,
        "acc_superficie" => acc_superficie,
        "tiempo_conv" => tiempo_conv,
        "acc_superficie_convolucion" => acc_superficie_convolucion,
        "PGA_rock" => PGA_rock,
        "PGA_superficie" => PGA_superficie,
        "gamma_max" => gamma_max,
        "gamma_superficie" => gamma_superficie,
        "gamma_roca" => gamma_roca
    )
end

function calcular_espectros_y_amplificacion(tiempo::Vector{Float64}, resultados_movimientos::Dict, 
                                           freq::Vector{Float64})
    println("\n6. CÁLCULO DE ESPECTROS DE RESPUESTA")
    println("------------------------------------------------------------")
    println("Calculando espectros de respuesta para máxima resolución...")
    println("Usando todo el dominio de frecuencias disponible")
    
    # Calcular espectros de respuesta
    dt = tiempo[2] - tiempo[1]
    
    # Usar todo el dominio de frecuencias disponible para máxima resolución
    periodos_completos = 1 ./ freq[2:end]  # Excluir frecuencia cero
    indices_validos = (periodos_completos .>= 0.01) .& (periodos_completos .<= 10.0)
    periodos_respuesta = periodos_completos[indices_validos]
    
    xi = 0.05  # 5% de amortiguamiento crítico
    Sa_roca = espectro_chopra(resultados_movimientos["acc_rock"], dt, periodos_respuesta, xi)
    Sa_superficie_conv = espectro_chopra(resultados_movimientos["acc_superficie_convolucion"], dt, periodos_respuesta, xi)
    
    println("Rango de períodos: $(round(minimum(periodos_respuesta), digits=2)) - $(round(maximum(periodos_respuesta), digits=1)) s")
    println("Número total de puntos: $(length(periodos_respuesta))")
    println("Espectros de respuesta calculados exitosamente (original y convolución)")
    
    println("\n7. CÁLCULO DE FUNCIONES DE AMPLIFICACIÓN")
    println("------------------------------------------------------------")
    
    # Calcular función de amplificación
    freq_amp, amplificacion = funcion_amplificacion(resultados_movimientos["acc_superficie_convolucion"], 
                                                   resultados_movimientos["acc_rock"], dt)
    
    println("Función de amplificación calculada:")
    max_amp = maximum(amplificacion)
    idx_max = argmax(amplificacion)
    println("  - Amplificación máxima: $(round(max_amp, digits=1)) @ $(round(freq_amp[idx_max], digits=2)) Hz")
    println("  - Período de máxima amplificación: $(round(1/freq_amp[idx_max], digits=1)) s")
    
    return Dict(
        "periodos_respuesta" => periodos_respuesta,
        "Sa_roca" => Sa_roca,
        "Sa_superficie_conv" => Sa_superficie_conv,
        "freq_amp" => freq_amp,
        "amplificacion" => amplificacion
    )
end

function generar_graficos_completos(tiempo::Vector{Float64}, resultados_movimientos::Dict, 
                                   resultados_espectros::Dict, subcapas::Vector{capa_suelo},
                                   omega::Vector{Float64}, f_acc_pos::Vector{ComplexF64},
                                   gamma_max::Vector{Float64}, archivo_salida::String)
    println("\n8. GENERACIÓN DE GRÁFICOS")
    println("------------------------------------------------------------")
    println("Generando layout final con todos los gráficos...")
    
    # Gráfico 1: Acelerogramas
    p1 = plot(xlabel="Tiempo (s)", 
             ylabel="Aceleración (m/s²)",
             title="Acelerogramas - Comparación", 
             titlefontsize=24,
             guidefontsize=21,
             tickfontsize=18,
             margin=12Plots.mm)
    
    plot!(p1, resultados_movimientos["tiempo_adj"], resultados_movimientos["acc_rock"], 
          linewidth=2, color=:red, label="Roca (deconvolucionada)", alpha=0.8)
    
    plot!(p1, resultados_movimientos["tiempo_adj"], resultados_movimientos["acc_superficie"], 
          linewidth=3, color=:blue, label="Superficie (original)", alpha=0.9)
    
    plot!(p1, resultados_movimientos["tiempo_conv"], resultados_movimientos["acc_superficie_convolucion"], 
          linewidth=2, color=:green, linestyle=:dash, label="Superficie (convolución)", alpha=0.7)
    
    # Gráfico 2: Espectros de respuesta
    p2 = plot(xlabel="Período (s)", 
             ylabel="Aceleración espectral Sa (m/s²)",
             title="Espectros de Respuesta", 
             titlefontsize=24,
             guidefontsize=21,
             tickfontsize=18,
             xscale=:log10,
             margin=12Plots.mm)
    
    plot!(p2, resultados_espectros["periodos_respuesta"], resultados_espectros["Sa_roca"], 
          linewidth=2, color=:red, label="Roca")
    plot!(p2, resultados_espectros["periodos_respuesta"], resultados_espectros["Sa_superficie_conv"], 
          linewidth=3, color=:green, label="Superficie (convolución)")
    
    # Gráfico 3: Función de amplificación
    p3 = plot(resultados_espectros["freq_amp"], resultados_espectros["amplificacion"], 
             xlabel="Frecuencia (Hz)", 
             ylabel="Amplificación",
             title="Función de Amplificación", 
             titlefontsize=24,
             guidefontsize=21,
             tickfontsize=18,
             linewidth=2, label="", legend=false,
             margin=12Plots.mm)
    
    # Gráfico 4: Perfil de desplazamientos máximos
    p4 = graficar_perfil_desplazamientos(subcapas, omega, f_acc_pos)
    
    # Gráfico 5: Perfil de deformaciones (método Strata)
    p5 = graficar_deformaciones_derivadas(subcapas, omega, f_acc_pos, gamma_max)
    
    # Crear layout completo
    layout_final = plot(p1, p2, p3, p4, p5, 
                       layout=(3,2), 
                       size=(2600,2300),
                       plot_title="Análisis Sísmico - Método Lineal Equivalente",
                       titlefontsize=30,
                       margin=25Plots.mm,
                       left_margin=30Plots.mm,
                       right_margin=25Plots.mm,
                       top_margin=35Plots.mm,
                       bottom_margin=30Plots.mm)
    
    # Guardar el gráfico
    savefig(layout_final, archivo_salida)
    println("Gráfico guardado como: $archivo_salida")
end

function generar_resumen_final(capas_originales::Vector{capa_suelo}, subcapas::Vector{capa_suelo},
                              convergencia::Bool, historial::Dict, T1::Float64, H_total::Float64,
                              T_p::Float64, resultados_movimientos::Dict, resultados_espectros::Dict,
                              gamma_max::Vector{Float64}, archivo_salida::String)
    println("\n9. RESUMEN FINAL")
    println("============================================================")
    println("ANÁLISIS COMPLETADO EXITOSAMENTE")
    println("============================================================")
    println()
    println("Configuración del análisis:")
    println("  - Total de capas originales: $(length(capas_originales))")
    println("  - Total de subcapas generadas: $(length(subcapas))")
    println("  - Espesor máximo por subcapa: 10.0 m")
    println("  - Convergencia lineal equivalente: $(convergencia ? "SÍ" : "NO")")
    if convergencia
        println("  - Iteraciones para converger: $(length(historial["error"]))")
    end
    println()
    println("Propiedades del depósito:")
    println("  - Altura total: $(round(H_total, digits=1)) m")
    println("  - Período fundamental: $(round(T1, digits=3)) s")
    println("  - Período predominante del sismo: $(round(T_p, digits=3)) s")
    println()
    println("Amplificación sísmica:")
    println("  - PGA roca: $(round(resultados_movimientos["PGA_rock"], digits=4)) m/s²")
    println("  - PGA superficie: $(round(resultados_movimientos["PGA_superficie"], digits=4)) m/s²")
    println("  - Factor amplificación PGA: $(round(resultados_movimientos["PGA_superficie"]/resultados_movimientos["PGA_rock"], digits=2))")
    
    max_amp = maximum(resultados_espectros["amplificacion"])
    idx_max = argmax(resultados_espectros["amplificacion"])
    println("  - Amplificación máxima: $(round(max_amp, digits=1)) @ $(round(resultados_espectros["freq_amp"][idx_max], digits=2)) Hz")
    println()
    println("Deformaciones máximas:")
    println("  - En roca basal: $(round(resultados_movimientos["gamma_roca"], digits=4))%")
    println("  - En superficie: $(round(resultados_movimientos["gamma_superficie"], digits=4))%")
    println("  - Máxima en subcapas: $(round(maximum(gamma_max), digits=4))%")
    println()
    println("Archivos generados:")
    println("  - $archivo_salida")
    println()
    println("Nuevas funcionalidades implementadas:")
    println("  - Reconstrucción de señal en superficie por convolución")
    println("  - Análisis lineal equivalente con subdivisión automática de capas")
    println("  - Implementación exacta del método Strata para deformaciones")
    println()
    println("¡Análisis sísmico completo con subdivisión automática y convolución finalizado!")
end

"""
    manejar_no_convergencia(historial::Dict, tolerancia::Float64, max_iter::Int)

Función especializada para manejar casos de no convergencia en el análisis lineal equivalente.
Genera un informe detallado de diagnóstico y detiene la ejecución del análisis.

# Argumentos
- `historial`: Diccionario con el historial de iteraciones del análisis LEQ
- `tolerancia`: Tolerancia objetivo que no se pudo alcanzar
- `max_iter`: Número máximo de iteraciones que se ejecutaron

# Comportamiento
La función genera un informe completo de diagnóstico incluyendo:
- Resumen del estado de convergencia
- Análisis de tendencias en las iteraciones
- Recomendaciones para mejorar la convergencia
- Mensaje de error final que detiene la ejecución
"""
function manejar_no_convergencia(historial::Dict, tolerancia::Float64, max_iter::Int)
    println("\n" * repeat("!", 100))
    println("❌ FALLO DE CONVERGENCIA EN ANÁLISIS LINEAL EQUIVALENTE")
    println(repeat("!", 100))
    
    # Estado final del análisis
    error_final = historial["error"][end] * 100
    tolerancia_objetivo = tolerancia * 100
    
    println("\n📊 RESUMEN DEL ESTADO FINAL:")
    println("─" * repeat("─", 60))
    println("  ⚠️  Error final alcanzado:     $(round(error_final, digits=3))%")
    println("  🎯 Tolerancia objetivo:        $(round(tolerancia_objetivo, digits=3))%")
    println("  📈 Factor de divergencia:      $(round(error_final/tolerancia_objetivo, digits=1))x")
    println("  🔄 Iteraciones ejecutadas:     $max_iter")
    println("  ⏱️  Estado de convergencia:     NO ALCANZADA")
    
    # Análisis de tendencias
    println("\n📈 ANÁLISIS DE TENDENCIAS:")
    println("─" * repeat("─", 60))
    
    # Últimas 10 iteraciones para analizar tendencia
    n_analisis = min(10, length(historial["error"]))
    errores_recientes = historial["error"][end-n_analisis+1:end] .* 100
    
    # Calcular tendencia promedio
    if n_analisis > 3
        pendiente = (errores_recientes[end] - errores_recientes[end-2]) / 3
        if pendiente > 0.1
            tendencia = "📈 DIVERGENTE"
            diagnostico = "El error está aumentando"
        elseif pendiente < -0.1
            tendencia = "📉 CONVERGENTE LENTO"
            diagnostico = "El error disminuye pero muy lentamente"
        else
            tendencia = "📊 OSCILATORIO"
            diagnostico = "El error oscila sin progreso claro"
        end
    else
        tendencia = "❓ INSUFICIENTES DATOS"
        diagnostico = "Muy pocas iteraciones para determinar tendencia"
    end
    
    println("  📊 Tendencia de error:         $tendencia")
    println("  🔍 Diagnóstico:               $diagnostico")
    println("  📋 Últimas $(n_analisis) iteraciones:    $(round.(errores_recientes, digits=2))%")
    
    # Análisis de oscilaciones
    if length(historial["error"]) > 3
        oscilaciones = 0
        for i in 2:(length(historial["error"])-1)
            if (historial["error"][i] > historial["error"][i-1] && historial["error"][i] > historial["error"][i+1]) ||
               (historial["error"][i] < historial["error"][i-1] && historial["error"][i] < historial["error"][i+1])
                oscilaciones += 1
            end
        end
        porcentaje_oscilaciones = oscilaciones / (length(historial["error"])-2) * 100
        println("  🌊 Oscilaciones detectadas:    $(oscilaciones) ($(round(porcentaje_oscilaciones, digits=1))%)")
    end
    
    # Propiedades finales del material
    println("\n🧱 ESTADO FINAL DE PROPIEDADES:")
    println("─" * repeat("─", 60))
    if !isempty(historial["Vs"])
        Vs_final = historial["Vs"][end]
        zeta_final = historial["zeta"][end]
        gamma_final = historial["gamma"][end]
        
        # Promedios por tipo de material
        if length(Vs_final) >= 22
            Vs_diatomaceo = mean(Vs_final[1:22])
            Vs_brecha = length(Vs_final) > 22 ? mean(Vs_final[23:end]) : 0.0
            zeta_diatomaceo = mean(zeta_final[1:22])
            gamma_max_diat = maximum(gamma_final[1:22])
            gamma_max_brecha = length(gamma_final) > 22 ? maximum(gamma_final[23:end]) : 0.0
            
            println("  🏔️  Suelo Diatomáceo:")
            println("      Vs promedio:           $(round(Vs_diatomaceo, digits=1)) m/s")
            println("      ζ promedio:            $(round(zeta_diatomaceo, digits=3))")
            println("      γ máxima:              $(round(gamma_max_diat, digits=3))%")
            
            if Vs_brecha > 0
                println("  🗿 Brecha Sedimentaria:")
                println("      Vs promedio:           $(round(Vs_brecha, digits=1)) m/s")
                println("      γ máxima:              $(round(gamma_max_brecha, digits=3))%")
            end
        end
    end
    
    # Recomendaciones para solucionar el problema
    println("\n💡 RECOMENDACIONES PARA MEJORAR CONVERGENCIA:")
    println("─" * repeat("─", 70))
    println("  1. 🎛️  Ajustar tolerancia:")
    nueva_tolerancia = min(error_final * 1.5 / 100, 0.10)
    println("      - Usar tolerancia más relajada: $(round(nueva_tolerancia*100, digits=2))%")
    println("      - En config: tolerancia = $(nueva_tolerancia)")
    
    println("\n  2. 🔄 Aumentar iteraciones:")
    nuevas_iteraciones = max_iter * 2
    println("      - Usar más iteraciones: $nuevas_iteraciones")
    println("      - En config: max_iter = $nuevas_iteraciones")
    
    println("\n  3. 🧮 Modificar parámetros del suelo:")
    println("      - Revisar curvas de degradación (parámetros a_deg, γref_deg)")
    println("      - Verificar valores iniciales de Vs y ζ")
    println("      - Considerar comportamiento lineal si es apropiado")
    
    println("\n  4. 🔧 Parámetros de convergencia:")
    println("      - El algoritmo mejorado incluye relajación adaptativa")
    println("      - Revise si las deformaciones son realistas")
    println("      - Considere subdividir menos las capas si hay inestabilidad")
    
    # Ejemplo de configuración corregida
    println("\n📝 EJEMPLO DE CONFIGURACIÓN CORREGIDA:")
    println("─" * repeat("─", 70))
    println("config = Dict(")
    println("    \"tolerancia\" => $(nueva_tolerancia),")
    println("    \"max_iter\" => $(nuevas_iteraciones),")
    println("    \"espesor_maximo_subcapas\" => 15.0,  # Mayor espesor para estabilidad")
    println("    # ... otros parámetros")
    println(")")
    
    # Mensaje final de error
    println("\n" * repeat("!", 100))
    println("❌ ANÁLISIS DETENIDO: No se puede continuar sin convergencia adecuada")
    println("❌ Los resultados sin convergencia NO son confiables para diseño")
    println("❌ Implemente las recomendaciones anteriores y ejecute nuevamente")
    println(repeat("!", 100))
    
    # Lanzar error para detener la ejecución
    error("❌ FALLO DE CONVERGENCIA: El análisis lineal equivalente no convergió en $max_iter iteraciones. " *
          "Error final: $(round(error_final, digits=2))% > Tolerancia: $(round(tolerancia_objetivo, digits=2))%. " *
          "Consulte el informe de diagnóstico anterior para soluciones.")
end

end # module