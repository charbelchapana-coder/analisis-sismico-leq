"""
Lib_Espectros.jl - Librería para cálculo de espectros de respuesta
Etapa 6: Cálculo de espectros de respuesta y análisis espectral
"""

module Lib_Espectros

using LinearAlgebra, FFTW, Statistics

export calcular_espectros_y_amplificacion, espectro_chopra, espectro_fourier

"""
    espectro_fourier(señal::Vector{Float64}, dt::Float64)

Calcula el espectro de Fourier de una señal.
"""
function espectro_fourier(señal::Vector{Float64}, dt::Float64)
    N = length(señal)
    
    # FFT de la señal
    Y = fft(señal)
    
    # Magnitudes (solo mitad positiva del espectro)
    magnitudes = abs.(Y[1:div(N,2)])
    
    # Frecuencias correspondientes
    fs = 1.0 / dt
    frecuencias = [i * fs / N for i in 0:(div(N,2)-1)]
    
    return frecuencias, magnitudes
end

"""
    espectro_chopra(aceleracion, dt, periodos, xi)

Calcula el espectro de respuesta usando el método de Chopra.
Implementación estándar basada en integración paso a paso de Newmark.
"""
function espectro_chopra(aceleracion, dt, periodos, xi)
    Sa = zeros(length(periodos))
    Sv = zeros(length(periodos))
    Sd = zeros(length(periodos))
    m = 1.0
    
    for k in eachindex(periodos)
        T = periodos[k]
        if T == 0
            Sa[k] = maximum(abs.(aceleracion))
            Sv[k] = 0.0
            Sd[k] = 0.0
            continue
        end
        
        w_n = 2 * π / T
        K = w_n^2 * m
        w_d = w_n * sqrt(1 - xi^2)
        e_xi_wn_dt = exp(-xi * w_n * dt)
        sin_wd_dt = sin(w_d * dt)
        cos_wd_dt = cos(w_d * dt)

        # Constantes de integración de Newmark (método de Chopra)
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

        # Integración paso a paso
        for i in 1:(length(aceleracion) - 1)
            u[i+1] = A * u[i] + B * u_dot[i] + C * p[i] + D * p[i+1]
            u_dot[i+1] = Aprime * u[i] + Bprime * u_dot[i] + Cprime * p[i] + Dprime * p[i+1]
        end

        # Cálculo de la aceleración total del sistema SDOF
        acel_total = -(2 * xi * w_n .* u_dot + w_n^2 .* u)
        
        # Máximos de respuesta
        Sa[k] = maximum(abs.(acel_total))  # Pseudo-aceleración (método correcto de Chopra)
        Sv[k] = w_n * maximum(abs.(u))     # Pseudo-velocidad
        Sd[k] = maximum(abs.(u))           # Desplazamiento
    end
    
    return Sa, Sv, Sd
end

"""
    calcular_espectros_y_amplificacion(tiempo::Vector{Float64}, resultados_movimientos::Dict, 
                                      freq::Vector{Float64})

Función principal para el cálculo de espectros de respuesta.
Implementa la Etapa 6 del análisis con optimizaciones anidadas.
"""
function calcular_espectros_y_amplificacion(tiempo::Vector{Float64}, resultados_movimientos::Dict, 
                                           freq::Vector{Float64})
    println("\n6. CÁLCULO DE ESPECTROS DE RESPUESTA")
    println("------------------------------------------------------------")
    println("Calculando espectros de respuesta para máxima resolución...")
    println("Usando todo el dominio de frecuencias disponible")
    
    # Función anidada para configuración de períodos optimizada
    function configurar_periodos_analisis(freq_vector)
        periodos_completos = 1 ./ freq_vector[2:end]  # Excluir frecuencia cero
        indices_validos = (periodos_completos .>= 0.01) .& (periodos_completos .<= 10.0)
        periodos_filtrados = periodos_completos[indices_validos]
        
        # Agregar períodos críticos adicionales
        periodos_criticos = [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0]
        periodos_union = sort(unique(vcat(periodos_filtrados, periodos_criticos)))
        
        return periodos_union
    end
    
    # Función anidada para cálculo optimizado de espectros múltiples
    function calcular_espectros_multiples(tiempo_vec, acc_rock, acc_superficie, acc_conv, acc_alt, periodos)
        dt = tiempo_vec[2] - tiempo_vec[1]
        xi = 0.05  # 5% amortiguamiento crítico
        
        # Calcular espectros en paralelo conceptual (secuencial optimizado)
        println("Calculando espectro de respuesta para señal en roca...")
        Sa_rock, Sv_rock, Sd_rock = espectro_chopra(acc_rock, dt, periodos, xi)
        
        println("Calculando espectro de respuesta para señal en superficie...")
        Sa_superficie, Sv_superficie, Sd_superficie = espectro_chopra(acc_superficie, dt, periodos, xi)
        
        println("Calculando espectro de respuesta para señal de convolución...")
        Sa_convolucion, Sv_convolucion, Sd_convolucion = espectro_chopra(acc_conv, dt, periodos, xi)
        
        println("Calculando espectro de respuesta para señal de convolución alternativa...")
        Sa_alternativa, Sv_alternativa, Sd_alternativa = espectro_chopra(acc_alt, dt, periodos, xi)
        
        return (Sa_rock, Sv_rock, Sd_rock), (Sa_superficie, Sv_superficie, Sd_superficie), (Sa_convolucion, Sv_convolucion, Sd_convolucion), (Sa_alternativa, Sv_alternativa, Sd_alternativa)
    end
    
    # Configurar períodos de análisis
    periodos_analisis = configurar_periodos_analisis(freq)
    println("Rango de períodos: $(minimum(periodos_analisis)) - $(maximum(periodos_analisis)) s")
    println("Número total de puntos: $(length(periodos_analisis))")
    
    # Calcular todos los espectros
    espectros_roca, espectros_superficie, espectros_convolucion, espectros_alternativa = calcular_espectros_multiples(
        resultados_movimientos["tiempo_adj"],
        resultados_movimientos["acc_rock"],
        resultados_movimientos["acc_superficie"],
        resultados_movimientos["acc_superficie_convolucion"],
        resultados_movimientos["acc_superficie_alternativa"],
        periodos_analisis
    )
    
    println("Espectros de respuesta calculados exitosamente (original, convolución y alternativa)")
    
    return Dict(
        "periodos" => periodos_analisis,
        "Sa_rock" => espectros_roca[1],
        "Sv_rock" => espectros_roca[2],
        "Sd_rock" => espectros_roca[3],
        "Sa_superficie" => espectros_superficie[1],
        "Sv_superficie" => espectros_superficie[2],
        "Sd_superficie" => espectros_superficie[3],
        "Sa_convolucion" => espectros_convolucion[1],
        "Sv_convolucion" => espectros_convolucion[2],
        "Sd_convolucion" => espectros_convolucion[3],
        "Sa_alternativa" => espectros_alternativa[1],
        "Sv_alternativa" => espectros_alternativa[2],
        "Sd_alternativa" => espectros_alternativa[3]
    )
end

end # module