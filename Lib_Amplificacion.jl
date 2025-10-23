"""
Lib_Amplificacion.jl - Librería para funciones de amplificación
Etapa 7: Cálculo de funciones de amplificación y análisis espectral
"""

module Lib_Amplificacion

using LinearAlgebra, FFTW, Statistics
using ..Lib_TransferenciaCore: capa_suelo, parametros_degradacion

export calcular_funciones_amplificacion, funcion_amplificacion

"""
    funcion_amplificacion(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                         capa_referencia::Int, capa_objetivo::Int)

Calcula la función de amplificación entre dos capas usando análisis de transferencia.
"""
function funcion_amplificacion(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                              capa_referencia::Int, capa_objetivo::Int)
    
    # Función anidada para cálculo optimizado de transferencia
    function calcular_transferencia_capa(capas_vec, omega_vec, capa_idx)
        if capa_idx == length(capas_vec)  # Roca
            return ones(ComplexF64, length(omega_vec))
        else
            return funcion_transferencia(capas_vec[1:capa_idx], omega_vec)
        end
    end
    
    # Función de transferencia hasta capa de referencia
    H_ref = calcular_transferencia_capa(capas, omega, capa_referencia)
    
    # Función de transferencia hasta capa objetivo
    H_obj = calcular_transferencia_capa(capas, omega, capa_objetivo)
    
    # Función de amplificación
    FA = H_obj ./ H_ref
    
    return FA
end

"""
    funcion_amplificacion(acc_superficie::Vector{Float64}, acc_roca::Vector{Float64}, dt::Float64)

Calcula la función de amplificación a partir de señales de aceleración en el dominio del tiempo.
Versión optimizada con funciones anidadas para procesamiento eficiente.
"""
function funcion_amplificacion(acc_superficie::Vector{Float64}, acc_roca::Vector{Float64}, dt::Float64)
    
    # Función anidada para validación y preprocesamiento
    function validar_y_preparar_señales(señal1, señal2, dt_val)
        if length(señal1) != length(señal2)
            error("Las señales de superficie y roca deben tener la misma longitud")
        end
        
        N = length(señal1)
        fs = 1.0 / dt_val
        freq = fftfreq(N, fs)
        n_freq = div(N, 2)
        
        return N, fs, freq[1:n_freq]
    end
    
    # Función anidada para cálculo FFT optimizado
    function calcular_ffts_y_amplificacion(superficie, roca, freq_positivas)
        # Calcular FFT de ambas señales
        fft_superficie = fft(superficie)
        fft_roca = fft(roca)
        
        # Tomar solo la mitad positiva del espectro
        n_freq = length(freq_positivas)
        fft_sup_pos = fft_superficie[1:n_freq]
        fft_roca_pos = fft_roca[1:n_freq]
        
        # Calcular función de amplificación con protección contra división por cero
        amplificacion = similar(fft_sup_pos, Float64)
        
        for i in eachindex(amplificacion)
            if abs(fft_roca_pos[i]) > 1e-10
                amplificacion[i] = abs(fft_sup_pos[i]) / abs(fft_roca_pos[i])
            else
                amplificacion[i] = 1.0
            end
        end
        
        return freq_positivas, amplificacion
    end
    
    # Función anidada para filtrado de frecuencias
    function filtrar_frecuencias_validas(freq_vec, amp_vec)
        # Evitar la frecuencia cero y muy bajas
        indices_validos = freq_vec .> 0.01  # Frecuencias > 0.01 Hz
        
        return freq_vec[indices_validos], amp_vec[indices_validos]
    end
    
    # Ejecutar procesamiento optimizado
    N, fs, freq_pos = validar_y_preparar_señales(acc_superficie, acc_roca, dt)
    freq_result, amplificacion_raw = calcular_ffts_y_amplificacion(acc_superficie, acc_roca, freq_pos)
    freq_final, amplificacion_final = filtrar_frecuencias_validas(freq_result, amplificacion_raw)
    
    return freq_final, amplificacion_final
end

"""
    calcular_funciones_amplificacion(tiempo::Vector{Float64}, resultados_movimientos::Dict, 
                                    freq::Vector{Float64})

Función principal para el cálculo de funciones de amplificación.
Implementa la Etapa 7 del análisis con optimizaciones integradas.
"""
function calcular_funciones_amplificacion(tiempo::Vector{Float64}, resultados_movimientos::Dict, 
                                         freq::Vector{Float64})
    println("\n7. CÁLCULO DE FUNCIONES DE AMPLIFICACIÓN")
    println("------------------------------------------------------------")
    
    # Función anidada para análisis estadístico de amplificación
    function analizar_amplificacion(freq_amp, amp_valores)
        max_amp = maximum(amp_valores)
        idx_max = argmax(amp_valores)
        freq_max = freq_amp[idx_max]
        periodo_max = 1.0 / freq_max
        
        return max_amp, freq_max, periodo_max
    end
    
    # Función anidada para cálculo múltiple de amplificaciones
    function calcular_amplificaciones_multiples(resultados)
        dt = resultados["tiempo_adj"][2] - resultados["tiempo_adj"][1]
        
        # Amplificación principal (superficie vs roca)
        freq_amp_principal, amp_principal = funcion_amplificacion(
            resultados["acc_superficie"], 
            resultados["acc_rock"], 
            dt
        )
        
        # Amplificación de verificación (convolución vs roca)
        freq_amp_verificacion, amp_verificacion = funcion_amplificacion(
            resultados["acc_superficie_convolucion"], 
            resultados["acc_rock"], 
            dt
        )
        
        return (freq_amp_principal, amp_principal), (freq_amp_verificacion, amp_verificacion)
    end
    
    # Calcular amplificaciones
    (freq_amp, amplificacion), (freq_amp_verif, amp_verif) = calcular_amplificaciones_multiples(resultados_movimientos)
    
    # Análisis estadístico
    max_amp, freq_max, periodo_max = analizar_amplificacion(freq_amp, amplificacion)
    
    println("Función de amplificación calculada:")
    println("  - Amplificación máxima: $(round(max_amp, digits=1)) @ $(round(freq_max, digits=2)) Hz")
    println("  - Período de máxima amplificación: $(round(periodo_max, digits=1)) s")
    
    # Verificación de consistencia
    if length(amp_verif) > 0
        max_amp_verif = maximum(amp_verif)
        diferencia_relativa = abs(max_amp - max_amp_verif) / max_amp * 100
        
        if diferencia_relativa < 10.0
            println("  ✓ Verificación exitosa: Diferencia entre métodos < 10%")
        else
            println("  ⚠️ Verificación: Diferencia entre métodos = $(round(diferencia_relativa, digits=1))%")
        end
    end
    
    return Dict(
        "freq_amp" => freq_amp,
        "amplificacion" => amplificacion,
        "freq_amp_verificacion" => freq_amp_verif,
        "amplificacion_verificacion" => amp_verif,
        "max_amplificacion" => max_amp,
        "freq_max_amplificacion" => freq_max,
        "periodo_max_amplificacion" => periodo_max
    )
end

end # module