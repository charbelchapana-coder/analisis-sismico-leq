"""
Lib_Preparacion.jl - Librería para preparación del análisis
Etapa 3: Preparación de parámetros para el análisis
"""

module Lib_Preparacion

using LinearAlgebra, FFTW

export preparar_analisis

"""
    preparar_analisis(tiempo::Vector{Float64}, aceleracion_escalada::Vector{Float64}, dt::Float64)

Función principal para la preparación del análisis.
Implementa la Etapa 3 del análisis.
"""
function preparar_analisis(tiempo::Vector{Float64}, aceleracion_escalada::Vector{Float64}, dt::Float64)
    println("\n3. PREPARACIÓN PARA ANÁLISIS")
    println("------------------------------------------------------------")
    
    # Configurar rango de frecuencias para análisis
    fs = 1.0 / dt
    N = length(aceleracion_escalada)
    
    # Frecuencias para análisis (usar toda la resolución disponible)
    freq_min = 0.01
    freq_max = min(99.99, fs/2 - 0.01)
    
    # Usar la resolución natural de la FFT para preservar toda la información
    freq_natural = [i * fs / N for i in 0:(div(N,2)-1)]
    
    # Filtrar frecuencias dentro del rango de interés
    indices_validos = (freq_natural .>= freq_min) .& (freq_natural .<= freq_max)
    freq = freq_natural[indices_validos]
    omega = 2π * freq
    
    # Transformada de Fourier de la aceleración
    f_acc = fft(aceleracion_escalada)
    f_acc_pos = f_acc[1:div(N,2)]
    
    # Usar solo las frecuencias válidas (sin interpolación para preservar información)
    f_acc_interp = f_acc_pos[indices_validos]
    
    println("Parámetros de análisis:")
    println("  - Rango de frecuencias: $(freq_min) - $(freq_max) Hz")
    println("  - Número de frecuencias: $(length(freq))")
    println("  - Resolución en frecuencia: $(round((freq_max - freq_min) / length(freq), digits=4)) Hz")
    
    return freq, omega, f_acc_interp
end

end # module