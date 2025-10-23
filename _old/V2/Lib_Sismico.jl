"""
Lib_Sismico.jl - Librería para procesamiento de datos sísmicos
Etapa 1: Lectura y procesamiento de datos sísmicos
"""

module Lib_Sismico

using LinearAlgebra, FFTW, Statistics, DelimitedFiles

export leer_sismo, escalar_sismo, periodo_predominante, procesar_datos_sismicos

"""
    leer_sismo(archivo::String)

Lee un archivo de sismo y retorna tiempo, aceleración y dt.
"""
function leer_sismo(archivo::String)
    # Leer archivo con coma como delimitador de columnas y punto como separador decimal
    datos = readdlm(archivo, ',')
    
    if size(datos, 2) == 1
        aceleracion = Float64.(datos[:, 1])
        dt = 0.005  # Valor por defecto
        tiempo = [Float64(i * dt) for i in 0:(length(aceleracion)-1)]
    else
        tiempo = Float64.(datos[:, 1])
        aceleracion = Float64.(datos[:, 2])
        dt = tiempo[2] - tiempo[1]
    end
    
    return tiempo, aceleracion, dt
end

"""
    escalar_sismo(dt::Float64, aceleracion::Vector{Float64}, factor::Float64)

Escala una señal sísmica por un factor dado y calcula el período predominante.
Implementa la Etapa 1 del análisis con optimización de escalado y análisis espectral.
"""
function escalar_sismo(dt::Float64, aceleracion::Vector{Float64}, factor::Float64)
    
    # Función anidada para validación de parámetros
    function validar_parametros_escalado(dt_val, acel_vec, factor_val)
        if dt_val <= 0
            throw(ArgumentError("El intervalo de tiempo dt debe ser positivo"))
        end
        if isempty(acel_vec)
            throw(ArgumentError("El vector de aceleración no puede estar vacío"))
        end
        if factor_val <= 0
            throw(ArgumentError("El factor de escala debe ser positivo"))
        end
    end
    
    # Función anidada para cálculo de período predominante
    function calcular_periodo_predominante(acel, dt_calc)
        N = length(acel)
        if N < 4  # Mínimo para FFT útil
            return 1.0  # Valor por defecto
        end
        
        # FFT optimizada
        acel_fft = fft(acel .- mean(acel))  # Remover media
        freqs = fftfreq(N, 1/dt_calc)
        
        # Solo frecuencias positivas
        n_half = div(N, 2)
        mag = abs.(acel_fft[1:n_half])
        freqs_pos = freqs[1:n_half]
        
        # Encontrar frecuencia dominante (excluyendo DC)
        if length(mag) > 1
            idx_max = argmax(mag[2:end]) + 1  # +1 para compensar el offset
            freq_dom = freqs_pos[idx_max]
            return freq_dom > 0 ? 1.0 / freq_dom : 1.0
        else
            return 1.0
        end
    end
    
    # Validar parámetros
    validar_parametros_escalado(dt, aceleracion, factor)
    
    println("   Escalando señal sísmica...")
    println("   Factor de escala: $factor")
    println("   Puntos originales: $(length(aceleracion))")
    
    # Escalar aceleración
    aceleracion_escalada = aceleracion * factor
    
    # Calcular período predominante
    T_p = calcular_periodo_predominante(aceleracion_escalada, dt)
    
    println("   Período predominante: $(round(T_p, digits=3)) s")
    
    return dt, aceleracion_escalada, T_p
end

"""
    periodo_predominante(aceleracion::Vector{Float64}, dt::Float64)

Calcula el período predominante de una señal sísmica usando FFT.
"""
function periodo_predominante(aceleracion::Vector{Float64}, dt::Float64)
    N = length(aceleracion)
    
    # Calcular FFT
    A_fft = fft(aceleracion)
    A_mag = abs.(A_fft[1:div(N,2)])
    
    # Frecuencias
    fs = 1.0 / dt
    freq = [i * fs / N for i in 0:(div(N,2)-1)]
    
    # Encontrar frecuencia dominante (excluyendo frecuencia cero)
    idx_max = argmax(A_mag[2:end]) + 1
    freq_dom = freq[idx_max]
    
    # Período dominante
    T_dom = 1.0 / freq_dom
    
    return T_dom
end

"""
    procesar_datos_sismicos(archivo_sismo::String, factor_escala::Float64)

Función principal para el procesamiento completo de datos sísmicos.
Implementa la Etapa 1 del análisis.
"""
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
    
    # PGA inicial
    PGA_inicial = maximum(abs.(aceleracion))
    println("  - PGA inicial: $(round(PGA_inicial, digits=4)) m/s²")
    
    # Escalar sismo
    println("  - Factor de escala: $factor_escala")
    aceleracion_escalada = escalar_sismo(aceleracion, factor_escala)
    
    # PGA escalada
    PGA_escalada = maximum(abs.(aceleracion_escalada))
    println("  - PGA escalada: $(round(PGA_escalada, digits=4)) m/s²")
    
    # Período predominante
    T_p = periodo_predominante(aceleracion_escalada, dt)
    println("  - Período predominante: $(round(T_p, digits=3)) segundos")
    
    return tiempo, aceleracion_escalada, dt, T_p
end

end # module