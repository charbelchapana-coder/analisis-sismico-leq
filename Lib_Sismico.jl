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
Detecta automáticamente el separador del archivo (coma, tab, punto y coma, espacio).
"""
function leer_sismo(archivo::String)
    # Función auxiliar para detectar el separador
    function detectar_separador(linea::String)
        separadores = ['\t', ',', ';', ' ']
        conteos = [count(x -> x == sep, linea) for sep in separadores]
        max_idx = argmax(conteos)
        
        # Si el separador más común tiene al menos 1 ocurrencia, usarlo
        if conteos[max_idx] > 0
            return separadores[max_idx]
        else
            # Si no hay separadores obvios, asumir espacio
            return ' '
        end
    end
    
    # Leer las primeras líneas para detectar el formato
    lineas_muestra = []
    open(archivo, "r") do file
        for i in 1:min(10, countlines(archivo))  # Leer máximo 10 líneas o menos si el archivo es pequeño
            linha = readline(file)
            if !isempty(strip(linha))  # Ignorar líneas vacías
                push!(lineas_muestra, linha)
            end
        end
    end
    
    if isempty(lineas_muestra)
        throw(ArgumentError("El archivo está vacío o no contiene datos válidos"))
    end
    
    # Detectar separador basándose en la primera línea no vacía
    separador = detectar_separador(lineas_muestra[1])
    
    println("  - Separador detectado: '$(separador == '\t' ? "\\t (tab)" : separador == ' ' ? "espacio" : string(separador))'")
    
    # Leer archivo con el separador detectado
    try
        datos = readdlm(archivo, separador)
        
        # Verificar que los datos se leyeron correctamente
        if size(datos, 1) == 0
            throw(ArgumentError("No se pudieron leer datos del archivo"))
        end
        
        # Convertir a matriz de números si es necesario
        if eltype(datos) == Any
            # Intentar convertir cada elemento a Float64
            datos_numericos = similar(datos, Float64)
            for i in 1:size(datos, 1), j in 1:size(datos, 2)
                try
                    if isa(datos[i, j], AbstractString)
                        # Limpiar string y convertir
                        str_limpio = replace(string(datos[i, j]), r"[^\d.eE+-]" => "")
                        datos_numericos[i, j] = parse(Float64, str_limpio)
                    else
                        datos_numericos[i, j] = Float64(datos[i, j])
                    end
                catch
                    throw(ArgumentError("Error al convertir dato en fila $i, columna $j: '$(datos[i, j])'"))
                end
            end
            datos = datos_numericos
        end
        
        if size(datos, 2) == 1
            # Si solo hay una columna, asumir que son solo aceleraciones
            aceleracion = Float64.(datos[:, 1])
            dt = 0.005  # Valor por defecto
            tiempo = [Float64(i * dt) for i in 0:(length(aceleracion)-1)]
            println("  - Formato detectado: Solo aceleraciones (dt = $dt s)")
        else
            # Si hay dos o más columnas: tiempo y aceleración
            tiempo = Float64.(datos[:, 1])
            aceleracion = Float64.(datos[:, 2])
            dt = tiempo[2] - tiempo[1]
            println("  - Formato detectado: Tiempo y aceleración (dt = $(round(dt, digits=6)) s)")
        end
        
        return tiempo, aceleracion, dt
        
    catch e
        if isa(e, ArgumentError)
            rethrow(e)
        else
            throw(ArgumentError("Error al leer el archivo con separador '$separador': $e"))
        end
    end
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