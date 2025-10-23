# Script LEQ usando la biblioteca original que funciona
using FFTW, Plots, DelimitedFiles, Statistics

# Incluir la biblioteca original
include("_old/Lib_LEQ.jl")
using .Lib_LEQ

# Configurar backend para gráficos
gr()

# ===== CONFIGURACIÓN DEL ANÁLISIS =====

# Configuración del sismo
archivo_sismo = "sismo_recortado.txt"
factor_escala = 1.0

# Configuración del modelo de suelo (3 capas)
capas_modelo = [
    (1400, 500, 0.05, 220, "Suelo diatomáceo"),     # ρ, Vs, ζ, h, nombre
    (2000, 850, 0.03, 280, "Brecha sedimentaria"),  # ρ, Vs, ζ, h, nombre  
    (2500, 2050, 0.02, 0, "Roca basal")             # ρ, Vs, ζ, h, nombre (h=0 para roca)
]

# Configuración de subdivisión (espesor máximo por subcapa)
espesor_max_subcapa = 10.0  # metros

# Configuración de análisis lineal equivalente
tolerancia_convergencia = 0.1  # 0.1%
max_iteraciones = 15
factor_relajacion = 0.5

println("INICIANDO ANÁLISIS SÍSMICO CON BIBLIOTECA ORIGINAL")
println("============================================================")

# ===== 1. CARGAR Y PROCESAR SISMO =====
println("\\n1. CARGA Y PROCESAMIENTO DEL SISMO")
println("------------------------------------------------------------")

# Cargar el sismo
try
    t, acc_superficie, dt = leer_sismo(archivo_sismo)
    
    # Escalar el sismo
    acc_superficie_escalada = escalar_sismo(acc_superficie, factor_escala)
    
    # Calcular período predominante
    T_p = periodo_predominante(acc_superficie_escalada, dt)
    
    println("  Archivo: $archivo_sismo")
    println("  Factor de escala: $factor_escala")
    println("  Puntos originales: $(length(acc_superficie_escalada))")
    println("  Período predominante: $(round(T_p, digits=3)) s")
    println("  PGA señal de entrada: $(round(maximum(abs.(acc_superficie_escalada)), digits=4)) m/s²")
    
    global tiempo = t
    global acc_entrada = acc_superficie_escalada
    global dt_sismo = dt
    
catch e
    error("Error al cargar el sismo: $e")
end

# ===== 2. CREAR MODELO DE SUELO =====
println("\\n2. PROCESAMIENTO DEL MODELO DE SUELO")
println("====================================================")

# Crear capas iniciales
capas_iniciales = capa_suelo[]  # Vector tipado
for (i, (ρ, Vs, ζ, h, nombre)) in enumerate(capas_modelo)
    # Parámetros de degradación para suelos (no para roca)
    if i < length(capas_modelo)  # No es roca basal
        params = parametros_degradacion(0.5, 0.01, 1.0, 0.3, 0.02, "equivalente")
    else
        params = parametros_degradacion(0.0, 0.01, 0.0, 0.0, 0.0, "lineal")  # Roca = lineal
    end
    
    capa = capa_suelo(ρ, Vs, ζ, h, params, nombre)
    push!(capas_iniciales, capa)
    
    if h > 0
        println("  Capa $i: $nombre - h=$(h)m, Vs=$(Vs)m/s, ζ=$(ζ)")
    else
        println("  Capa $i: $nombre - Roca basal")
    end
end

# Subdividir capas
println("\\n  Subdividiendo capas...")
subcapas = subdividir_capas(capas_iniciales, espesor_max_subcapa)

total_subcapas = length(subcapas) - 1  # Excluir roca basal
altura_total = sum(capa.h for capa in subcapas[1:end-1])

println("  Total de subcapas: $total_subcapas")
println("  Altura total del depósito: $(altura_total) m")

# Calcular período fundamental
T_fundamental, H_total, Vs_prom = periodo_fundamental(subcapas)
println("  Período fundamental: $(round(T_fundamental, digits=3)) s")

# ===== 3. PREPARACIÓN PARA ANÁLISIS =====
println("\\n3. PREPARACIÓN PARA ANÁLISIS")
println("------------------------------------------------------------")

# Calcular parámetros de frecuencia
freq_max = 100.0  # Hz
freq_min = 0.01   # Hz
n_freq = length(acc_entrada) ÷ 2 + 1
freq = range(freq_min, freq_max, length=n_freq)
omega = 2π .* freq

println("  - Rango de frecuencias: $(freq_min) - $(freq_max) Hz")
println("  - Número de frecuencias: $(length(omega))")
println("  - Resolución en frecuencia: $(round((freq_max-freq_min)/length(omega), digits=4)) Hz")

# ===== 4. ANÁLISIS LINEAL EQUIVALENTE =====
println("\\n4. ANÁLISIS LINEAL EQUIVALENTE")
println("------------------------------------------------------------")

try
    # Preparar datos para la función original
    # Calcular FFT de la señal de entrada
    n_puntos = length(acc_entrada)
    aceleracion_fft = fft(acc_entrada)
    omega = 2π * fftfreq(n_puntos, 1/dt_sismo)
    omega = omega[1:div(n_puntos,2)+1]  # Solo frecuencias positivas
    aceleracion_base_fft = aceleracion_fft[1:div(n_puntos,2)+1]
    
    # Llamar función original con parámetros correctos
    capas_finales, convergencia, historial = analisis_lineal_equivalente(subcapas, omega, aceleracion_base_fft, 
                                                                          max_iteraciones, tolerancia_convergencia)
    
    println("✅ ¡Análisis convergió exitosamente!")
    println("  Convergencia alcanzada: $convergencia")
    println("  Iteraciones realizadas: $(length(historial["Vs"]))")
    
    global capas_resultado = capas_finales
    
catch e
    error("Error en análisis lineal equivalente: $e")
end

# ===== 5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES =====
println("\\n5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES")
println("------------------------------------------------------------")

try
    # Calcular deformaciones de corte
    gamma_max, gamma_superficie, gamma_roca, gamma_interfaces = deformaciones_corte_shake(
        capas_resultado, omega, acc_entrada, dt_sismo
    )
    
    println("DEFORMACIONES DE CORTE CALCULADAS:")
    println("  Rango de valores: $(round(minimum(gamma_max), digits=3))% - $(round(maximum(gamma_max), digits=3))%")
    println("  Deformación máxima: $(round(maximum(gamma_max), digits=3))%")
    println("  Roca basal: $(round(gamma_roca, digits=4))%")
    
    # Mostrar primeras deformaciones
    for i in 1:min(10, length(gamma_max))
        println("  $(capas_resultado[i].id): $(round(gamma_max[i], digits=4))%")
    end
    
    global deformaciones_resultado = gamma_max
    
catch e
    println("Error en cálculo de deformaciones: $e")
    global deformaciones_resultado = zeros(length(capas_resultado)-1)
end

# ===== 6. GENERAR GRÁFICOS =====
println("\\n6. GENERACIÓN DE GRÁFICOS")
println("------------------------------------------------------------")

try
    # Crear gráfico de perfil de deformaciones usando la función original
    archivo_salida = "analisis_original_completo.png"
    
    p_deform = graficar_perfil_deformaciones_detallado(
        capas_resultado, deformaciones_resultado, gamma_interfaces, omega, acc_entrada
    )
    
    # Guardar gráfico
    savefig(p_deform, archivo_salida)
    println("  ✅ Gráfico guardado como: $archivo_salida")
    
catch e
    println("Error en generación de gráficos: $e")
    println("Continuando sin gráficos...")
end

println("\\n============================================================")
println("ANÁLISIS COMPLETADO EXITOSAMENTE")
println("============================================================")