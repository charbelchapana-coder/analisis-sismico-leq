# ===================================================================
# SCRIPT ANÁLISIS SÍSMICO - VERSIÓN RECUPERADA CON LIBRERÍA ORIGINAL
# ===================================================================

# Cargar la librería original
include("_old/Lib_LEQ.jl")
using .Lib_LEQ

using FFTW, Statistics, Printf

println("INICIANDO ANÁLISIS SÍSMICO CON BIBLIOTECA ORIGINAL")
println(repeat("=", 60))

# ===== CONFIGURACIÓN =====
archivo_sismo = "sismo_recortado.txt"
factor_escala = 1.0
max_iteraciones = 15
tolerancia_convergencia = 0.1  # 10%

# ===== 1. CARGA Y PROCESAMIENTO DEL SISMO =====
println("\\n1. CARGA Y PROCESAMIENTO DEL SISMO")
println(repeat("-", 60))

tiempo, acc_entrada, dt_sismo = leer_sismo(archivo_sismo)
acc_entrada = escalar_sismo(acc_entrada, factor_escala)
T_p = periodo_predominante(acc_entrada, dt_sismo)

println("  Archivo: $archivo_sismo")
println("  Factor de escala: $factor_escala")
println("  Puntos originales: $(length(acc_entrada))")
println("  Período predominante: $(round(T_p, digits=3)) s")
println("  PGA señal de entrada: $(round(maximum(abs.(acc_entrada)), digits=4)) m/s²")

# ===== 2. PROCESAMIENTO DEL MODELO DE SUELO =====
println("\\n2. PROCESAMIENTO DEL MODELO DE SUELO")
println(repeat("=", 52))

# Parámetros de degradación para suelo diatomáceo
params_diat = parametros_degradacion(0.5, 0.01, 1.0, 0.3, 0.02, "equivalente")

# Parámetros de degradación para brecha (comportamiento más rígido)
params_brecha = parametros_degradacion(0.5, 0.01, 1.0, 0.3, 0.02, "equivalente")

# Crear capas principales simplificadas
subcapas = [
    capa_suelo(1400.0, 500.0, 0.05, 220.0, params_diat, "Suelo diatomáceo"),
    capa_suelo(2000.0, 850.0, 0.03, 280.0, params_brecha, "Brecha sedimentaria"),
    capa_suelo(2500.0, 2050.0, 0.02, 0.0, nothing, "Roca basal")  # Roca basal
]

# Subdividir capas para mejor resolución
subcapas = subdividir_capas(subcapas, 10.0)  # Subcapas de máximo 10m

# Mostrar información de capas
for (i, capa) in enumerate(subcapas)
    if capa.h > 0
        println("  Capa $i: $(capa.id) - h=$(capa.h)m, Vs=$(capa.Vs_inicial)m/s, ζ=$(capa.ζ_inicial)")
    else
        println("  Capa $i: $(capa.id) - Roca basal")
    end
end

println("\\n  Subdividiendo capas...")
total_subcapas = length(subcapas) - 1  # Excluir roca basal
altura_total = sum(capa.h for capa in subcapas[1:end-1])

println("  Total de subcapas: $total_subcapas")
println("  Altura total del depósito: $(altura_total) m")

# Calcular período fundamental
T_fundamental, H_total, Vs_prom = periodo_fundamental(subcapas)
println("  Período fundamental: $(round(T_fundamental, digits=3)) s")

# ===== 3. PREPARACIÓN PARA ANÁLISIS =====
println("\\n3. PREPARACIÓN PARA ANÁLISIS")
println(repeat("-", 60))

# Preparar datos de frecuencia
N = length(acc_entrada)
ω = 2π * fftfreq(N, 1/dt_sismo) |> Vector{Float64}

# Filtrar frecuencias positivas y ajustar rango
ω_pos = ω[1:div(N,2)+1]
ω_pos[1] = 0.01  # Evitar frecuencia cero

println("  - Rango de frecuencias: 0.01 - $(round(ω_pos[end]/(2π), digits=1)) Hz")
println("  - Número de frecuencias: $(length(ω_pos))")
println("  - Resolución en frecuencia: $(round((ω_pos[end]-ω_pos[1])/(2π*length(ω_pos)), digits=4)) Hz")

# FFT del sismo
aceleracion_base_fft = fft(acc_entrada)
aceleracion_base_fft_pos = aceleracion_base_fft[1:div(N,2)+1]

# ===== 4. ANÁLISIS LINEAL EQUIVALENTE =====
println("\\n4. ANÁLISIS LINEAL EQUIVALENTE")
println(repeat("-", 60))

try
    # Llamar al análisis con parámetros correctos
    capas_finales, convergencia, historial = analisis_lineal_equivalente(
        subcapas, 
        ω_pos, 
        aceleracion_base_fft_pos, 
        max_iteraciones, 
        tolerancia_convergencia
    )
    
    println(repeat("-", 100))
    if convergencia
        println("✅ ¡Análisis convergió exitosamente!")
        println("  Convergencia alcanzada: Sí")
        println("  Iteraciones realizadas: $(length(historial["error"]))")
    else
        println("⚠️  Advertencia: No se alcanzó convergencia en $max_iteraciones iteraciones")
        if !isempty(historial["error"])
            println("   Error final: $(round(historial["error"][end]*100, digits=1))% (tolerancia: $(tolerancia_convergencia*100)%)")
        end
        println("  Convergencia alcanzada: No")
        println("  Iteraciones realizadas: $(length(historial["error"]))")
    end
    println(repeat("=", 100))
    
    global resultado_completo = (
        capas_finales = capas_finales,
        convergencia = convergencia,
        historial = historial
    )
    
catch e
    println("❌ Error en análisis lineal equivalente: $e")
    global resultado_completo = nothing
end

# ===== 5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES =====
println("\\n5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES")
println(repeat("-", 60))

if resultado_completo !== nothing
    try
        # Usar las capas finales del análisis
        capas_finales = resultado_completo.capas_finales
        
        # Calcular deformaciones de corte
        gamma_max, gamma_superficie, gamma_roca, gamma_interfaces = deformaciones_corte_shake(
            capas_finales, ω_pos, aceleracion_base_fft_pos
        )
        
        println("✅ Deformaciones calculadas exitosamente")
        println("  Deformación máxima en superficie: $(round(gamma_superficie, digits=3))%")
        println("  Deformación máxima en roca: $(round(gamma_roca, digits=3))%")
        println("  Deformación máxima en capas: $(round(maximum(gamma_max), digits=3))%")
        
        global gamma_resultados = (
            gamma_max = gamma_max,
            gamma_superficie = gamma_superficie, 
            gamma_roca = gamma_roca,
            gamma_interfaces = gamma_interfaces
        )
        
    catch e
        println("❌ Error en cálculo de deformaciones: $e")
        global gamma_resultados = nothing
    end
else
    println("⚠️  Omitiendo cálculo de deformaciones debido a error en análisis")
    global gamma_resultados = nothing
end

# ===== 6. GENERACIÓN DE GRÁFICOS =====
println("\\n6. GENERACIÓN DE GRÁFICOS")
println(repeat("-", 60))

if resultado_completo !== nothing && gamma_resultados !== nothing
    try
        # Los gráficos pueden requerir funciones adicionales no disponibles en la librería original
        println("⚠️  Funciones de graficación no disponibles en librería original")
        println("Datos disponibles para graficación externa:")
        println("  - resultado_completo (capas finales, convergencia, historial)")
        println("  - gamma_resultados (deformaciones)")
        
    catch e
        println("❌ Error en generación de gráficos: $e")
    end
else
    println("⚠️  Omitiendo gráficos debido a errores anteriores")
end

println("\\n" * repeat("=", 60))
println("ANÁLISIS COMPLETADO")
println(repeat("=", 60))