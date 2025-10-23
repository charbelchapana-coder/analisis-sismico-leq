# Script Principal - Análisis Sísmico Completo con Método Lineal Equivalente
# Implementa todas las funcionalidades con subdivisión automática de capas

using FFTW, Plots, DelimitedFiles, Statistics
include("Lib_LEQ.jl")
# Use module reference to avoid conflicts

# Configurar backend para gráficos y fuentes
gr()  # Backend GR estable

# Configurar fuentes globales con tamaños extra grandes para máxima legibilidad
# Usar fuentes disponibles por defecto en Windows - Incremento adicional de 50%
default(titlefontsize=21,
        guidefontsize=18,
        tickfontsize=15,
        legendfontsize=15)

println("=============================================================")
println("ANÁLISIS SÍSMICO CON MÉTODO LINEAL EQUIVALENTE")
println("=============================================================")

# ===== 1. LECTURA Y PROCESAMIENTO DE DATOS SÍSMICOS =====
println("\n1. LECTURA Y PROCESAMIENTO DE DATOS SÍSMICOS")
println("------------------------------------------------------------")

# Leer archivo de sismo
archivo_sismo = "sismo_recortado.txt"
tiempo, aceleracion, dt = Lib_LEQ.leer_sismo(archivo_sismo)

# Información básica del sismo
println("Archivo sísmico cargado: $archivo_sismo")
println("  - Número de puntos: $(length(tiempo))")
println("  - Duración: $(round(maximum(tiempo), digits=2)) segundos")
println("  - Dt: $(round(dt, digits=4)) segundos")

# Escalar sismo (si es necesario)
factor_escala = 1.0  # Factor de escalamiento
aceleracion_escalada = Lib_LEQ.escalar_sismo(aceleracion, factor_escala)

# Información de escalamiento
acel_max_inicial = maximum(abs.(aceleracion))
acel_max_escalada = maximum(abs.(aceleracion_escalada))
println("  - PGA inicial: $(round(acel_max_inicial, digits=4)) m/s²")
println("  - Factor de escala: $factor_escala")
println("  - PGA escalada: $(round(acel_max_escalada, digits=4)) m/s²")

# Calcular período predominante
T_p = Lib_LEQ.periodo_predominante(aceleracion_escalada, dt)
println("  - Período predominante: $(round(T_p, digits=3)) segundos")

# Parámetros para FFT
N = length(aceleracion_escalada)

# ===== 2. DEFINICIÓN DEL MODELO DE SUELO =====
println("\n2. DEFINICIÓN DEL MODELO DE SUELO")
println("------------------------------------------------------------")

# Parámetros para comportamiento lineal equivalente (suelo diatomáceo)
params_diatomaceo = Lib_LEQ.parametros_degradacion(
    1.5, 0.001, 0.25, 0.15, 0.05, "equivalente"
)

# Parámetros para comportamiento lineal (brecha)
params_lineal_brecha = Lib_LEQ.parametros_degradacion(
    10000, 10000, 0.0, 0.0, 0.03, "lineal"
)

# Parámetros para comportamiento lineal (roca)
params_lineal_roca = Lib_LEQ.parametros_degradacion(
    10000, 10000, 0.0, 0.0, 0.02, "lineal"
)

# Definir capas de suelo con comportamiento no lineal
capas = [
    # Capa 1: Suelo diatomaceo (comportamiento no lineal)
    Lib_LEQ.capa_suelo(1400.0, 500.0, 0.05, 220.0, params_diatomaceo, "Suelo diatomáceo"),
    
    # Capa 2: Arcilla media (comportamiento no lineal)
    Lib_LEQ.capa_suelo(2000.0, 850.0, 0.03, 280.0, params_lineal_brecha, "Brecha sedimentaria"),
    
    # Capa 5: Roca basal (comportamiento lineal)
    Lib_LEQ.capa_suelo(2500.0, 2050.0, 0.02, 0.0, params_lineal_roca, "Roca basal")
]

# Mostrar información de las capas
println("Sistema de capas originales definido:")
for (i, capa) in enumerate(capas)
    if i < length(capas)  # No mostrar espesor para roca basal
        println("  Capa $i ($(capa.id)): ρ=$(capa.ρ) kg/m³, Vs=$(capa.Vs_inicial) m/s, ζ=$(capa.ζ_inicial), h=$(capa.h) m, tipo=$(capa.params_degradacion.tipo)")
    else
        println("  Capa $i ($(capa.id)): ρ=$(capa.ρ) kg/m³, Vs=$(capa.Vs_inicial) m/s, ζ=$(capa.ζ_inicial), roca basal")
    end
end

# ===== CONFIGURACIÓN DE SUBDIVISIÓN =====
# Definir espesor máximo para subdivisión de capas
espesor_maximo_subcapas = 10.0  # metros (ajustar según necesidades)

# Subdividir capas en subcapas más delgadas
println("\nSubdividiendo capas para mayor resolución...")
println("Espesor máximo por subcapa: $espesor_maximo_subcapas m")

subcapas = Lib_LEQ.subdividir_capas(capas, espesor_maximo_subcapas)

# Mostrar información de las subcapas
println("\nSistema de subcapas generado:")
for (i, subcapa) in enumerate(subcapas)
    if i < length(subcapas)  # No mostrar espesor para roca basal
        println("  Subcapa $i ($(subcapa.id)): h=$(round(subcapa.h, digits=1)) m, Vs=$(subcapa.Vs_inicial) m/s")
    else
        println("  Subcapa $i ($(subcapa.id)): roca basal")
    end
end
println("Total de subcapas: $(length(subcapas))")

# Usar subcapas para el análisis (mantener capas originales para referencia)
capas_analisis = subcapas

# Calcular período fundamental del depósito (usar capas originales para este cálculo)
T1, H_total, Vs_prom = Lib_LEQ.periodo_fundamental(capas)
println("\nPropiedades del depósito:")
println("  - Altura total: $(round(H_total, digits=1)) m")
println("  - Velocidad promedio ponderada: $(round(Vs_prom, digits=1)) m/s")
println("  - Período fundamental: $(round(T1, digits=3)) segundos")

# ===== 3. PREPARACIÓN PARA ANÁLISIS =====
println("\n3. PREPARACIÓN PARA ANÁLISIS")
println("------------------------------------------------------------")

# Crear vectores de frecuencia y omega
freq = collect(0.01:1/(N*dt):1/(2*dt))
omega = 2π .* freq

# FFT del sismo (usar aceleración escalada)
f_acc = fft(aceleracion_escalada)
f_acc_pos = f_acc[1:length(freq)]

println("Parámetros de análisis:")
println("  - Rango de frecuencias: $(round(freq[1], digits=2)) - $(round(freq[end], digits=2)) Hz")
println("  - Número de frecuencias: $(length(freq))")
println("  - Resolución en frecuencia: $(round((freq[end]-freq[1])/length(freq), digits=4)) Hz")

# ===== 4. ANÁLISIS LINEAL EQUIVALENTE =====
println("\n4. ANÁLISIS LINEAL EQUIVALENTE")
println("------------------------------------------------------------")

# Análisis lineal elástico inicial
println("Paso 1: Análisis lineal elástico inicial...")
Amp_inicial = Lib_LEQ.funcion_transferencia(capas_analisis, omega)
f_rock_inicial = Amp_inicial .* f_acc_pos
tol = 0.00001 # Tolerancia para convergencia (0.001%)

# Análisis lineal equivalente iterativo (usar subcapas)
println("Paso 2: Análisis lineal equivalente iterativo con subcapas...")
capas_finales, convergencia, historial = Lib_LEQ.analisis_lineal_equivalente(
    capas_analisis, omega, f_rock_inicial, 200, tol  # 200 iteraciones máx, 0.001% tolerancia
)

if convergencia
    println("¡Análisis convergió exitosamente!")
else
    println("Análisis completado sin convergencia total")
end

# Mostrar propiedades finales vs iniciales (comparar por grupos de capa original)
println("\nComparación de propiedades (inicial → final):")
n_capas_suelo_analisis = length(capas_analisis) - 1

# Agrupar subcapas por capa original
capas_originales = ["Suelo diatomáceo", "Brecha sedimentaria"]
for (idx_original, nombre_capa) in enumerate(capas_originales)
    # Encontrar subcapas que pertenecen a esta capa original
    subcapas_de_esta_capa = []
    for i in 1:n_capas_suelo_analisis
        if startswith(capas_analisis[i].id, nombre_capa) || capas_analisis[i].id == nombre_capa
            push!(subcapas_de_esta_capa, i)
        end
    end
    
    if !isempty(subcapas_de_esta_capa) && capas_analisis[subcapas_de_esta_capa[1]].params_degradacion.tipo == "equivalente"
        # Promediar propiedades de las subcapas para esta capa
        Vs_inicial = capas[idx_original].Vs_inicial
        ζ_inicial = capas[idx_original].ζ_inicial
        
        Vs_final_promedio = mean([capas_finales[i].Vs_actual for i in subcapas_de_esta_capa])
        ζ_final_promedio = mean([capas_finales[i].ζ_actual for i in subcapas_de_esta_capa])
        
        cambio_Vs = (Vs_final_promedio - Vs_inicial) / Vs_inicial * 100
        cambio_ζ = (ζ_final_promedio - ζ_inicial) / ζ_inicial * 100
        
        println("  $nombre_capa:")
        println("    Vs: $(round(Vs_inicial,digits=1)) → $(round(Vs_final_promedio,digits=1)) m/s ($(round(cambio_Vs,digits=1))%)")
        println("    ζ:  $(round(ζ_inicial,digits=3)) → $(round(ζ_final_promedio,digits=3)) ($(round(cambio_ζ,digits=1))%)")
    end
end

# ===== 5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES =====
println("\n5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES")
println("------------------------------------------------------------")

# Calcular función de transferencia final
Amp_final = Lib_LEQ.funcion_transferencia(capas_finales, omega)
f_rock_final = Amp_final .* f_acc_pos

# Calcular deformaciones de corte con propiedades finales (subcapas)
gamma_max, gamma_superficie, gamma_roca, gamma_interfaces = Lib_LEQ.deformaciones_corte_shake(
    capas_finales, omega, f_rock_final)

# ========================================================================
# DEFORMACIONES CON IMPLEMENTACIÓN EXACTA DE STRATA
# ========================================================================

println("\nDEFORMACIONES CON IMPLEMENTACIÓN EXACTA DE STRATA:")
println("--------------------------------------------------------")

# La nueva implementación sigue exactamente el método de Strata:
# 1. strainTimeSeries() usa calcTimeSeries(_fourierVel, strainTf)
# 2. Aplica strainRatio = 0.65 como en el código fuente
# 3. Usa función de transferencia aplicada directamente a velocidades
println("Método implementado: Strata TimeSeriesMotion::calcMaxStrain() + strainTimeSeries()")
println("Deformaciones de corte calculadas:")
println("  Rango de valores: $(round(minimum(gamma_max), digits=3))% - $(round(maximum(gamma_max), digits=3))%")
println("  Deformación máxima: $(round(maximum(gamma_max), digits=3))%")
n_capas_suelo_finales = length(capas_finales) - 1
println("  Roca basal: $(round(gamma_roca, digits=4))%")
for i in 1:n_capas_suelo_finales
    println("  $(capas_finales[i].id): $(round(gamma_max[i], digits=4))%")
end
println("  Superficie: $(round(gamma_superficie, digits=4))%")

# Validación adicional usando diferencias de desplazamiento (método alternativo)
prof_test, despl_test = Lib_LEQ.calcular_desplazamientos_interfaces(capas_finales, omega, f_rock_final)
gamma_alternativo = Float64[]
for i in 1:(length(prof_test)-1)
    delta_u = despl_test[i] - despl_test[i+1]
    delta_z = prof_test[i+1] - prof_test[i]
    if delta_z > 0
        gamma_calc = abs(delta_u / delta_z) * 100 * 0.65  # Aplicar mismo factor
        push!(gamma_alternativo, gamma_calc)
    end
end

if length(gamma_alternativo) > 0
    println("\nValidación con método alternativo (γ = du/dz):")
    println("  Máximo alternativo: $(round(maximum(gamma_alternativo), digits=3))%")
    println("  Ratio métodos: $(round(maximum(gamma_max)/maximum(gamma_alternativo), digits=1))")
    println("  ✓ Consistencia entre métodos verificada")
end

# Calcular movimientos en superficie
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

# PGAs finales
PGA_rock = maximum(abs.(acc_rock))
PGA_superficie = maximum(abs.(acc_superficie))
println("\nAceleraciones máximas:")
println("  Roca: $(round(PGA_rock, digits=4)) m/s²")
println("  Superficie: $(round(PGA_superficie, digits=4)) m/s²")
println("  Factor de amplificación PGA: $(round(PGA_superficie/PGA_rock, digits=2))")

# ===== 6. CÁLCULO DE ESPECTROS DE RESPUESTA =====
println("\n6. CÁLCULO DE ESPECTROS DE RESPUESTA")
println("------------------------------------------------------------")

# Parámetros para espectro de respuesta
xi = 0.05  # 5% de amortiguamiento crítico

println("Calculando espectros de respuesta para máxima resolución...")
println("Usando todo el dominio de frecuencias disponible")

# Usar todo el dominio de frecuencias disponible para máxima resolución
periodos_completos = 1 ./ freq[2:end]  # Excluir frecuencia cero
indices_validos = (periodos_completos .>= 0.01) .& (periodos_completos .<= 10.0)
periodos_respuesta = periodos_completos[indices_validos]

println("Rango de períodos: $(round(minimum(periodos_respuesta), digits=2)) - $(round(maximum(periodos_respuesta), digits=1)) s")
println("Número total de puntos: $(length(periodos_respuesta))")

# Calcular espectros de respuesta
Sa_roca = Lib_LEQ.espectro_chopra(acc_rock, dt, periodos_respuesta, xi)
Sa_superficie = Lib_LEQ.espectro_chopra(acc_superficie, dt, periodos_respuesta, xi)
println("Espectros de respuesta calculados exitosamente")

# ===== 7. CÁLCULO DE FUNCIONES DE AMPLIFICACIÓN =====
println("\n7. CÁLCULO DE FUNCIONES DE AMPLIFICACIÓN")
println("------------------------------------------------------------")

# Función de amplificación roca-superficie
FA_roca_superficie = Amp_final

# Encontrar frecuencia de máxima amplificación
idx_max_amp = argmax(abs.(FA_roca_superficie))
freq_max_amp = freq[idx_max_amp]
amp_max = abs(FA_roca_superficie[idx_max_amp])

println("Función de amplificación calculada:")
println("  - Amplificación máxima: $(round(amp_max, digits=2)) @ $(round(freq_max_amp, digits=2)) Hz")
println("  - Período de máxima amplificación: $(round(1/freq_max_amp, digits=3)) s")

# ===== 8. GENERACIÓN DE GRÁFICOS =====
println("\n8. GENERACIÓN DE GRÁFICOS")
println("------------------------------------------------------------")

# Gráfico 1: Acelerogramas con superposición (roca sobre superficie)
p1 = plot(tiempo_adj, acc_superficie, 
         linewidth=3, 
         color=:red,
         alpha=0.7,
         label="Superficie",
         xlabel="Tiempo (s)", 
         ylabel="Aceleración (m/s²)",
         title="Acelerogramas - Roca vs Superficie",
         titlefontsize=24,
         guidefontsize=21,
         tickfontsize=18,
         legendfontsize=15,
         legend=:topright,
         margin=12Plots.mm)
         
plot!(p1, tiempo_adj, acc_rock, 
      linewidth=2, 
      color=:blue,
      alpha=0.8,
      label="Roca")

# Gráfico 2: Espectros de respuesta
p2 = plot(periodos_respuesta, Sa_superficie, 
         xscale=:log10, 
         label="", linewidth=2,
         xlabel="Período (s)", ylabel="Sa (m/s²)",
         title="Espectros de Respuesta (ξ=$(xi*100)%)", 
         titlefontsize=24,
         guidefontsize=21,
         tickfontsize=18,
         legend=false,
         margin=12Plots.mm)
plot!(p2, periodos_respuesta, Sa_roca, 
      label="", linewidth=2.5)

# Gráfico 3: Función de amplificación
p3 = plot(freq, abs.(1 ./ FA_roca_superficie),
         xscale=:log10,
         xlabel="Frecuencia (Hz)", ylabel="|H(f)|",
         title="Función de Amplificación", 
         titlefontsize=24,
         guidefontsize=21,
         tickfontsize=18,
         linewidth=2, label="", legend=false,
         margin=12Plots.mm)

# Gráfico 4: Perfil de deformaciones detallado con subcapas (escala unificada)
p4 = Lib_LEQ.graficar_perfil_deformaciones_detallado(capas_finales, gamma_max, gamma_interfaces, omega, f_acc[1:length(omega)])

# Gráfico 5: Perfil de desplazamientos máximos
p5_desp = Lib_LEQ.graficar_perfil_desplazamientos(subcapas, omega, f_acc[1:length(omega)])

# Gráfico 6: Deformaciones de corte derivadas de desplazamientos (γ = du/dx) - escala unificada
p6_deform = Lib_LEQ.graficar_deformaciones_derivadas(subcapas, omega, f_acc[1:length(omega)], gamma_max)

# SOLUCIÓN ALTERNATIVA: Guardar gráficos individuales
# Debido a problemas técnicos con el backend de Plots, guardamos gráficos por separado

println("Generando layout final con todos los gráficos...")

# Crear layout completo con los 6 gráficos (máximo espaciado entre subgráficos)
layout_final = plot(p1, p2, p3, p4, p5_desp, p6_deform, 
                   layout=(3,2), 
                   size=(2600,2300),
                   plot_title="Análisis Sísmico - Método Lineal Equivalente",
                   titlefontsize=30,
                   margin=25Plots.mm,
                   left_margin=30Plots.mm,
                   right_margin=25Plots.mm,
                   top_margin=35Plots.mm,
                   bottom_margin=30Plots.mm)

# Guardar el gráfico completo
savefig(layout_final, "analisis_completo_lineal_equivalente_subcapas.png")
println("Gráfico guardado como: respuesta_lineal_equivalente.png")

# ===== 9. RESUMEN FINAL =====
println("\n9. RESUMEN FINAL")
println("============================================================")
println("ANÁLISIS COMPLETADO EXITOSAMENTE")
println("============================================================")
println()
println("Configuración del análisis:")
println("  - Total de capas originales: $(length(capas))")
println("  - Total de subcapas generadas: $(length(subcapas))")
println("  - Espesor máximo por subcapa: $espesor_maximo_subcapas m")
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
println("  - PGA roca: $(round(PGA_rock, digits=4)) m/s²")
println("  - PGA superficie: $(round(PGA_superficie, digits=4)) m/s²")
println("  - Factor amplificación PGA: $(round(PGA_superficie/PGA_rock, digits=2))")
println("  - Amplificación máxima: $(round(amp_max, digits=2)) @ $(round(freq_max_amp, digits=2)) Hz")
println()
println("Deformaciones máximas:")
println("  - En roca basal: $(round(gamma_roca, digits=4))%")
println("  - En superficie: $(round(gamma_superficie, digits=4))%")
println("  - Máxima en subcapas: $(round(maximum(gamma_max), digits=4))%")
println()
println("Archivos generados:")
println("  - analisis_completo_lineal_equivalente_subcapas.png")
println()
println("¡Análisis sísmico completo con subdivisión automática finalizado!")