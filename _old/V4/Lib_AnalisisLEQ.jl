"""
Lib_AnalisisLEQ.jl - Librería para análisis lineal equivalente
Etapa 4: Análisis lineal equivalente con convergencia avanzada
"""

module Lib_AnalisisLEQ

using LinearAlgebra, Printf, Statistics
using ..Lib_TransferenciaCore: capa_suelo, curva_degradacion, curva_amortiguamiento,
                                matriz_T, funcion_transferencia, actualizar_propiedades_dinamicas!,
                                deformaciones_corte_shake

export ejecutar_analisis_leq, manejar_no_convergencia

"""
    ejecutar_analisis_leq(subcapas::Vector{capa_suelo}, omega::Vector{Float64}, 
                         f_acc_pos::Vector{ComplexF64}, max_iter::Int, tolerancia::Float64)

Función principal para el análisis lineal equivalente.
Implementa la Etapa 4 del análisis con algoritmo de convergencia mejorado.
"""
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
        mostrar_comparacion_propiedades(subcapas, capas_finales)
    else
        # 🚨 MANEJO CRÍTICO DE NO CONVERGENCIA
        println("❌ FALLO CRÍTICO: El análisis no convergió")
        manejar_no_convergencia(historial, tolerancia, max_iter)
        Amp_final = Amp_inicial  # Código de respaldo
    end
    
    return capas_finales, convergencia, historial, Amp_inicial, Amp_final
end

"""
    analisis_lineal_equivalente(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                               aceleracion_base_fft::Vector{ComplexF64}, 
                               max_iter::Int=10, tol::Float64=0.05)

Algoritmo de convergencia lineal equivalente con relajación adaptativa.
"""
function analisis_lineal_equivalente(capas::Vector{capa_suelo}, omega::Vector{Float64}, 
                                   aceleracion_base_fft::Vector{ComplexF64}, 
                                   max_iter::Int=10, tol::Float64=0.05)
    
    # Hacer copia de las capas para modificar
    capas_trabajo = deepcopy(capas)
    n_capas = length(capas_trabajo) - 1  # Excluir roca
    
    # Historial para tracking
    historial = Dict{String, Vector{Any}}(
        "iteracion" => Int[],
        "error" => Float64[],
        "factor_relajacion" => Float64[],
        "Vs" => Vector{Float64}[],
        "zeta" => Vector{Float64}[],
        "gamma" => Vector{Float64}[]
    )
    
    # Variables de control del algoritmo mejorado
    factor_relajacion = 0.5  # Factor inicial conservador
    contador_estabilidad = 0
    errors_recientes = Float64[]
    error_anterior = Inf
    
    # ===== ALGORITMO MEJORADO DE CONVERGENCIA =====
    # Obtener nombres de las capas originales para las columnas
    capas_nombres = []
    if length(capas) >= 1 && !isempty(capas[1].id)
        push!(capas_nombres, capas[1].id)
    else
        push!(capas_nombres, "Capa_1")
    end
    
    if length(capas) >= 2 && !isempty(capas[2].id)
        push!(capas_nombres, capas[2].id)
    else
        push!(capas_nombres, "Capa_2")
    end
    
    # Truncar nombres largos para que quepan en la tabla
    nombre_capa1 = length(capas_nombres[1]) > 8 ? capas_nombres[1][1:5] * "..." : capas_nombres[1]
    nombre_capa2 = length(capas_nombres[2]) > 8 ? capas_nombres[2][1:5] * "..." : capas_nombres[2]
    
    println("\nTabla de convergencia del análisis lineal equivalente (ALGORITMO MEJORADO):")
    println(repeat("=", 120))
    println("| Iter |   Error  | Factor_Relaj |   Vs_$nombre_capa1   |  Vs_$nombre_capa2  | γ_max_$nombre_capa1 | γ_max_$nombre_capa2 |   ζ_$nombre_capa1   | Estado |")
    println("|      |    (%)   |      (-)     |     (m/s)     |   (m/s)   |     (%)    |      (%)     |      (-)     |        |")
    println(repeat("=", 120))
    
    convergencia = false
    
    for iter in 1:max_iter
        # Calcular deformaciones actuales
        gamma_max, _, _, _ = deformaciones_corte_shake(capas_trabajo, omega, aceleracion_base_fft)
        
        # Actualizar propiedades solo para capas equivalentes
        Vs_nuevos = Vector{Float64}(undef, n_capas)
        zeta_nuevos = Vector{Float64}(undef, n_capas)
        
        for i in 1:n_capas
            if capas_trabajo[i].params_degradacion.tipo == "equivalente"
                # Usar γ_efectiva = 0.65 * γ_max
                gamma_efectiva = gamma_max[i] * 0.65 / 100  # Convertir a fracción
                
                # Calcular nuevas propiedades
                G_Go = curva_degradacion(gamma_efectiva, capas_trabajo[i].params_degradacion)
                
                G_inicial = capas_trabajo[i].ρ * capas_trabajo[i].Vs_inicial^2
                G_nuevo = G_inicial * G_Go
                Vs_nuevo = sqrt(G_nuevo / capas_trabajo[i].ρ)
                
                zeta_nuevo = curva_amortiguamiento(gamma_efectiva, capas_trabajo[i].params_degradacion)
                
                Vs_nuevos[i] = Vs_nuevo
                zeta_nuevos[i] = max(zeta_nuevo, 0.02)  # Mínimo 2%
            else
                # Capas lineales mantienen sus propiedades
                Vs_nuevos[i] = capas_trabajo[i].Vs_actual
                zeta_nuevos[i] = capas_trabajo[i].ζ_actual
            end
        end
        
        # Aplicar relajación adaptativa
        for i in 1:n_capas
            if capas_trabajo[i].params_degradacion.tipo == "equivalente"
                capas_trabajo[i].Vs_actual = (1 - factor_relajacion) * capas_trabajo[i].Vs_actual + 
                                           factor_relajacion * Vs_nuevos[i]
                capas_trabajo[i].ζ_actual = (1 - factor_relajacion) * capas_trabajo[i].ζ_actual + 
                                          factor_relajacion * zeta_nuevos[i]
                capas_trabajo[i].G_actual = capas_trabajo[i].ρ * capas_trabajo[i].Vs_actual^2
            end
        end
        
        # Calcular errores
        error_vs = maximum([abs(capas_trabajo[i].Vs_actual - Vs_nuevos[i]) / Vs_nuevos[i] 
                           for i in 1:n_capas if capas_trabajo[i].params_degradacion.tipo == "equivalente"]) * 100
        
        error_zeta = maximum([abs(capas_trabajo[i].ζ_actual - zeta_nuevos[i]) / zeta_nuevos[i] 
                             for i in 1:n_capas if capas_trabajo[i].params_degradacion.tipo == "equivalente"]) * 100
        
        error_combinado = sqrt(error_vs^2 + error_zeta^2)
        
        # Guardar historial
        push!(historial["iteracion"], iter)
        push!(historial["error"], error_combinado)
        push!(historial["factor_relajacion"], factor_relajacion)
        push!(historial["Vs"], [capa.Vs_actual for capa in capas_trabajo[1:n_capas]])
        push!(historial["zeta"], [capa.ζ_actual for capa in capas_trabajo[1:n_capas]])
        push!(historial["gamma"], copy(gamma_max))
        
        # Métricas para reporte - dinámicas basadas en las primeras capas
        # Simplificamos mostrando las primeras 2 subcapas como representativas
        gamma_capa1_max = gamma_max[1]
        gamma_capa2_max = length(gamma_max) >= 2 ? gamma_max[2] : 0.0
        Vs_capa1 = capas_trabajo[1].Vs_actual
        Vs_capa2 = length(capas_trabajo) >= 2 ? capas_trabajo[2].Vs_actual : 0.0
        zeta_capa1 = capas_trabajo[1].ζ_actual
        
        # === ALGORITMO DE CONTROL ADAPTATIVO ===
        estado = "Normal"
        
        # Detectar estabilidad/oscilaciones
        push!(errors_recientes, error_combinado)
        if length(errors_recientes) > 5
            errors_recientes = errors_recientes[end-4:end]
        end
        
        # Detectar patrones de estabilidad
        if length(errors_recientes) >= 4
            variacion = std(errors_recientes)
            if variacion < 0.05 * mean(errors_recientes)  # Baja variación
                contador_estabilidad += 1
                estado = "Estable"
                
                # Incrementar factor de relajación gradualmente si es estable
                if contador_estabilidad >= 2 && factor_relajacion < 0.8
                    factor_relajacion = min(0.8, factor_relajacion + 0.06)
                end
            else
                contador_estabilidad = 0
            end
        end
        
        # Imprimir fila de la tabla
        @printf("|  %3d | %7.2f%% |     %7.2f |     %7.1f |   %7.1f |    %6.3f |      %6.3f |      %6.3f | %6s |\n",
                iter, error_combinado, factor_relajacion, Vs_capa1, Vs_capa2, 
                gamma_capa1_max, gamma_capa2_max, zeta_capa1, estado)
        
        # Verificar convergencia con criterios múltiples
        if error_vs < tol && error_zeta < tol && error_combinado < tol
            println(repeat("-", 120))
            println("CONVERGENCIA ALCANZADA en iteracion ", iter, "!")
            println("  - Error Vs: ", round(error_vs, digits=3), "% < ", tol, "%")
            println("  - Error zeta: ", round(error_zeta, digits=1), "% < ", tol, "%")
            println("  - Error combinado: ", round(error_combinado, digits=3), "% < ", tol, "%")
            println(repeat("=", 120))
            convergencia = true
            break
        end
        
        error_anterior = error_combinado
    end
    
    println(repeat("=", 120))
    
    # Calcular funciones de transferencia inicial y final
    Amp_inicial = funcion_transferencia(capas, omega)
    if convergencia
        Amp_final = funcion_transferencia(capas_trabajo, omega)
    else
        Amp_final = Amp_inicial
    end
    
    return capas_trabajo, convergencia, historial, Amp_inicial, Amp_final
end

"""
    mostrar_comparacion_propiedades(subcapas_inicial, capas_finales)

Muestra la comparación de propiedades entre estado inicial y final.
"""
function mostrar_comparacion_propiedades(subcapas_inicial, capas_finales)
    println("\nComparación de propiedades (inicial → final):")
    for i in 1:min(length(subcapas_inicial)-1, 22)  # Solo mostrar primeras capas representativas
        capa_inicial = subcapas_inicial[i]
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
end

"""
    manejar_no_convergencia(historial::Dict, tolerancia::Float64, max_iter::Int)

Maneja casos de no convergencia con análisis detallado y recomendaciones.
"""
function manejar_no_convergencia(historial::Dict, tolerancia::Float64, max_iter::Int)
    println("\n" * "="^80)
    println("🚨 ANÁLISIS DE NO CONVERGENCIA - DIAGNÓSTICO DETALLADO")
    println("="^80)
    
    n_iter = length(historial["iteracion"])
    
    if n_iter == 0
        println("❌ ERROR CRÍTICO: No se completó ninguna iteración")
        println("   Posibles causas:")
        println("   • Error en configuración de parámetros")
        println("   • Problema con función de transferencia")
        println("   • Datos de entrada inválidos")
        error("Análisis terminado por falta de datos de iteración")
    end
    
    error_inicial = historial["error"][1]
    error_final = historial["error"][end]
    reduccion_error = ((error_inicial - error_final) / error_inicial) * 100
    
    println("📊 ESTADÍSTICAS DE CONVERGENCIA:")
    println("   • Iteraciones completadas: $n_iter de $max_iter")
    println("   • Tolerancia objetivo: $tolerancia%")
    println("   • Error inicial: $(round(error_inicial, digits=2))%")
    println("   • Error final: $(round(error_final, digits=2))%") 
    println("   • Reducción del error: $(round(reduccion_error, digits=1))%")
    
    # Análisis de tendencia
    if n_iter >= 5
        errors_ultimos_5 = historial["error"][end-4:end]
        tendencia = mean(diff(errors_ultimos_5))
        
        if tendencia > 0.01
            println("   ⚠️  TENDENCIA: Divergencia detectada (+$(round(tendencia, digits=3))%/iter)")
        elseif abs(tendencia) < 0.001
            println("   ⚠️  TENDENCIA: Estancamiento detectado (±$(round(tendencia, digits=4))%/iter)")
        else
            println("   ✓  TENDENCIA: Convergencia lenta (-$(round(abs(tendencia), digits=3))%/iter)")
        end
    end
    
    println("\n🔧 RECOMENDACIONES DE AJUSTE:")
    
    if error_final > 10.0
        println("   • ERROR ALTO (>10%): Problema fundamental de configuración")
        println("     - Revisar parámetros de curvas de degradación")
        println("     - Verificar propiedades iniciales del suelo")
        println("     - Considerar subdivisión más fina de capas")
    elseif error_final > 5.0
        println("   • ERROR MODERADO (5-10%): Ajuste de algoritmo necesario")
        println("     - Aumentar número máximo de iteraciones a $(max_iter * 2)")
        println("     - Reducir factor de relajación inicial a 0.3")
        println("     - Incrementar tolerancia temporalmente a $(tolerancia * 2)%")
    else
        println("   • ERROR BAJO (<5%): Convergencia casi alcanzada")
        println("     - Aumentar iteraciones máximas a $(max_iter + 20)")
        println("     - Mantener parámetros actuales")
        println("     - Considerar usar resultado actual (precisión aceptable)")
    end
    
    println("\n" * "="^80)
    error("ANÁLISIS LINEAL EQUIVALENTE NO CONVERGIÓ - Aplicar recomendaciones y reintentar")
end

end # module