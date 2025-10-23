"""
Lib_Resumen.jl - Librería para generación de resumen final
Etapa 9: Resumen final del análisis con estadísticas optimizadas
"""

module Lib_Resumen

using Statistics, Printf
using ..Lib_TransferenciaCore: capa_suelo

export generar_resumen_final

"""
    generar_resumen_final(capas_originales::Vector{capa_suelo}, subcapas::Vector{capa_suelo}, 
                         convergencia::Bool, historial::Dict, T1::Float64, H_total::Float64, 
                         T_p::Float64, resultados_movimientos::Dict, resultados_espectros::Dict, 
                         gamma_max::Vector{Float64}, archivo_salida::String)

Función principal para la generación del resumen final.
Implementa la Etapa 9 del análisis con análisis estadístico optimizado.
"""
function generar_resumen_final(capas_originales::Vector{capa_suelo}, subcapas::Vector{capa_suelo}, 
                              convergencia::Bool, historial::Dict, T1::Float64, H_total::Float64, 
                              T_p::Float64, resultados_movimientos::Dict, resultados_espectros::Dict, 
                              gamma_max::Vector{Float64}, archivo_salida::String, espesor_maximo::Float64)
    println("\n9. RESUMEN FINAL")
    println("============================================================")
    
    # Función anidada para análisis estadístico de convergencia
    function analizar_convergencia_estadisticas(hist_dict, convergio)
        if !convergio
            return "❌ NO CONVERGIÓ", 0, 0.0, 0.0
        end
        
        n_iter = length(hist_dict["iteracion"])
        error_inicial = hist_dict["error"][1]
        error_final = hist_dict["error"][end]
        reduccion = ((error_inicial - error_final) / error_inicial) * 100
        
        estado = n_iter <= 50 ? "✅ RÁPIDA" : n_iter <= 100 ? "✅ NORMAL" : "✅ LENTA"
        
        return estado, n_iter, error_inicial, reduccion
    end
    
    # Función anidada para análisis de propiedades del suelo
    function analizar_propiedades_suelo(capas_orig, subcapas_finales)
        # Cambios en velocidades
        cambios_vs = Float64[]
        cambios_zeta = Float64[]
        
        for (i, capa_orig) in enumerate(capas_orig[1:end-1])  # Excluir roca
            if capa_orig.params_degradacion.tipo == "equivalente"
                # Encontrar subcapas correspondientes
                subcapas_capa = [sc for sc in subcapas_finales if startswith(sc.id, capa_orig.id)]
                
                if !isempty(subcapas_capa)
                    vs_inicial = capa_orig.Vs_inicial
                    vs_final = mean([sc.Vs_actual for sc in subcapas_capa])
                    zeta_inicial = capa_orig.ζ_inicial
                    zeta_final = mean([sc.ζ_actual for sc in subcapas_capa])
                    
                    push!(cambios_vs, (vs_final - vs_inicial) / vs_inicial * 100)
                    push!(cambios_zeta, (zeta_final - zeta_inicial) / zeta_inicial * 100)
                end
            end
        end
        
        return cambios_vs, cambios_zeta
    end
    
    # Función anidada para análisis de respuesta sísmica
    function analizar_respuesta_sismica(res_movimientos, res_espectros)
        # PGAs y factores de amplificación
        pga_roca = res_movimientos["PGA_rock"]
        pga_superficie = res_movimientos["PGA_superficie"]
        factor_amp_pga = pga_superficie / pga_roca
        
        # Espectros máximos
        sa_max_roca = maximum(res_espectros["Sa_rock"])
        sa_max_superficie = maximum(res_espectros["Sa_superficie"])
        factor_amp_sa = sa_max_superficie / sa_max_roca
        
        # Períodos de máxima respuesta
        idx_max_sa = argmax(res_espectros["Sa_superficie"])
        periodo_max_sa = res_espectros["periodos"][idx_max_sa]
        
        return pga_roca, pga_superficie, factor_amp_pga, sa_max_roca, sa_max_superficie, factor_amp_sa, periodo_max_sa
    end
    
    # Función anidada para análisis de deformaciones
    function analizar_deformaciones_detallado(gamma_vector, subcapas_vec)
        gamma_max_global = maximum(gamma_vector)
        gamma_promedio = mean(gamma_vector)
        gamma_mediana = median(gamma_vector)
        
        # Ubicación de máxima deformación
        idx_max = argmax(gamma_vector)
        capa_max_def = idx_max <= length(subcapas_vec) ? subcapas_vec[idx_max].id : "N/A"
        
        # Profundidad de máxima deformación
        prof_max = sum([subcapas_vec[i].h for i in 1:min(idx_max-1, length(subcapas_vec)-1)])
        
        # Distribución por rangos
        rango_bajo = count(x -> x < 0.1, gamma_vector) / length(gamma_vector) * 100
        rango_medio = count(x -> 0.1 <= x < 0.5, gamma_vector) / length(gamma_vector) * 100
        rango_alto = count(x -> x >= 0.5, gamma_vector) / length(gamma_vector) * 100
        
        return gamma_max_global, gamma_promedio, gamma_mediana, capa_max_def, prof_max, rango_bajo, rango_medio, rango_alto
    end
    
    # Ejecutar análisis optimizado
    estado_conv, n_iter, error_ini, reduccion = analizar_convergencia_estadisticas(historial, convergencia)
    cambios_vs, cambios_zeta = analizar_propiedades_suelo(capas_originales, subcapas)
    pga_r, pga_s, f_amp_pga, sa_max_r, sa_max_s, f_amp_sa, T_max_sa = analizar_respuesta_sismica(resultados_movimientos, resultados_espectros)
    gamma_max_g, gamma_prom, gamma_med, capa_max, prof_max, r_bajo, r_medio, r_alto = analizar_deformaciones_detallado(gamma_max, subcapas)
    
    # REPORTE FINAL OPTIMIZADO
    println("ANÁLISIS SÍSMICO COMPLETADO - MÉTODO LINEAL EQUIVALENTE")
    println("============================================================")
    
    println("\n📊 ESTADÍSTICAS DE CONVERGENCIA:")
    println("   Estado: $estado_conv")
    if convergencia
        println("   Iteraciones: $n_iter")
        println("   Error inicial: $(round(error_ini, digits=2))%")
        println("   Reducción de error: $(round(reduccion, digits=1))%")
    end
    
    println("\n🏗️ MODELO DE SUELO:")
    println("   Altura total del depósito: $(round(H_total, digits=1)) m")
    println("   Período fundamental: $(round(T1, digits=3)) s")
    println("   Período predominante del sismo: $(round(T_p, digits=3)) s")
    println("   Número de subcapas: $(length(subcapas)-1)")
    println("   Espesor máximo por subcapa: $(espesor_maximo) m")
    
    if !isempty(cambios_vs)
        println("\n🔄 CAMBIOS EN PROPIEDADES:")
        println("   Velocidad Vs: $(round(mean(cambios_vs), digits=1))% (rango: $(round(minimum(cambios_vs), digits=1))% a $(round(maximum(cambios_vs), digits=1))%)")
        if !isempty(cambios_zeta)
            println("   Amortiguamiento ζ: $(round(mean(cambios_zeta), digits=1))% (rango: $(round(minimum(cambios_zeta), digits=1))% a $(round(maximum(cambios_zeta), digits=1))%)")
        end
    end
    
    println("\n🌊 RESPUESTA SÍSMICA:")
    println("   PGA roca: $(round(pga_r, digits=4)) m/s²")
    println("   PGA superficie: $(round(pga_s, digits=4)) m/s²")
    println("   Factor amplificación PGA: $(round(f_amp_pga, digits=2))")
    println("   Sa máx roca: $(round(sa_max_r, digits=2)) m/s²")
    println("   Sa máx superficie: $(round(sa_max_s, digits=2)) m/s²")
    println("   Factor amplificación Sa: $(round(f_amp_sa, digits=2))")
    println("   Período de máxima respuesta: $(round(T_max_sa, digits=3)) s")
    
    println("\n📏 DEFORMACIONES DE CORTE:")
    println("   Deformación máxima: $(round(gamma_max_g, digits=3))% en $capa_max")
    println("   Profundidad de máx. deformación: $(round(prof_max, digits=1)) m")
    println("   Deformación promedio: $(round(gamma_prom, digits=3))%")
    println("   Deformación mediana: $(round(gamma_med, digits=3))%")
    println("   Distribución: Baja (<0.1%): $(round(r_bajo, digits=1))%, Media (0.1-0.5%): $(round(r_medio, digits=1))%, Alta (≥0.5%): $(round(r_alto, digits=1))%")
    
    println("\n💾 ARCHIVOS GENERADOS:")
    println("   Gráficos: $archivo_salida")
    
    # Evaluación de calidad del análisis
    println("\n🎯 EVALUACIÓN DE CALIDAD:")
    
    calidad_convergencia = convergencia && n_iter <= 100 ? "EXCELENTE" : convergencia ? "BUENA" : "DEFICIENTE"
    println("   Convergencia: $calidad_convergencia")
    
    calidad_amplificacion = 1.2 <= f_amp_pga <= 5.0 ? "REALISTA" : f_amp_pga > 5.0 ? "ALTA" : "BAJA"
    println("   Amplificación: $calidad_amplificacion (Factor PGA: $(round(f_amp_pga, digits=2)))")
    
    calidad_deformaciones = gamma_max_g < 2.0 ? "ACEPTABLE" : gamma_max_g < 5.0 ? "ALTA" : "CRÍTICA"
    println("   Deformaciones: $calidad_deformaciones (Máx: $(round(gamma_max_g, digits=3))%)")
    
    println("\n" * "="^60)
    println("✅ ANÁLISIS COMPLETADO EXITOSAMENTE")
    println("="^60)
    return nothing
end

end # module