module Lib_Graficos

using Plots
using FFTW  
using Statistics
using ..Lib_TransferenciaCore: capa_suelo, funcion_transferencia

export generar_graficos_completos, generar_grafico_velocidades_corte, generar_grafico_amortiguamiento

function generar_graficos_completos(tiempo::Vector{Float64}, resultados_movimientos::Dict, 
                                   resultados_espectros::Dict, subcapas::Vector{capa_suelo}, 
                                   omega::Vector{Float64}, f_acc_pos::Vector{ComplexF64}, 
                                   gamma_max::Vector{Float64}, archivo_salida::String, 
                                   capas_finales::Union{Vector{capa_suelo}, Nothing}=nothing,
                                   capas_originales::Union{Vector{capa_suelo}, Nothing}=nothing)
    
    println("\n8. GENERACION DE GRAFICOS")
    println("------------------------------------------------------------")
    
    try
        gr()
        # Layout 3x2 exacto: 3 filas, 2 columnas
        layout_3x2 = @layout([
            [a; b] c
            d e  
            f g
        ])
        p = plot(layout=layout_3x2, size=(1400, 1200), dpi=150)
        
        # Posición 1 (1,1a): Acelerograma en roca
        println("  Generando acelerogramas...")
        
        try
            plot!(p[1], 
                  title="Acelerogramas Roca y Superficie",
                  xlabel="Tiempo (s)", 
                  ylabel="Aceleración (m/s²)", 
                  legend=:topright,
                  grid=true)
            
            # Primero graficar superficie (por debajo)
            if haskey(resultados_movimientos, "tiempo_adj") && haskey(resultados_movimientos, "acc_superficie")
                tiempo_adj = resultados_movimientos["tiempo_adj"]
                acc_superficie = resultados_movimientos["acc_superficie"]
                
                if length(tiempo_adj) > 0 && length(acc_superficie) > 0
                    plot!(p[1], tiempo_adj, acc_superficie, 
                          label="Superficie", 
                          color=:blue, 
                          lw=1.5,
                          alpha=0.7)
                end
            end
            
            # Luego graficar roca (por encima)
            if haskey(resultados_movimientos, "tiempo_adj") && haskey(resultados_movimientos, "acc_rock")
                tiempo_adj = resultados_movimientos["tiempo_adj"]
                acc_rock = resultados_movimientos["acc_rock"]
                
                if length(tiempo_adj) > 0 && length(acc_rock) > 0
                    plot!(p[1], tiempo_adj, acc_rock, 
                          label="Roca", 
                          color=:red, 
                          lw=1.5)
                end
            end
            
            # Si no hay datos
            if !haskey(resultados_movimientos, "tiempo_adj")
                plot!(p[1], [0, 10], [0, 0.1], label="No disponible")
            end
        catch e
            println("    Error en acelerogramas: ", e)
            plot!(p[1], [0, 10], [0, 0.1], 
                  title="Acelerogramas Roca y Superficie",
                  xlabel="Tiempo (s)", 
                  ylabel="Aceleración (m/s²)", 
                  label="Error")
        end
        
        # Posición 2 (1,1b): Acelerogramas Superficie
        try
            plot!(p[2], 
                  title="Acelerogramas Superficie",
                  xlabel="Tiempo (s)", 
                  ylabel="Aceleración (m/s²)", 
                  legend=:topright,
                  grid=true)
            
            if haskey(resultados_movimientos, "tiempo_adj") && haskey(resultados_movimientos, "acc_superficie")
                tiempo_adj = resultados_movimientos["tiempo_adj"]
                acc_superficie = resultados_movimientos["acc_superficie"]
                
                if length(tiempo_adj) > 0 && length(acc_superficie) > 0
                    plot!(p[2], tiempo_adj, acc_superficie, 
                          label="Superficie",
                          color=:blue,
                          lw=1.5)
                end
            end
            
            if haskey(resultados_movimientos, "tiempo_conv") && haskey(resultados_movimientos, "acc_superficie_convolucion")
                tiempo_conv = resultados_movimientos["tiempo_conv"]
                acc_reconstruido = resultados_movimientos["acc_superficie_convolucion"]
                
                if length(tiempo_conv) > 0 && length(acc_reconstruido) > 0
                    plot!(p[2], tiempo_conv, acc_reconstruido, 
                          label="Reconstruida",
                          color=:green,
                          lw=1.5,
                          linestyle=:dash)
                end
            end
            
            if !haskey(resultados_movimientos, "tiempo_adj")
                plot!(p[2], [0, 10], [0, 0.1], label="No disponible")
            end
        catch e
            println("    Error en acelerogramas superficie: ", e)
        end

        # Posición 3 (1,2): Espectro de Respuesta
        println("  Generando espectros de respuesta...")
        try
            if haskey(resultados_espectros, "periodos")
                periodos = resultados_espectros["periodos"]
                
                plot!(p[3], title="Espectro de Respuesta",
                      xlabel="Período (s)",
                      ylabel="Aceleración Espectral (m/s²)",
                      legend=:topright,
                      grid=true,
                      xscale=:log10,
                      xlims=(0.01, 10.0))
                
                if haskey(resultados_espectros, "Sa_superficie") && length(resultados_espectros["Sa_superficie"]) > 0
                    plot!(p[3], periodos, resultados_espectros["Sa_superficie"],
                          label="Superficie",
                          color=:blue,
                          lw=2)
                end
                
                if haskey(resultados_espectros, "Sa_roca") && length(resultados_espectros["Sa_roca"]) > 0
                    plot!(p[3], periodos, resultados_espectros["Sa_roca"],
                          label="Roca",
                          color=:red,
                          lw=2)
                end
                
                if haskey(resultados_espectros, "Sa_reconstruido") && length(resultados_espectros["Sa_reconstruido"]) > 0
                    plot!(p[3], periodos, resultados_espectros["Sa_reconstruido"],
                          label="Reconstruido",
                          color=:green,
                          lw=2,
                          linestyle=:dash)
                end
            else
                plot!(p[3], [0.01, 10], [0, 1],
                      title="Espectro de Respuesta",
                      xlabel="Período (s)",
                      ylabel="Aceleración Espectral (m/s²)",
                      xscale=:log10,
                      label="No disponible")
            end
        catch e
            println("    Error en espectros: ", e)
        end

        # Posición 4 (2,1): Función de Transferencia
        println("  Generando función de transferencia...")
        try
            if haskey(resultados_movimientos, "H_transferencia") && length(omega) > 0
                H_transferencia = resultados_movimientos["H_transferencia"]
                freqs = omega ./ (2 * π)
                H_transfer = abs.(H_transferencia)
                
                indices_validos = freqs .<= 100.0
                freqs_plot = freqs[indices_validos]
                H_plot = H_transfer[indices_validos]
                
                plot!(p[4], freqs_plot, H_plot,
                      title="Función de Transferencia",
                      xlabel="Frecuencia (Hz)",
                      ylabel="Amplificación",
                      legend=false,
                      grid=true,
                      lw=2,
                      color=:purple)
                
                xlims!(p[4], (0, 100))
            else
                plot!(p[4], [0, 100], [0, 5],
                      title="Función de Transferencia",
                      xlabel="Frecuencia (Hz)",
                      ylabel="Amplificación",
                      label="No disponible")
            end
        catch e
            println("    Error en función de transferencia: ", e)
        end

        # Posición 5 (2,2): Perfil de deformaciones
        println("  Generando perfil de deformaciones...")
        try
            if haskey(resultados_movimientos, "gamma_max") && length(subcapas) > 0
                gamma_max_datos = resultados_movimientos["gamma_max"]
                profundidades = Float64[]
                deformaciones = Float64[]
                
                z_acum = 0.0
                for (i, gamma) in enumerate(gamma_max_datos)
                    # Punto al inicio de la subcapa
                    push!(profundidades, z_acum)
                    push!(deformaciones, gamma)  # gamma_max ya está en porcentaje
                    
                    # Actualizar profundidad
                    if i <= length(subcapas)
                        z_acum += subcapas[i].h
                    else
                        z_acum += 10.0
                    end
                    
                    # Punto al final de la subcapa
                    push!(profundidades, z_acum)
                    push!(deformaciones, gamma)  # gamma_max ya está en porcentaje
                end
                
                plot!(p[5], deformaciones, profundidades,
                      title="Perfil de Deformaciones",
                      xlabel="Deformación Cortante γ (%)",
                      ylabel="Profundidad (m)",
                      yflip=true,
                      lw=3,
                      color=:purple,
                      legend=false,
                      grid=true)
                
                if length(deformaciones) > 0
                    def_max = maximum(abs.(deformaciones))
                    # Asegurar que el rango X sea apropiado para deformaciones < 1%
                    x_max = max(def_max*1.2, 0.1)  # Mínimo 0.1% para visualización
                    xlims!(p[5], (0, x_max))
                    ylims!(p[5], (0, maximum(profundidades)*1.05))
                end
            else
                plot!(p[5], [0, 0.1], [0, 500],
                      title="Perfil de Deformaciones",
                      xlabel="Deformación Cortante γ (%)",
                      ylabel="Profundidad (m)",
                      yflip=true,
                      label="No hay datos disponibles")
            end
        catch e
            println("    Error en deformaciones: ", e)
        end

        # Posición 6 (3,1): Velocidades de Corte
        println("  Generando gráfico de velocidades de corte...")
        try
            if capas_finales !== nothing && length(subcapas) > 0
                profundidades = Float64[]
                velocidades = Float64[]
                
                z_acum = 0.0
                for (i, subcapa) in enumerate(subcapas)
                    vs_final = subcapa.Vs_inicial
                    if i <= length(capas_finales)
                        vs_final = capas_finales[i].Vs_actual
                    end
                    
                    push!(profundidades, z_acum)
                    push!(velocidades, vs_final)
                    
                    z_acum += subcapa.h
                    
                    push!(profundidades, z_acum)
                    push!(velocidades, vs_final)
                end
                
                if length(velocidades) > 2
                    plot!(p[6], velocidades, profundidades, 
                          title="Perfil de Velocidades de Onda de Corte",
                          xlabel="Velocidad Vs (m/s)",
                          ylabel="Profundidad (m)",
                          yflip=true,
                          lw=3, 
                          color=:blue,
                          legend=false,
                          grid=true)
                    
                    v_min = 0.0
                    v_max = maximum(velocidades) * 1.05
                    xlims!(p[6], (v_min, v_max))
                    ylims!(p[6], (0, maximum(profundidades)*1.05))
                else
                    plot!(p[6], [0, 1000], [0, 500], 
                          title="Perfil de Velocidades de Onda de Corte",
                          xlabel="Velocidad Vs (m/s)",
                          ylabel="Profundidad (m)",
                          yflip=true,
                          label="No disponible")
                end
            else
                plot!(p[6], [0, 1000], [0, 500], 
                      title="Perfil de Velocidades de Onda de Corte",
                      xlabel="Velocidad Vs (m/s)",
                      ylabel="Profundidad (m)",
                      yflip=true,
                      label="No disponible")
            end
        catch e
            println("    Error en velocidades: ", e)
        end

        # Posición 7 (3,2): Amortiguamiento
        println("  Generando gráfico de amortiguamiento...")
        try
            if capas_finales !== nothing && length(subcapas) > 0
                profundidades = Float64[]
                amortiguamientos = Float64[]
                
                z_acum = 0.0
                for (i, subcapa) in enumerate(subcapas)
                    zeta_final = subcapa.ζ_inicial
                    if i <= length(capas_finales)
                        zeta_final = capas_finales[i].ζ_actual
                    end
                    
                    zeta_pct = zeta_final * 100.0
                    
                    push!(profundidades, z_acum)
                    push!(amortiguamientos, zeta_pct)
                    
                    z_acum += subcapa.h
                    
                    push!(profundidades, z_acum)
                    push!(amortiguamientos, zeta_pct)
                end
                
                if length(amortiguamientos) > 2
                    plot!(p[7], amortiguamientos, profundidades, 
                          title="Perfil de Amortiguamiento",
                          xlabel="Amortiguamiento ζ (%)",
                          ylabel="Profundidad (m)",
                          yflip=true,
                          lw=4,
                          color=:red,
                          legend=false,
                          grid=true)
                    
                    a_min = 0.0
                    a_max = max(maximum(amortiguamientos) * 1.2, 5.0)
                    xlims!(p[7], (a_min, a_max))
                    ylims!(p[7], (0, maximum(profundidades)*1.05))
                else
                    plot!(p[7], [0, 5], [0, 500], 
                          title="Perfil de Amortiguamiento",
                          xlabel="Amortiguamiento ζ (%)",
                          ylabel="Profundidad (m)",
                          yflip=true,
                          label="No disponible")
                end
            else
                plot!(p[7], [0, 5], [0, 500], 
                      title="Perfil de Amortiguamiento",
                      xlabel="Amortiguamiento ζ (%)",
                      ylabel="Profundidad (m)",
                      yflip=true,
                      label="No disponible")
            end
        catch e
            println("    Error en amortiguamiento: ", e)
        end
        
        # Guardar el gráfico
        println("  Guardando gráfico combinado...")
        savefig(p, archivo_salida)
        println("   Gráfico guardado: $archivo_salida")
        println()
        
        return p
        
    catch e
        println("Error al generar gráficos: $e")
        return nothing
    end
end

# Funciones auxiliares
function generar_grafico_velocidades_corte(subcapas, capas_finales=nothing, capas_originales=nothing)
    println("Función auxiliar de velocidades disponible")
    return nothing
end

function generar_grafico_amortiguamiento(subcapas, capas_finales=nothing, capas_originales=nothing)
    println("Función auxiliar de amortiguamiento disponible")
    return nothing
end

end # module
