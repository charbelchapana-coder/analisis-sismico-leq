module Lib_Graficos

using Plots
using FFTW  
using Statistics
using ..Lib_TransferenciaCore: capa_suelo, funcion_transferencia

export generar_graficos_completos

"""
Función auxiliar para calcular deformaciones de corte usando derivada de desplazamientos (γ = du/dz)
"""
function calcular_deformaciones_por_derivada_desplazamientos(capas_finales, omega, f_rock_final)
    try
        n_capas = length(capas_finales) - 1  # Excluir roca basal
        gamma_derivada = zeros(n_capas)
        prof_gamma_derivada = zeros(n_capas)
        
        # Evitar división por cero
        omega_safe = [max(w, 1e-6) for w in omega]
        
        # Calcular desplazamientos en cada interface
        despl_interfaces = []
        prof_interfaces = [0.0]  # Superficie
        
        # Desplazamiento en superficie
        H_superficie = funcion_transferencia(capas_finales, omega)
        despl_superficie_fft = H_superficie .* (-f_rock_final ./ (omega_safe.^2))
        despl_superficie_fft[1] = 0.0
        despl_superficie_full = [despl_superficie_fft; conj(reverse(despl_superficie_fft[2:end-1]))]
        despl_superficie_t = real(FFTW.ifft(despl_superficie_full))
        despl_max_superficie = maximum(abs.(despl_superficie_t))
        push!(despl_interfaces, despl_max_superficie)
        
        # Desplazamientos en interfaces entre capas
        z_acum = 0.0
        for i in 1:n_capas
            z_acum += capas_finales[i].h
            push!(prof_interfaces, z_acum)
            
            # Función de transferencia hasta esta profundidad
            H_interface = funcion_transferencia(capas_finales[i+1:end], omega)
            despl_interface_fft = H_interface .* (-f_rock_final ./ (omega_safe.^2))
            despl_interface_fft[1] = 0.0
            despl_interface_full = [despl_interface_fft; conj(reverse(despl_interface_fft[2:end-1]))]
            despl_interface_t = real(FFTW.ifft(despl_interface_full))
            despl_max_interface = maximum(abs.(despl_interface_t))
            push!(despl_interfaces, despl_max_interface)
        end
        
        # Calcular deformaciones por diferenciación numérica: γ = du/dz
        for i in 1:n_capas
            if i < length(despl_interfaces)
                du = despl_interfaces[i] - despl_interfaces[i+1]  # Diferencia de desplazamientos
                dz = prof_interfaces[i+1] - prof_interfaces[i]     # Espesor de la capa
                
                if dz > 0
                    gamma_derivada[i] = abs(du / dz) * 100.0  # Convertir a porcentaje
                end
                
                # Profundidad al centro de la capa
                prof_gamma_derivada[i] = (prof_interfaces[i] + prof_interfaces[i+1]) / 2
            end
        end
        
        return gamma_derivada, prof_gamma_derivada
        
    catch e
        println("Error en cálculo de deformaciones por derivada: $e")
        return zeros(length(capas_finales)-1), zeros(length(capas_finales)-1)
    end
end

function generar_graficos_completos(tiempo::Vector{Float64}, resultados_movimientos::Dict, 
                                   resultados_espectros::Dict, subcapas::Vector{capa_suelo}, 
                                   omega::Vector{Float64}, f_acc_pos::Vector{ComplexF64}, 
                                   gamma_max::Vector{Float64}, archivo_salida::String, 
                                   capas_finales::Union{Vector{capa_suelo}, Nothing}=nothing)
    
    println("\\n8. GENERACION DE GRAFICOS")
    println("------------------------------------------------------------")
    
    try
        gr()
        p = plot(layout=(2,2), size=(1200, 900), dpi=150)
        
        # 1. Acelerogramas
        println("  Generando acelerogramas...")
        try
            if haskey(resultados_movimientos, "acc_rock") && haskey(resultados_movimientos, "acc_superficie")
                acc_rock = resultados_movimientos["acc_rock"]
                acc_superficie = resultados_movimientos["acc_superficie"]
                
                n_total = min(length(tiempo), length(acc_rock), length(acc_superficie))
                if n_total > 10
                    step = max(1, div(n_total, 2000))
                    idx = 1:step:n_total
                    
                    plot!(p[1], tiempo[idx], acc_superficie[idx], 
                          title="Acelerogramas", label="Superficie", lw=1.5, color=:red)
                    plot!(p[1], tiempo[idx], acc_rock[idx], 
                          label="Roca", lw=1.5, color=:blue)
                    
                    if haskey(resultados_movimientos, "acc_superficie_convolucion")
                        acc_conv = resultados_movimientos["acc_superficie_convolucion"]
                        if length(acc_conv) >= n_total
                            plot!(p[1], tiempo[idx], acc_conv[idx], 
                                  label="Reconstruida", lw=1.5, color=:green, linestyle=:dash)
                        end
                    end
                    
                    if haskey(resultados_movimientos, "acc_superficie_alternativa")
                        acc_alt = resultados_movimientos["acc_superficie_alternativa"]
                        if length(acc_alt) >= n_total
                            plot!(p[1], tiempo[idx], acc_alt[idx], 
                                  label="Alternativa", lw=1.5, color=:cyan, linestyle=:dot)
                        end
                    end
                    
                    xlabel!(p[1], "Tiempo (s)")
                    ylabel!(p[1], "Aceleracion (m/s2)")
                    xlims!(p[1], (0, tiempo[end]))
                end
            else
                plot!(p[1], [0], [0], title="Acelerogramas", label="No disponible")
            end
        catch e
            println("    Error en acelerogramas: ", e)
            plot!(p[1], [0], [0], title="Acelerogramas", label="Error")
        end
        
        # 2. Funcion de transferencia
        println("  Generando funcion de transferencia...")
        try
            if haskey(resultados_movimientos, "H_transferencia") && length(omega) > 0
                H_transfer = resultados_movimientos["H_transferencia"]
                if length(H_transfer) > 0
                    freq_hz = omega ./ (2π)
                    H_direct = abs.(H_transfer)
                    
                    n_max = min(length(freq_hz), length(H_direct))
                    mask = (freq_hz[1:n_max] .>= 0.1) .& (freq_hz[1:n_max] .<= 50.0) .& (.!isinf.(H_direct[1:n_max]))
                    
                    if sum(mask) > 5
                        plot!(p[2], freq_hz[1:n_max][mask], H_direct[1:n_max][mask], 
                              title="Funcion de Transferencia", xscale=:log10, 
                              legend=false, lw=2, color=:blue)
                        xlabel!(p[2], "Frecuencia (Hz)")
                        ylabel!(p[2], "|H(w)|")
                        
                        plot!(p[2], [0.1, 50], [1, 1], color=:gray, linestyle=:dot, 
                              legend=false, alpha=0.7)
                    else
                        plot!(p[2], [1, 10], [1, 1], title="Funcion de Transferencia", 
                              label="Datos insuficientes", xscale=:log10)
                    end
                else
                    plot!(p[2], [1, 10], [1, 1], title="Funcion de Transferencia", 
                          label="No disponible", xscale=:log10)
                end
            else
                plot!(p[2], [1, 10], [1, 1], title="Funcion de Transferencia", 
                      label="No disponible", xscale=:log10)
            end
        catch e
            println("    Error en funcion de transferencia: ", e)
            plot!(p[2], [1, 10], [1, 1], title="Funcion de Transferencia", 
                  label="Error", xscale=:log10)
        end
        
        # 3. Espectros de respuesta
        println("  Generando espectros de respuesta...")
        try
            if haskey(resultados_espectros, "periodos")
                periodos = resultados_espectros["periodos"]
                
                if haskey(resultados_espectros, "Sa_superficie")
                    sa_sup = resultados_espectros["Sa_superficie"]
                    if length(periodos) > 0 && length(sa_sup) > 0
                        n_max = min(length(periodos), length(sa_sup))
                        mask = (sa_sup[1:n_max] .> 0) .& (.!isnan.(sa_sup[1:n_max])) .& (periodos[1:n_max] .> 0)
                        
                        if sum(mask) > 5
                            plot!(p[3], periodos[1:n_max][mask], sa_sup[1:n_max][mask], 
                                  title="Espectros de Respuesta", xlabel="Periodo (s)", ylabel="Sa (m/s2)",
                                  xscale=:log10, label="Superficie", lw=2, color=:red)
                        end
                    end
                end
                
                if haskey(resultados_espectros, "Sa_rock")
                    sa_rock = resultados_espectros["Sa_rock"]
                    if length(periodos) > 0 && length(sa_rock) > 0
                        n_max = min(length(periodos), length(sa_rock))
                        mask = (sa_rock[1:n_max] .> 0) .& (.!isnan.(sa_rock[1:n_max])) .& (periodos[1:n_max] .> 0)
                        
                        if sum(mask) > 5
                            plot!(p[3], periodos[1:n_max][mask], sa_rock[1:n_max][mask], 
                                  label="Roca (Deconvolución)", lw=2, color=:blue, linestyle=:dash)
                            println("    ✅ Espectro de roca deconvolucionada incluido en el gráfico")
                        end
                    end
                end
                
                if haskey(resultados_espectros, "Sa_convolucion")
                    sa_conv = resultados_espectros["Sa_convolucion"]
                    if length(periodos) > 0 && length(sa_conv) > 0
                        n_max = min(length(periodos), length(sa_conv))
                        mask = (sa_conv[1:n_max] .> 0) .& (.!isnan.(sa_conv[1:n_max])) .& (periodos[1:n_max] .> 0)
                        
                        if sum(mask) > 5
                            plot!(p[3], periodos[1:n_max][mask], sa_conv[1:n_max][mask], 
                                  label="Reconstruida", lw=2, color=:green, linestyle=:dot)
                        end
                    end
                end
                
                if haskey(resultados_espectros, "Sa_alternativa")
                    sa_alt = resultados_espectros["Sa_alternativa"]
                    if length(periodos) > 0 && length(sa_alt) > 0
                        n_max = min(length(periodos), length(sa_alt))
                        mask = (sa_alt[1:n_max] .> 0) .& (.!isnan.(sa_alt[1:n_max])) .& (periodos[1:n_max] .> 0)
                        
                        if sum(mask) > 5
                            plot!(p[3], periodos[1:n_max][mask], sa_alt[1:n_max][mask], 
                                  label="Alternativa", lw=2, color=:cyan, linestyle=:dashdot)
                        end
                    end
                end
            else
                plot!(p[3], [0.1, 10], [0.1, 0.1], title="Espectros de Respuesta", 
                      label="No disponible", xscale=:log10)
            end
        catch e
            println("    Error en espectros: ", e)
            plot!(p[3], [0.1, 10], [0.1, 0.1], title="Espectros de Respuesta", 
                  label="Error", xscale=:log10)
        end
        
        # 4. Perfil de deformaciones
        println("  Generando perfil de deformaciones...")
        try
            if length(gamma_max) > 0 && length(subcapas) > 0
                n_capas = min(length(gamma_max), length(subcapas))
                if n_capas > 1
                    # Crear perfil suavizado desde superficie hacia profundidad
                    profundidades_plot = [0.0]  # Superficie
                    gamma_plot = [0.0]          # Deformación cero en superficie
                    
                    # Calcular profundidades en el centro de cada capa
                    z_acum = 0.0
                    for i in 1:n_capas
                        if !isnan(gamma_max[i]) && gamma_max[i] > 0
                            z_centro = z_acum + subcapas[i].h / 2  # Centro de la capa
                            push!(profundidades_plot, z_centro)
                            push!(gamma_plot, gamma_max[i])  # gamma_max ya está en %
                        end
                        z_acum += subcapas[i].h
                    end
                    
                    # Agregar punto en la base con deformación reducida gradualmente
                    if length(profundidades_plot) > 1
                        push!(profundidades_plot, z_acum)  # Profundidad total
                        push!(gamma_plot, gamma_plot[end] * 0.5)  # Reducir en la base
                    end
                    
                    if length(profundidades_plot) > 2
                        # Asegurar valores mínimos para escala logarítmica
                        gamma_plot_safe = [max(g, 0.001) for g in gamma_plot]
                        
                        # Crear gráfico con curva verde (SHAKE/STRATA) sin leyenda
                        plot!(p[4], gamma_plot_safe, profundidades_plot, 
                              title="Perfil de Deformaciones", yflip=true, 
                              xlabel="gamma (%)", ylabel="Profundidad (m)", 
                              lw=2, color=:green, marker=:circle, markersize=3,
                              xscale=:log10, xlims=(0.001, 1.0),
                              xticks=([0.001, 0.01, 0.1, 1.0], ["0.001", "0.01", "0.1", "1.0"]),
                              legend=false)
                    else
                        plot!(p[4], [0.001, 0.1], [0, 100], title="Perfil de Deformaciones", yflip=true, 
                              xlabel="gamma (%)", ylabel="Profundidad (m)",
                              xscale=:log10, xlims=(0.001, 1.0),
                              xticks=([0.001, 0.01, 0.1, 1.0], ["0.001", "0.01", "0.1", "1.0"]),
                              legend=false)
                    end
                else
                    plot!(p[4], [0.001, 0.1], [0, 100], title="Perfil de Deformaciones", yflip=true, 
                          xlabel="gamma (%)", ylabel="Profundidad (m)",
                          xscale=:log10, xlims=(0.001, 1.0),
                          xticks=([0.001, 0.01, 0.1, 1.0], ["0.001", "0.01", "0.1", "1.0"]),
                          legend=false)
                end
            else
                plot!(p[4], [0.001, 0.1], [0, 100], title="Perfil de Deformaciones", yflip=true, 
                      xlabel="gamma (%)", ylabel="Profundidad (m)",
                      xscale=:log10, xlims=(0.001, 1.0),
                      xticks=([0.001, 0.01, 0.1, 1.0], ["0.001", "0.01", "0.1", "1.0"]),
                      legend=false)
            end
        catch e
            println("    Error en perfil: ", e)
            plot!(p[4], [0.001, 0.1], [0, 100], title="Perfil de Deformaciones", yflip=true, 
                  xlabel="gamma (%)", ylabel="Profundidad (m)",
                  xscale=:log10, xlims=(0.001, 1.0),
                  xticks=([0.001, 0.01, 0.1, 1.0], ["0.001", "0.01", "0.1", "1.0"]),
                  legend=false)
        end
        
        println("  Guardando grafico como: ", archivo_salida)
        savefig(p, archivo_salida)
        println("  ✅ Grafico guardado exitosamente")
        
        return p
        
    catch e
        println("ERROR en generacion de graficos: ", e)
        return nothing
    end
end

end # module