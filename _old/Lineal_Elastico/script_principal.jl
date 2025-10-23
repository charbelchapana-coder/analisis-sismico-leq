using Revise
using FFTW, Plots, DelimitedFiles
includet("Lib_LElastico.jl")
import .Lib_LElastico: capa_suelo, funcion_transferencia, espectro_chopra, deformaciones_corte_shake

# Configurar backend para gráficos
gr()

## Ejemplo de cálculo

# Ingreso de sismo
cd("C:/Users/charb/OneDrive/Desktop/Modelos Plaxis/Deconvolucion/Julia")
input = readdlm("sismo_recortado.txt", ',')

tiempo = input[:,1]
aceleracion = input[:,2]
dt = tiempo[2] - tiempo[1]
N = length(aceleracion)
freq = collect(0.01:1/(N*dt):1/(2*dt))
omega = 2π .* freq

capas = [capa_suelo(1400, 500, 0.05, 220),
         capa_suelo(2000, 850, 0.03, 280),
         capa_suelo(2500, 2050, 0.02, 0)] # Último vector es la capa de roca, con espesor cero

# Calcular función de transferencia
Amp = funcion_transferencia(capas, omega)
H = Amp.^(-1)

# Aplicar a la FFT del sismo en superficie
f_acc = fft(aceleracion)
f_acc_pos = f_acc[1:length(freq)]
f_rock = Amp .* f_acc_pos
f_recons = H .* f_rock

# Reconstituyendo la señal de forma simétrica para la IFFT
f_rock_full = [f_rock; conj(reverse(f_rock[2:end-1]))]
f_recons_full = [f_recons; conj(reverse(f_recons[2:end-1]))]   

# Obteniendo las señales en el dominio del tiempo
acc_rock = real(ifft(f_rock_full))
acc_recons = real(ifft(f_recons_full))

# Asegurar que los vectores tengan la misma longitud que el tiempo
min_len = min(length(acc_rock), length(tiempo))
acc_rock = acc_rock[1:min_len]
tiempo_adj = tiempo[1:min_len]

min_len2 = min(length(acc_recons), length(tiempo_adj))
acc_recons = acc_recons[1:min_len2]
tiempo_adj = tiempo_adj[1:min_len2]
aceleracion_adj = aceleracion[1:min_len2]

# Espectros de respuesta
periodo = freq.^(-1)
xi = 0.05
Sa_superficie = espectro_chopra(aceleracion_adj, dt, periodo, xi)
Sa_roca = espectro_chopra(acc_rock, dt, periodo, xi)
Sa_recons = espectro_chopra(acc_recons, dt, periodo, xi)

# Calcular PGA (Peak Ground Acceleration)
PGA_superficie = maximum(abs.(aceleracion_adj))
PGA_roca = maximum(abs.(acc_rock))

# Calcular deformaciones de corte máximas usando método SHAKE
gamma_max, gamma_superficie, gamma_roca, gamma_interfaces = deformaciones_corte_shake(capas, omega, f_rock)

# Mostrar resultados de deformaciones de corte
println("\n=== DEFORMACIONES DE CORTE MÁXIMAS ===")
println("Roca basal: γ_max = $(round(gamma_roca, digits=4))%")
for i in 1:(length(capas)-1)
    println("Capa $i (Vs=$(capas[i].Vs) m/s, h=$(capas[i].h) m): γ_max = $(round(gamma_max[i], digits=4))%")
    if i < length(capas)-1  # Si no es la última capa, mostrar interface
        println("  Interface $i-$(i+1): γ_max = $(round(gamma_interfaces[i], digits=4))%")
    end
end
println("Superficie: γ_max = $(round(gamma_superficie, digits=4))%")

# Gráficos señales en el tiempo en superficie, en roca y reconstruida
p1 = plot(tiempo_adj,aceleracion_adj,
         xlabel="Tiempo (s)",ylabel="Aceleración (m/s²)",
         title="Señales en el tiempo", 
         label = "Señal en superficie", 
         legend=:topright,
         titlefontsize=10, guidefontsize=9, tickfontsize=8, legendfontsize=8)  
p2 = plot!(tiempo_adj,acc_recons, label = "Señal reconstruida")
p3 = plot!(tiempo_adj,acc_rock, label = "Señal en roca")

# Gráficos de espectros de aceleración con etiquetas de PGA
p4 = plot(periodo,Sa_superficie,
         xscale=:log10,xlabel="Periodo (s)",ylabel="Aceleración (m/s²)", 
         title = "Espectros de aceleración", 
         label = "Señal en superficie", 
         legend=:topright,
         titlefontsize=10, guidefontsize=9, tickfontsize=8, legendfontsize=8)
p5 = plot!(periodo,Sa_recons, label = "Señal reconstruida")
p6 = plot!(periodo,Sa_roca, label = "Señal en roca")

# Agregar puntos y etiquetas de PGA
# Para T=0 (PGA), usamos un período muy pequeño para la visualización
T_pga = 0.01  # Período muy pequeño para representar PGA
plot!(p4, [T_pga], [PGA_superficie], seriestype=:scatter, markersize=6, markercolor=:blue, label="")
plot!(p4, [T_pga], [PGA_roca], seriestype=:scatter, markersize=6, markercolor=:green, label="")

# Posicionar automáticamente las etiquetas de PGA
# Calcular posiciones para evitar solapamiento
y_max = max(maximum(Sa_superficie), maximum(Sa_roca), maximum(Sa_recons))
y_offset_superficie = PGA_superficie * 1.3
y_offset_roca = PGA_roca * 0.7

annotate!(p4, T_pga, y_offset_superficie, text("PGA superficie = $(round(PGA_superficie, digits=3)) m/s²", :blue, :left, 8))
annotate!(p4, T_pga, y_offset_roca, text("PGA roca = $(round(PGA_roca, digits=3)) m/s²", :green, :left, 8))

# Gráfico de Función de Transferencia
p7 = plot(periodo,abs.(Amp.^(-1)),
         xscale=:log10,xlabel="Periodo (s)",ylabel="|H(ω)|", 
         title = "Función de Transferencia",
         legend=false,
         titlefontsize=10, guidefontsize=9, tickfontsize=8)

# Gráfico de perfil de deformaciones de corte
# Calcular profundidades en el centro de cada capa
n_capas = length(capas) - 1
profundidades = zeros(n_capas)

# Calcular profundidades usando una función auxiliar
function calcular_profundidades(capas_suelo, n)
    prof = zeros(n)
    acumulado = 0.0
    for i in 1:n
        acumulado += capas_suelo[i].h / 2  # Centro de la capa
        prof[i] = acumulado
        acumulado += capas_suelo[i].h / 2  # Final de la capa
    end
    return prof
end

profundidades = calcular_profundidades(capas, n_capas)

# Agregar profundidades para superficie, interfaces y roca basal
prof_superficie = 0.0
prof_roca_basal = sum(capa.h for capa in capas[1:n_capas])

# Calcular profundidades de las interfaces entre capas
function calcular_interfaces(capas_suelo, n)
    prof_int = zeros(n - 1)
    acumulado = 0.0
    for i in 1:(n - 1)
        acumulado += capas_suelo[i].h
        prof_int[i] = acumulado
    end
    return prof_int
end

prof_interfaces = calcular_interfaces(capas, n_capas)

# Combinar TODOS los datos para el gráfico (superficie, capas, interfaces, roca basal)
# Ordenar por profundidad para una visualización coherente
gamma_completo = [gamma_superficie; gamma_max; gamma_interfaces; gamma_roca]
prof_completa = [prof_superficie; profundidades; prof_interfaces; prof_roca_basal]

# Ordenar por profundidad
indices_orden = sortperm(prof_completa)
gamma_completo = gamma_completo[indices_orden]
prof_completa = prof_completa[indices_orden]

# Crear perfil de deformaciones de corte con marcadores diferenciados
p8 = plot(xlabel="Deformación de corte máxima (%)", 
         ylabel="Profundidad (m)", 
         title="Perfil de Deformaciones de Corte",
         legend=false,
         grid=true,
         minorgrid=true,
         yflip=true,  # Invertir eje Y para mostrar profundidad hacia abajo
         titlefontsize=10, guidefontsize=9, tickfontsize=8)

# Agregar superficie
scatter!(p8, [gamma_superficie], [prof_superficie], 
         markersize=10, markercolor=:blue, markershape=:star8,
         markerstrokewidth=2, markerstrokecolor=:darkblue, 
         label="Superficie")

# Agregar capas
scatter!(p8, gamma_max, profundidades, 
         markersize=8, markercolor=:red, markershape=:circle,
         markerstrokewidth=2, markerstrokecolor=:darkred,
         label="Capas")

# Agregar interfaces
scatter!(p8, gamma_interfaces, prof_interfaces, 
         markersize=6, markercolor=:orange, markershape=:diamond,
         markerstrokewidth=2, markerstrokecolor=:darkorange,
         label="Interfaces")

# Agregar roca basal
scatter!(p8, [gamma_roca], [prof_roca_basal], 
         markersize=10, markercolor=:green, markershape=:square,
         markerstrokewidth=2, markerstrokecolor=:darkgreen,
         label="Roca basal")

# Conectar puntos con líneas
plot!(p8, gamma_completo, prof_completa, 
      linewidth=2, 
      linecolor=:red, 
      linestyle=:dash)

# Agregar líneas horizontales para mostrar las separaciones entre capas
x_max = maximum([maximum(gamma_completo) * 1.3, 0.02])  # Asegurar que las líneas sean visibles

# Función para calcular profundidades de separación entre capas
function calcular_separaciones(capas_suelo, n)
    prof_sep = [0.0]  # Superficie
    acumulado = 0.0
    for i in 1:n
        acumulado += capas_suelo[i].h
        push!(prof_sep, acumulado)
    end
    return prof_sep
end

prof_separaciones = calcular_separaciones(capas, n_capas)

# Dibujar líneas de separación entre capas
for i in 1:length(prof_separaciones)
    plot!(p8, [0, x_max], [prof_separaciones[i], prof_separaciones[i]], 
          linewidth=1, linecolor=:black, linestyle=:solid, alpha=0.5)
end

# Agregar etiquetas con valores de deformación
# Superficie
if gamma_superficie > 0.001
    annotate!(p8, gamma_superficie * 1.1, prof_superficie, 
             text("$(round(gamma_superficie, digits=4))%", :red, :left, 8))
end

# Capas
for i in 1:length(gamma_max)
    if gamma_max[i] > 0.001  # Solo mostrar etiquetas para valores significativos
        annotate!(p8, gamma_max[i] * 1.1, profundidades[i], 
                 text("$(round(gamma_max[i], digits=4))%", :red, :left, 8))
    end
end

# Roca basal
if gamma_roca > 0.001
    annotate!(p8, gamma_roca * 1.1, prof_roca_basal, 
             text("$(round(gamma_roca, digits=4))%", :red, :left, 8))
end

# Agregar etiquetas de identificación de capas en el lado derecho
# Superficie
annotate!(p8, x_max * 0.8, prof_superficie, 
         text("Superficie", :blue, :center, 7))

# Capas
for i in 1:n_capas
    # Posición en el centro de cada capa
    prof_centro = profundidades[i]
    
    # Texto con información de la capa
    texto_capa = "Capa $i\nVs=$(Int(capas[i].Vs)) m/s\nh=$(Int(capas[i].h)) m"
    
    # Colocar etiqueta en el lado derecho del gráfico
    annotate!(p8, x_max * 0.8, prof_centro, 
             text(texto_capa, :blue, :center, 7))
end

# Roca basal
annotate!(p8, x_max * 0.8, prof_roca_basal, 
         text("Roca basal\nVs=$(Int(capas[end].Vs)) m/s", :blue, :center, 7))

# Agregar etiquetas de profundidad de medición
# Superficie
annotate!(p8, x_max * 0.05, prof_superficie, 
         text("$(Int(prof_superficie))m", :black, :left, 6))

# Capas
for i in 1:length(profundidades)
    annotate!(p8, x_max * 0.05, profundidades[i], 
             text("$(Int(profundidades[i]))m", :black, :left, 6))
end

# Roca basal
annotate!(p8, x_max * 0.05, prof_roca_basal, 
         text("$(Int(prof_roca_basal))m", :black, :left, 6))

# Mantener el eje Y numérico sin reemplazar con nombres
# (removemos la línea de yticks personalizada)

# Gráfico final con 4 subplots - Configuración mejorada para visualización
plot(p1, p4, p7, p8, 
     layout=(2,2), 
     size=(1600,1200),  # Tamaño mayor para mejor legibilidad
     margin=15Plots.mm,  # Márgenes generosos
     left_margin=20Plots.mm,  # Margen izquierdo extra para etiquetas Y
     bottom_margin=15Plots.mm,  # Margen inferior para etiquetas X
     top_margin=15Plots.mm,  # Margen superior para títulos
     titlefontsize=12,  # Tamaño de fuente para títulos
     guidefontsize=10,  # Tamaño de fuente para etiquetas de ejes
     tickfontsize=8,  # Tamaño de fuente para números en ejes
     legendfontsize=9,  # Tamaño de fuente para leyenda
     plot_title="Análisis de Respuesta Sísmica y Deformaciones de Corte",
     plot_titlefontsize=16)
savefig("analisis_respuesta_sismica.png")  # Guardar figura final
