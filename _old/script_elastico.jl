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

# Calcular deformaciones de corte máximas en el centro de cada capa
gamma_max = deformaciones_corte_shake(capas, omega, f_rock)

# Mostrar resultados de deformaciones de corte
println("\n=== DEFORMACIONES DE CORTE MÁXIMAS EN EL CENTRO DE CADA CAPA ===")
for i in 1:(length(capas)-1)
    println("Capa $i (Vs=$(capas[i].Vs) m/s, h=$(capas[i].h) m): γ_max = $(round(gamma_max[i], digits=4))%")
end

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

# Gráfico de perfil de deformaciones de corte tipo escalones
n_capas = length(capas) - 1

# Calcular profundidades de las interfaces entre capas (tope y base de cada capa)
prof_interfaces = zeros(n_capas + 1)  # Incluye superficie (0) y base de la última capa
prof_interfaces[1] = 0.0  # Superficie

for i in 1:n_capas
    prof_interfaces[i+1] = prof_interfaces[i] + capas[i].h
end

# Crear vectores para el gráfico tipo escalones
gamma_escalones = Float64[]
prof_escalones = Float64[]

# Para cada capa, crear el escalón correspondiente
for i in 1:n_capas
    prof_top = prof_interfaces[i]      # Tope de la capa
    prof_bottom = prof_interfaces[i+1] # Base de la capa
    gamma_capa = gamma_max[i]          # Deformación de la capa
    
    # Agregar puntos para formar el escalón
    push!(gamma_escalones, gamma_capa)   # Valor constante en el tope
    push!(prof_escalones, prof_top)
    
    push!(gamma_escalones, gamma_capa)   # Valor constante en la base
    push!(prof_escalones, prof_bottom)
end

# Filtrar valores positivos para escala logarítmica
gamma_escalones_pos = [max(val, 1e-6) for val in gamma_escalones]  # Evitar valores <= 0

# Definir los ticks del eje X en escala logarítmica (10^-5 a 10^0)
x_ticks = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0]
x_labels = ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰"]

# Crear el gráfico
p8 = plot(gamma_escalones_pos, prof_escalones,
         xlabel="Deformación de corte máxima (%)", 
         ylabel="Profundidad (m)", 
         title="Perfil de Deformaciones de Corte (Escalones)",
         legend=false,
         grid=true,
         minorgrid=false,  # Desactivar minorgrid para evitar problemas con escala log
         yflip=true,  # Invertir eje Y para mostrar profundidad hacia abajo
         xscale=:log10,  # Escala logarítmica base 10 para el eje X
         xticks=(x_ticks, x_labels),  # Ticks personalizados para el eje X
         xlims=(1e-5, 1e0),  # Límites del eje X
         linewidth=3,
         linecolor=:red,
         titlefontsize=10, guidefontsize=9, tickfontsize=8)

# Agregar líneas horizontales desde la curva hasta el eje vertical
x_max = maximum(gamma_escalones_pos) * 1.2
x_min_curva = minimum(gamma_escalones_pos)  # Límite izquierdo donde está la curva

for i in 1:length(prof_interfaces)
    # Dibujar líneas desde el punto más a la izquierda de la curva hasta el eje vertical (x_min)
    plot!(p8, [x_min_curva, 1e-5], [prof_interfaces[i], prof_interfaces[i]], 
          linewidth=1, linecolor=:gray, linestyle=:dash, alpha=0.7)
end

# Agregar etiquetas con valores de deformación en el centro de cada capa
for i in 1:n_capas
    prof_centro = (prof_interfaces[i] + prof_interfaces[i+1]) / 2
    gamma_capa = gamma_max[i]
    
    if gamma_capa > 0.001  # Solo mostrar etiquetas para valores significativos
        # Colocar texto rojo en la parte inferior del centro de la capa
        prof_texto_rojo = prof_centro + (prof_interfaces[i+1] - prof_interfaces[i]) * 0.2
        annotate!(p8, gamma_capa + x_max * 0.05, prof_texto_rojo, 
                 text("$(round(gamma_capa, digits=4))%", :red, :left, 9))
    end
end

# Agregar etiquetas de identificación de capas a la derecha de la curva
for i in 1:n_capas
    prof_centro = (prof_interfaces[i] + prof_interfaces[i+1]) / 2
    gamma_capa = gamma_max[i]
    
    # Texto con información de la capa
    texto_capa = "Capa $i\nVs=$(Int(capas[i].Vs)) m/s\nh=$(Int(capas[i].h)) m"
    
    # Colocar texto azul a la derecha de la curva roja con justificación izquierda
    prof_texto_azul = prof_centro - (prof_interfaces[i+1] - prof_interfaces[i]) * 0.2
    x_pos_texto = gamma_capa * 2.5  # Posición a la derecha del valor de deformación
    annotate!(p8, x_pos_texto, prof_texto_azul, 
             text(texto_capa, :blue, :left, 8))
end

# Agregar etiquetas de profundidad encima de las líneas grises
for i in 1:length(prof_interfaces)
    # Posición ligeramente a la derecha del eje vertical y encima de la línea
    x_pos_prof = 1.2e-5  # Ligeramente a la derecha del eje vertical
    y_pos_prof = prof_interfaces[i] - 5  # Encima de la línea (5m hacia arriba)
    annotate!(p8, x_pos_prof, y_pos_prof, 
             text("$(Int(prof_interfaces[i]))m", :black, :left, 7))
end

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
savefig("resultado_analisis_lineal.png")