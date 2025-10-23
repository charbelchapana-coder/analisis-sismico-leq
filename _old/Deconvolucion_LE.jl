using LinearAlgebra, FFTW, Plots, DelimitedFiles
cd("C:/Users/charb/OneDrive/Desktop/Modelos Plaxis/Deconvolucion/Julia")

## Definiciones previas de estructuras y funciones

# Definición de estructuras para capas de suelo y roca
struct capa_suelo
    ρ::Float64
    Vs::Float64
    ζ::Float64
    h::Float64
end

# Definición de función de degradacion de rigidez (G/Go)
function curva_degradacion(γ, a, γref)
    if γ > 1
        Ge = 1/(1+(1/γref)^a)
    else
        Ge = 1/(1+(γ/γref)^a)
    end
    return Ge
end

# Definición de función de amortiguamiento (D)
function curva_amortiguamiento(γ, a, b, c)
    if γ > 1
        De = a + b + c
    else
        De = a*γ^2 + b*γ + c
    end
    return De
end

# Función para calcular la matriz de transferencia de una capa de suelo a otra
function matriz_T(suelo_1::capa_suelo, suelo_2::capa_suelo, ω::Float64)
    # Parámetros de la primera capa
    ρ1, Vs1, ζ1, h1 = suelo_1.ρ, suelo_1.Vs, suelo_1.ζ, suelo_1.h
    # Parámetros de la segunda capa
    ρ2, Vs2, ζ2, h2 = suelo_2.ρ, suelo_2.Vs, suelo_2.ζ, suelo_2.h

    # Módulo de corte dinámico complejo
    G1 = ρ1 * Vs1^2 * (1 + 2im * ζ1)
    G2 = ρ2 * Vs2^2 * (1 + 2im * ζ2)

    # Número de onda complejo
    k1 = ρ1 * ω^2 / G1
    k2 = ρ2 * ω^2 / G2

    # Impedancia compleja
    αm = (k1*G1)/(k2*G2)

    # Funciones auxiliares
    epos = 0.5*exp(im*k1*h1)
    eneg = 0.5*exp(-im*k1*h1)

    # Matriz de transferencia 
    M = [(1+αm)*epos (1-αm)*eneg;
          (1-αm)*epos (1+αm)*eneg]
    return M
end

function funcion_transferencia(suelos::Vector{capa_suelo}, ω::Vector{Float64})
    H = zeros(ComplexF64, 2, length(ω)) # Inicializar vector para todos los resultados
    Amp = zeros(ComplexF64, length(ω))

    for i in 1:length(ω)
        # Inicialización de la matriz de transferencia total
        M_total = Matrix{ComplexF64}(I, 2, 2) # Corregido: matriz identidad compleja
        N = length(suelos)

        # Parámetros de la primera capa
        ρ1, Vs1, ζ1, h1 = suelos[1].ρ, suelos[1].Vs, suelos[1].ζ, suelos[1].h
        # Parámetros de la última capa
        ρn, Vsn, ζn, hn = suelos[end].ρ, suelos[end].Vs, suelos[end].ζ, suelos[end].h

        # Módulo de corte dinámico complejo
        G1 = ρ1 * Vs1^2 * (1 + 2im * ζ1)
        Gn = ρn * Vsn^2 * (1 + 2im * ζn)

        # Número de onda complejo
        k1 = ρ1 * ω[i]^2 / G1
        k2 = ρn * ω[i]^2 / Gn

        # Impedancia compleja
        αm = (k1*G1)/(k2*Gn)

        # Funciones auxiliares
        epos = 0.5*exp(im*k1*h1)
        eneg = 0.5*exp(-im*k1*h1)

        # Cálculo de la matriz de transferencia total a través de las capas de suelo
        for j in 1:(length(suelos)-1)
            M = matriz_T(suelos[j], suelos[j+1], ω[i])
            M_total = M * M_total
        end
        H[1,i] = M_total[1,1] + M_total[1,2]
        H[2,i] = M_total[2,1] + M_total[2,2]
        denom = 2*H[1,i]
        if abs(denom) < 1e-12
            Amp[i] = 0.0
        else
            Amp[i] = (H[1,i] + H[2,i]) / denom
        end
    end
    return Amp
end

function espectro_chopra(aceleracion, dt, periodos, xi)
    Sa = zeros(length(periodos))
    m = 1.0
    for k in 1:length(periodos)
        T = periodos[k]
        if T == 0
            Sa[k] = maximum(abs.(aceleracion))
            continue
        end
        w_n = 2 * π / T
        K = w_n^2 * m
        w_d = w_n * sqrt(1 - xi^2)
        e_xi_wn_dt = exp(-xi * w_n * dt)
        sin_wd_dt = sin(w_d * dt)
        cos_wd_dt = cos(w_d * dt)

        A = e_xi_wn_dt * (xi / sqrt(1 - xi^2) * sin_wd_dt + cos_wd_dt)
        B = e_xi_wn_dt * (1 / w_d * sin_wd_dt)
        C = (1 / K) * (2 * xi / (w_n * dt) + e_xi_wn_dt * (((1 - 2 * xi^2) / (w_d * dt) - xi / sqrt(1 - xi^2)) * sin_wd_dt - (1 + 2 * xi / (w_n * dt)) * cos_wd_dt))
        D = (1 / K) * (1 - 2 * xi / (w_n * dt) - e_xi_wn_dt * (((1 - 2 * xi^2) / (w_d * dt)) * sin_wd_dt - (2 * xi / (w_n * dt)) * cos_wd_dt))

        Aprime = -w_n / sqrt(1 - xi^2) * e_xi_wn_dt * sin_wd_dt
        Bprime = e_xi_wn_dt * (cos_wd_dt - xi / sqrt(1 - xi^2) * sin_wd_dt)
        Cprime = (1 / K) * (-1 / dt + e_xi_wn_dt * ((w_n / sqrt(1 - xi^2) + xi / (dt * sqrt(1 - xi^2))) * sin_wd_dt + (1 / dt) * cos_wd_dt))
        Dprime = (1 / (K * dt)) * (1 - e_xi_wn_dt * ((xi / sqrt(1 - xi^2)) * sin_wd_dt + cos_wd_dt))

        u = zeros(length(aceleracion))
        u_dot = zeros(length(aceleracion))
        p = -m .* aceleracion

        for i in 1:(length(aceleracion) - 1)
            u[i+1] = A * u[i] + B * u_dot[i] + C * p[i] + D * p[i+1]
            u_dot[i+1] = Aprime * u[i] + Bprime * u_dot[i] + Cprime * p[i] + Dprime * p[i+1]
        end

        acel_total = -(2 * xi * w_n .* u_dot + w_n^2 .* u)
        Sa[k] = maximum(abs.(acel_total))
    end
    return Sa
end

## Ejemplo de cálculo

input = readdlm("sismo_recortado.txt",',')

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
acc_rock = [0; acc_rock]
acc_recons = [0; acc_recons]

# Espectros de respuesta
periodo = freq.^(-1)
xi = 0.05
Sa_superficie = espectro_chopra(aceleracion, dt, periodo, xi)
Sa_roca = espectro_chopra(acc_rock, dt, periodo, xi)
Sa_recons = espectro_chopra(acc_recons, dt, periodo, xi)

# Gráficos señales en el tiempo en superficie, en roca y reconstruida
p1 = plot(tiempo,aceleracion,xlabel="Tiempo (s)",ylabel="Aceleración (m/s²)",title="Señales en el tiempo")  
p2 = plot!(tiempo,acc_recons)
p3 = plot!(tiempo,acc_rock)

# Gráficos de espectros de aceleración
p4 = plot(periodo,Sa_superficie,xscale=:log10,xlabel="Periodo (s)",ylabel="Aceleración (m/s²)", title = "Espectros de aceleración")
p5 = plot!(periodo,Sa_recons)
p6 = plot!(periodo,Sa_roca)

# Gráfico de Función de Transferencia
p7 = plot(freq,abs.(Amp.^(-1)),xscale=:log10,xlabel="Frecuencia (Hz)",ylabel="|H(ω)|", title = "Función de Transferencia",legend=false)

# Gráfico final
plot(p1,p4, p7,layout=(3,1),size=(800,900),font="Arial",label=["Señal en superficie" "Señal reconstruida" "Señal en roca"],minorgrid=true)
savefig("resultados.png")