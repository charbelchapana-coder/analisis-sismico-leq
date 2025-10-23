# Script Principal - Análisis Sísmico Completo con Método Lineal Equivalente
# Configuración modular optimizada - utiliza 9 librerías especializadas

using FFTW, Plots, DelimitedFiles, Statistics

# Incluir todas las librerías modulares optimizadas
include("Lib_TransferenciaCore.jl")  # Funciones auxiliares compartidas
include("Lib_Sismico.jl")           # Etapa 1: Datos sísmicos
include("Lib_Modelo.jl")            # Etapa 2: Modelo de suelo
include("Lib_Preparacion.jl")       # Etapa 3: Preparación
include("Lib_AnalisisLEQ.jl")       # Etapa 4: Análisis LEQ
include("Lib_Movimientos.jl")       # Etapa 5: Movimientos
include("Lib_Espectros.jl")         # Etapa 6: Espectros
include("Lib_Amplificacion.jl")     # Etapa 7: Amplificación
include("Lib_Graficos.jl")          # Etapa 8: Gráficos
include("Lib_Resumen.jl")           # Etapa 9: Resumen final

using .Lib_TransferenciaCore, .Lib_Sismico, .Lib_Modelo, .Lib_Preparacion
using .Lib_AnalisisLEQ, .Lib_Movimientos, .Lib_Espectros, .Lib_Amplificacion
using .Lib_Graficos, .Lib_Resumen

# Configurar backend para gráficos
gr()

# ===== CONFIGURACIÓN DEL ANÁLISIS =====

# Configuración del sismo
archivo_sismo = "sismo_recortado.txt"
factor_escala = 1.0

# Configuración de subdivisión de capas
espesor_maximo_subcapas = 10.0  # metros

# Configuración del análisis lineal equivalente
max_iter = 200
tolerancia = 0.1  # 0.1% para análisis más estable

# Configuración de la convolución iterativa
max_iter_conv = 200  # Máximo de iteraciones para convolución (configurable por usuario)
tolerancia_conv = 1e-4  # Tolerancia de convergencia para convolución (configurable por usuario)

# Archivo de salida
archivo_salida = "analisis_completo_lineal_equivalente_subcapas.png"

# ===== DEFINICIÓN DEL MODELO DE SUELO =====

# Parámetros para comportamiento lineal equivalente (suelo diatomáceo)
params_diatomaceo = parametros_degradacion(
    1.4, 0.14, -0.15, 0.27, 0.00773, "equivalente"
)

# Parámetros para comportamiento lineal (brecha)
params_lineal_brecha = parametros_degradacion(
    10000, 10000, 0.0, 0.0, 0.03, "lineal"
)

# Parámetros para comportamiento lineal (roca)
params_lineal_roca = parametros_degradacion(
    10000, 10000, 0.0, 0.0, 0.02, "lineal"
)

# Definir capas de suelo
capas = [
    # Capa 1: Suelo diatomaceo (comportamiento no lineal)
    capa_suelo(1400.0, 500.0, 0.05, 220.0, params_diatomaceo, "Suelo diatomáceo"),
    
    # Capa 2: Brecha sedimentaria (comportamiento lineal)
    capa_suelo(2000.0, 850.0, 0.03, 280.0, params_lineal_brecha, "Brecha sedimentaria"),
    
    # Capa 3: Roca basal (comportamiento lineal)
    capa_suelo(2500.0, 2050.0, 0.02, 0.0, params_lineal_roca, "Roca basal")
]

# ===== CONFIGURACIÓN COMPLETA =====
config = Dict(
    "archivo_sismo" => archivo_sismo,
    "factor_escala" => factor_escala,
    "capas" => capas,
    "espesor_maximo_subcapas" => espesor_maximo_subcapas,
    "max_iter" => max_iter,
    "tolerancia" => tolerancia,
    "max_iter_conv" => max_iter_conv,
    "tolerancia_conv" => tolerancia_conv,
    "archivo_salida" => archivo_salida
)

# ===== EJECUCIÓN DEL ANÁLISIS MODULAR =====

println("INICIANDO ANÁLISIS SÍSMICO CON ARQUITECTURA MODULAR OPTIMIZADA")
println("============================================================")

# ETAPA 1: Procesar datos sísmicos
tiempo, aceleraciones, dt = leer_sismo(archivo_sismo)
dt_escalado, aceleraciones_escaladas, periodo_predominante = escalar_sismo(dt, aceleraciones, factor_escala)

# DIAGNÓSTICO DE SEÑAL DE ENTRADA
println("\nDIAGNÓSTICO DE SEÑAL SÍSMICA DE ENTRADA:")
println("============================================================")
PGA_entrada = maximum(abs.(aceleraciones_escaladas))
println("  PGA señal de entrada: $(round(PGA_entrada, digits=4)) m/s²")
println("  Número de puntos: $(length(aceleraciones_escaladas))")
println("  Duración total: $(round(tiempo[end], digits=2)) s")
println("  Periodo predominante: $(round(periodo_predominante, digits=3)) s")

if PGA_entrada > 2.0
    println("  ⚠️  ADVERTENCIA: PGA de entrada parece muy alto (>2.0 m/s²)")
elseif PGA_entrada < 0.5
    println("  ⚠️  ADVERTENCIA: PGA de entrada parece muy bajo (<0.5 m/s²)")
else
    println("  ✅ PGA de entrada en rango esperado (0.5-2.0 m/s²)")
end

# ETAPA 2: Procesar modelo de suelo  
subcapas = subdividir_capas(capas, espesor_maximo_subcapas)
periodo_fundamental = periodo_fundamental_sitio(subcapas)

# ETAPA 3: Preparar análisis
frecuencias, omega, f_acc_pos = preparar_analisis(tiempo, aceleraciones_escaladas, dt_escalado)

# ETAPA 4: Análisis lineal equivalente
capas_finales, convergencia, historial, Amp_inicial, Amp_final = ejecutar_analisis_leq(
    subcapas, omega, f_acc_pos, max_iter, tolerancia
)

# ETAPA 5: Calcular movimientos y deformaciones
resultados_movimientos = calcular_movimientos_y_deformaciones(
    capas_finales, subcapas, omega, f_acc_pos, Amp_final, tiempo, dt, max_iter_conv, tolerancia_conv
)

# ETAPA 6: Calcular espectros de respuesta
resultados_espectros = calcular_espectros_y_amplificacion(
    tiempo, resultados_movimientos, frecuencias
)

# ETAPA 7: Calcular función de amplificación
fun_amplificacion = funcion_amplificacion(
    resultados_movimientos["acc_superficie"], 
    resultados_movimientos["acc_rock"], 
    dt
)

# ETAPA 8: Generar gráficos
generar_graficos_completos(
    tiempo, resultados_movimientos, resultados_espectros, 
    subcapas, omega, f_acc_pos, 
    resultados_movimientos["gamma_max"], archivo_salida,
    capas_finales  # Agregar capas finales para cálculo de derivada
)

# ETAPA 9: Generar resumen final
generar_resumen_final(
    capas, capas_finales, convergencia, historial, 
    periodo_fundamental, sum([c.h for c in capas if c.h > 0]),
    periodo_predominante, resultados_movimientos, resultados_espectros,
    resultados_movimientos["gamma_max"], archivo_salida
)