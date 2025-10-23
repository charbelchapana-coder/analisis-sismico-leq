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

# ===== CONFIGURACIÓN PRINCIPAL DEL ANÁLISIS =====
# Modifica estos parámetros según tus necesidades

# 1. CONFIGURACIÓN DEL SISMO
archivo_sismo = "sismo_recortado.txt"         # Archivo con el registro sísmico
factor_escala = 1.0                           # Factor de escala para el sismo (1.0 = sin escala)

# 2. CONFIGURACIÓN DE DISCRETIZACIÓN
espesor_maximo_subcapas = 1.0                 # Espesor máximo de subcapas en metros
                                              # Valores recomendados: 1.0-10.0 m
                                              # Menor valor = mayor precisión, mayor tiempo de cálculo

# 3. CONFIGURACIÓN DEL ANÁLISIS LINEAL EQUIVALENTE
max_iter = 200                                # Máximo número de iteraciones
tolerancia = 0.01                             # Tolerancia de convergencia en % (0.01% recomendado)

# 4. CONFIGURACIÓN DE LA CONVOLUCIÓN ITERATIVA
max_iter_conv = 200                           # Máximo de iteraciones para convolución
tolerancia_conv = 1e-4                        # Tolerancia de convergencia para convolución

# 5. CONFIGURACIÓN DE SALIDA
archivo_salida = "analisis_LEQ.png"           # Archivo de salida para gráficos

# ===== DEFINICIÓN DEL MODELO DE SUELO =====
# Modifica estos parámetros según tu proyecto específico

# PARÁMETROS DE COMPORTAMIENTO NO LINEAL (para suelos blandos)
# Usa estos parámetros para suelos que muestran degradación con la deformación
# Orden de parámetros: (a, γ_ref, b1, b2, b3, tipo)
params_diatomaceo = parametros_degradacion(
    1.4,        # a: parámetro para curva de degradación G/Go = 1/(1+(γ/γ_ref)^a)
    0.14,       # γ_ref: deformación de referencia (sin unidades, sin transformar a %)
    -0.15,      # b1: coeficiente cuadrático para curva de amortiguamiento ζ = b1*γ² + b2*γ + b3
    0.27,       # b2: coeficiente lineal para curva de amortiguamiento
    0.00773,    # b3: término independiente para curva de amortiguamiento
    "equivalente"  # tipo: comportamiento no lineal
)

# PARÁMETROS DE COMPORTAMIENTO LINEAL (para suelos duros/rocas)
# Usa estos parámetros para materiales que no muestran degradación significativa
params_lineal_brecha = parametros_degradacion(
    10000,      # a: valor alto para evitar degradación
    10000,      # γ_ref: valor alto para evitar degradación
    0.0,        # b1: sin componente cuadrática
    0.0,        # b2: sin componente lineal
    0.02,       # b3: amortiguamiento constante (2%)
    "lineal"    # tipo: comportamiento lineal
)

params_lineal_roca = parametros_degradacion(
    10000,      # a: valor alto para evitar degradación
    10000,      # γ_ref: valor alto para evitar degradación
    0.0,        # b1: sin componente cuadrática
    0.0,        # b2: sin componente lineal
    0.02,       # b3: amortiguamiento constante (2%)
    "lineal"    # tipo: comportamiento lineal
)

# DEFINICIÓN DE CAPAS DE SUELO
# Formato: capa_suelo(densidad, Vs_inicial, amortiguamiento, espesor, parámetros, nombre)
# - densidad: kg/m³
# - Vs_inicial: velocidad de corte inicial en m/s
# - amortiguamiento: amortiguamiento inicial (fracción, ej. 0.05 = 5%)
# - espesor: espesor de la capa en metros (0.0 para roca basal)
# - parámetros: usar params_diatomaceo para suelos blandos, params_lineal_* para duros
# - nombre: identificador de la capa (para resultados)

# NOTA: Para agregar más capas, simplemente añade más líneas siguiendo el formato anterior
# NOTA: La última capa debe tener espesor 0.0 (roca basal)
# NOTA: Para cambiar nombres, modifica el último parámetro de cada capa_suelo

capas = [
    # Capa 1: Suelo blando superior (comportamiento no lineal)
    capa_suelo(1400.0, 500.0, 0.05, 280.0, params_diatomaceo, "Suelo diatomáceo"),
    
    # Capa 2: Suelo duro intermedio (comportamiento lineal)
    capa_suelo(2000.0, 850.0, 0.02, 222.0, params_lineal_brecha, "Brecha sedimentaria"),
    
    # Capa 3: Roca basal (comportamiento lineal, sin límite de profundidad)
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
    capas_finales,  # Agregar capas finales para cálculo de derivada
    capas          # Agregar capas originales para interfaces
)

# ETAPA 9: Generar resumen final
generar_resumen_final(
    capas, capas_finales, convergencia, historial, 
    periodo_fundamental, sum([c.h for c in capas if c.h > 0]),
    periodo_predominante, resultados_movimientos, resultados_espectros,
    resultados_movimientos["gamma_max"], archivo_salida, espesor_maximo_subcapas
)