# Ejemplo Básico - Análisis Sísmico de 3 Capas
# Este script demuestra el uso básico del programa con un modelo simple

using FFTW, Plots, DelimitedFiles, Statistics

# Incluir todas las librerías
include("../Lib_TransferenciaCore.jl")
include("../Lib_Sismico.jl")
include("../Lib_Modelo.jl")
include("../Lib_Preparacion.jl")
include("../Lib_AnalisisLEQ.jl")
include("../Lib_Movimientos.jl")
include("../Lib_Espectros.jl")
include("../Lib_Amplificacion.jl")
include("../Lib_Graficos.jl")
include("../Lib_Resumen.jl")

using .Lib_TransferenciaCore, .Lib_Sismico, .Lib_Modelo, .Lib_Preparacion
using .Lib_AnalisisLEQ, .Lib_Movimientos, .Lib_Espectros, .Lib_Amplificacion
using .Lib_Graficos, .Lib_Resumen

# Configurar gráficos
gr()

# ===== CONFIGURACIÓN BÁSICA =====
println("EJEMPLO BÁSICO - ANÁLISIS SÍSMICO DE 3 CAPAS")
println("=" ^ 50)

# Configuración simplificada para demostración
config_basico = Dict(
    "archivo_sismo" => "../sismo_recortado.txt",  # Usar sismo de ejemplo
    "factor_escala" => 1.0,
    "espesor_maximo_subcapas" => 2.0,  # Discretización más gruesa para rapidez
    "max_iter" => 50,                   # Menos iteraciones para ejemplo
    "tolerancia" => 0.05,               # Tolerancia más relajada (5%)
    "max_iter_conv" => 100,
    "tolerancia_conv" => 1e-3,
    "archivo_salida" => "ejemplo_basico.png",
    "archivo_deconvolucion" => "ejemplo_deconvolucion.txt"
)

# ===== MODELO DE SUELO BÁSICO =====
# Modelo simple de 3 capas típico

# Parámetros para suelo blando (no lineal)
params_arcilla = parametros_degradacion(1.2, 0.1, -0.1, 0.25, 0.01, "equivalente")

# Parámetros para suelo duro (lineal)
params_arena = parametros_degradacion(10000, 10000, 0.0, 0.0, 0.03, "lineal")
params_roca = parametros_degradacion(10000, 10000, 0.0, 0.0, 0.01, "lineal")

# Definir capas básicas
capas_basico = [
    # Arcilla blanda superficial (20m)
    capa_suelo(1600.0, 200.0, 0.05, 20.0, params_arcilla, "Arcilla blanda"),
    
    # Arena densa intermedia (30m)
    capa_suelo(1900.0, 600.0, 0.03, 30.0, params_arena, "Arena densa"),
    
    # Roca basal
    capa_suelo(2400.0, 1500.0, 0.01, 0.0, params_roca, "Roca")
]

# ===== EJECUTAR ANÁLISIS =====

try
    # Procesar sismo
    tiempo, aceleraciones, dt = leer_sismo(config_basico["archivo_sismo"])
    dt_esc, acc_esc, T_pred = escalar_sismo(dt, aceleraciones, config_basico["factor_escala"])
    
    # Modelo de suelo
    subcapas = subdividir_capas(capas_basico, config_basico["espesor_maximo_subcapas"])
    T_fund = periodo_fundamental_sitio(subcapas)
    
    # Preparar análisis
    freq, omega, f_acc = preparar_analisis(tiempo, acc_esc, dt_esc)
    
    # Análisis LEQ
    capas_fin, conv, hist, Amp_i, Amp_f = ejecutar_analisis_leq(
        subcapas, omega, f_acc, 
        config_basico["max_iter"], config_basico["tolerancia"]
    )
    
    # Movimientos y deformaciones
    resultados = calcular_movimientos_y_deformaciones(
        capas_fin, subcapas, omega, f_acc, Amp_f, tiempo, dt, 
        config_basico["max_iter_conv"], config_basico["tolerancia_conv"],
        config_basico["archivo_deconvolucion"]
    )
    
    # Espectros
    espectros = calcular_espectros_y_amplificacion(tiempo, resultados, freq)
    
    # Amplificación
    fun_amp = funcion_amplificacion(resultados["acc_superficie"], resultados["acc_rock"], dt)
    
    # Gráficos
    generar_graficos_completos(
        tiempo, resultados, espectros, subcapas, omega, f_acc,
        resultados["gamma_max"], config_basico["archivo_salida"],
        capas_fin, capas_basico
    )
    
    # Resumen
    generar_resumen_final(
        capas_basico, capas_fin, conv, hist, T_fund, 50.0, T_pred,
        resultados, espectros, resultados["gamma_max"], 
        config_basico["archivo_salida"], config_basico["espesor_maximo_subcapas"]
    )
    
    println("\n✅ EJEMPLO BÁSICO COMPLETADO EXITOSAMENTE")
    println("📁 Revisa los archivos generados:")
    println("   - $(config_basico["archivo_salida"]) (gráficos)")
    println("   - $(config_basico["archivo_deconvolucion"]) (señal deconvolucionada)")
    
catch e
    println("\n❌ ERROR en el ejemplo básico:")
    println(e)
    println("\n💡 Verifica que:")
    println("   - El archivo de sismo existe")
    println("   - Las librerías están en el directorio correcto")
    println("   - Las dependencias de Julia están instaladas")
end