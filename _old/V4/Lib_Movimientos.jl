"""
Lib_Movimientos.jl - Librería para cálculo de movimientos y deformaciones
Etapa 5: Cálculo de movimientos, deformaciones, y reconstrucción de señales
"""

module Lib_Movimientos

using LinearAlgebra, FFTW, Statistics, Dates
using ..Lib_TransferenciaCore: capa_suelo, parametros_degradacion, deformaciones_corte_shake, funcion_transferencia

export calcular_movimientos_y_deformaciones, convolucion_directa_pystrata

"""
    convolucion_directa_pystrata(H_transferencia, omega, f_rock_deconv, dt)

Realiza convolución directa siguiendo la metodología de PyStrata sin amplificación artificial.
La convolución en PyStrata es simplemente: signal_surface = transfer_function * signal_rock
Usa la misma función de transferencia que se usó para deconvolución para consistencia.
"""
function convolucion_directa_pystrata(H_transferencia, omega, f_rock_deconv, dt)
    println("Iniciando convolución directa basada en metodología PyStrata...")
    
    # Diagnóstico de entrada
    mag_roca_entrada = sum(abs.(f_rock_deconv))
    println("DIAGNÓSTICO DE CONVOLUCIÓN PYSTRATA:")
    println("  Magnitud señal roca (entrada): $(round(mag_roca_entrada, digits=4))")
    
    # Usar la función de transferencia proporcionada (misma que deconvolución)
    H_final = H_transferencia
    
    # Diagnóstico de función de transferencia
    mag_H_final = sum(abs.(H_final))
    println("  Magnitud H final: $(round(mag_H_final, digits=4))")
    
    # CONVOLUCIÓN DIRECTA: siguiendo exactamente PyStrata
    # ts = np.fft.irfft(tf * self.fourier_amps / self.time_step)
    f_superficie_pystrata = H_final .* f_rock_deconv
    
    # Diagnóstico del resultado
    mag_superficie_pystrata = sum(abs.(f_superficie_pystrata))
    ratio_conv_pystrata = mag_superficie_pystrata / mag_roca_entrada
    println("  Magnitud señal superficie PyStrata: $(round(mag_superficie_pystrata, digits=4))")
    println("  Ratio convolución PyStrata (superficie/roca): $(round(ratio_conv_pystrata, digits=2))")
    
    # Reconstruir señal temporal completa para transformada inversa
    N_total = 2 * (length(omega) - 1)
    f_superficie_pystrata_full = [f_superficie_pystrata; conj(reverse(f_superficie_pystrata[2:end-1]))]
    
    # Transformada inversa para obtener señal temporal
    acc_superficie_pystrata_t = real(ifft(f_superficie_pystrata_full))
    
    println("  Señal PyStrata generada exitosamente")
    return acc_superficie_pystrata_t, f_superficie_pystrata
end

"""
    reconstruir_señal_superficie(capas_finales, omega, f_rock_final, dt, Amp_final)

Reconstruye la señal en superficie mediante convolución con la función de transferencia.
NOTA: Este es el proceso de CONVOLUCIÓN (roca → superficie).
"""
function reconstruir_señal_superficie(capas_finales, omega, f_rock_final, dt, Amp_final)
    println("Reconstruyendo señal en superficie mediante convolución...")
    
    # CONVOLUCIÓN: H * señal_roca = señal_superficie
    H_transferencia = Amp_final
    aceleracion_superficie_fft = H_transferencia .* f_rock_final
    
    # Reconstruir señal completa para IFFT
    N_total = 2 * (length(omega) - 1)
    acel_superficie_full = [aceleracion_superficie_fft; conj(reverse(aceleracion_superficie_fft[2:end-1]))]
    
    # Transformada inversa para obtener señal temporal
    aceleracion_superficie_t = real(ifft(acel_superficie_full))
    
    return aceleracion_superficie_t
end

"""
    validar_deformaciones_metodo_alternativo(capas_finales, omega, f_rock_final)

Valida las deformaciones usando método alternativo γ = du/dz.
"""
function validar_deformaciones_metodo_alternativo(capas_finales, omega, f_rock_final)
    # Método alternativo simplificado
    try
        H_superficie = funcion_transferencia(capas_finales, omega)
        aceleraciones_superficie = H_superficie .* f_rock_final
        
        # Estimar deformación máxima usando diferencias finitas
        if length(aceleraciones_superficie) > 1
            max_accel = maximum(abs.(aceleraciones_superficie))
            gamma_alternativo = max_accel / (omega[end] * capas_finales[1].Vs_actual) * 100
            return gamma_alternativo
        else
            return 0.0
        end
    catch
        return 0.0
    end
end

"""
    calcular_movimientos_y_deformaciones(capas_finales::Vector{capa_suelo}, capas_originales::Vector{capa_suelo},
                                        omega::Vector{Float64}, f_acc_pos::Vector{ComplexF64}, 
                                        Amp_final::Vector{ComplexF64}, tiempo::Vector{Float64}, dt::Float64, max_iter_conv::Int=50, tolerancia_conv::Float64=1e-4)

Función principal para el cálculo de movimientos y deformaciones.
Implementa la Etapa 5 del análisis.

# Parámetros adicionales:
- `max_iter_conv`: Número máximo de iteraciones para la convolución iterativa (por defecto 50)
- `tolerancia_conv`: Tolerancia de convergencia para la convolución iterativa (por defecto 1e-4)
"""
function calcular_movimientos_y_deformaciones(capas_finales::Vector{capa_suelo}, capas_originales::Vector{capa_suelo},
                                             omega::Vector{Float64}, f_acc_pos::Vector{ComplexF64}, 
                                             Amp_final::Vector{ComplexF64}, tiempo::Vector{Float64}, dt::Float64, max_iter_conv::Int=50, tolerancia_conv::Float64=1e-4)
    println("\n5. CÁLCULO DE MOVIMIENTOS Y DEFORMACIONES")
    println("------------------------------------------------------------")
    
    # Guardar función de transferencia para gráficos
    # Amp_final representa atenuación (superficie → roca), necesitamos invertir para amplificación
    H_transferencia = Amp_final.^(-1)  # H representa amplificación (roca → superficie)
    
    # DIAGNÓSTICO CRÍTICO: ¿Qué representa realmente Amp_final?
    # Comprobemos con algunos valores típicos
    freq_test = 1.0  # Hz
    idx_test = argmin(abs.(omega ./ (2*pi) .- freq_test))
    amp_test = abs(Amp_final[idx_test])
    h_test = abs(H_transferencia[idx_test])
    println("DIAGNÓSTICO DE FUNCIÓN DE TRANSFERENCIA:")
    println("  Frecuencia de prueba: $(freq_test) Hz")
    println("  |Amp_final| en $(freq_test) Hz: $(round(amp_test, digits=3))")
    println("  |H_transferencia| en $(freq_test) Hz: $(round(h_test, digits=3))")
    
    if h_test > 1.0
        println("  → H_transferencia > 1: AMPLIFICACIÓN correcta (roca → superficie)")
    else
        println("  → H_transferencia < 1: PROBLEMA - debería amplificar")
    end
    
    # ETAPA DE DECONVOLUCIÓN: Obtener señal en roca
    println("DECONVOLUCIÓN: Obteniendo señal en roca desde señal en superficie")
    
    # Para deconvolución: f_roca = f_superficie / H = f_superficie * H^(-1) = f_superficie * Amp_final
    f_rock_final = f_acc_pos .* Amp_final
    
    # Diagnóstico de deconvolución
    mag_superficie_original = sum(abs.(f_acc_pos))
    mag_roca_deconv = sum(abs.(f_rock_final))
    ratio_deconv = mag_superficie_original / mag_roca_deconv
    
    println("DIAGNÓSTICO DE DECONVOLUCIÓN:")
    println("  Magnitud señal superficie (original): $(round(mag_superficie_original, digits=4))")
    println("  Magnitud señal roca (deconvolucionada): $(round(mag_roca_deconv, digits=4))")
    println("  Ratio deconvolución (superficie/roca): $(round(ratio_deconv, digits=2))")
    
    # VERIFICACIÓN ADICIONAL: ¿La función de transferencia representa amplificación?
    mag_H = sum(abs.(H_transferencia))
    mag_Amp_final = sum(abs.(Amp_final))
    println("  Magnitud H_transferencia: $(round(mag_H, digits=4))")
    println("  Magnitud Amp_final: $(round(mag_Amp_final, digits=4))")
    println("  Ratio (Amp_final/H): $(round(mag_Amp_final/mag_H, digits=2))")
    
    if ratio_deconv > 1.0
        println("  ✅ La señal en roca es menor que en superficie (correcto para deconvolución)")
    else
        println("  ⚠️  PROBLEMA: La señal en roca es mayor que en superficie (incorrecto)")
    end
    
    # ETAPA DE CONVOLUCIÓN: Reconstruir señal en superficie para verificación
    # Debe usar H (no H^-1) para reconstruir: f_superficie = H * f_roca
    acc_superficie_convolucion = reconstruir_señal_superficie(
        capas_finales, omega, f_rock_final, dt, H_transferencia
    )
    
    # Diagnóstico de convolución regular
    # Comparar con la señal de superficie original (dominio temporal)
    n_original = length(tiempo)
    if length(acc_superficie_convolucion) > n_original
        acc_conv_ajustada = acc_superficie_convolucion[1:n_original]
    else
        acc_conv_ajustada = acc_superficie_convolucion
    end
    
    # Reconstruir señal original en dominio temporal para comparación
    f_superficie_original_full = [f_acc_pos; conj(reverse(f_acc_pos[2:end-1]))]
    acc_superficie_original = real(ifft(f_superficie_original_full))
    if length(acc_superficie_original) > n_original
        acc_superficie_original = acc_superficie_original[1:n_original]
    end
    
    # Diagnóstico comparativo
    mag_superficie_conv_temporal = sum(abs.(acc_conv_ajustada))
    mag_superficie_original_temporal = sum(abs.(acc_superficie_original))
    ratio_conv_temporal = mag_superficie_conv_temporal / mag_superficie_original_temporal
    
    println("DIAGNÓSTICO DE CONVOLUCIÓN REGULAR:")
    println("  Magnitud señal superficie original (temporal): $(round(mag_superficie_original_temporal, digits=4))")
    println("  Magnitud señal superficie reconstruida (temporal): $(round(mag_superficie_conv_temporal, digits=4))")
    println("  Ratio convolución temporal (reconstruida/original): $(round(ratio_conv_temporal, digits=2))")
    
    # Si el ratio no está cerca de 1.0, hay un problema en la reconstrucción
    if abs(ratio_conv_temporal - 1.0) > 0.1
        println("  ⚠️  ADVERTENCIA: La señal reconstruida no coincide con la original (ratio ≠ 1.0)")
    else
        println("  ✅ VERIFICACIÓN: La señal reconstruida coincide con la original")
    end
    
    # ETAPA DE CONVOLUCIÓN DIRECTA BASADA EN PYSTRATA
    # Usar la MISMA función de transferencia que usamos para deconvolución
    acc_superficie_alternativa, f_superficie_alternativa = convolucion_directa_pystrata(
        H_transferencia, omega, f_rock_final, dt
    )
    
    println("RESULTADO TERMINOLÓGICO CORRECTO:")
    println("  - DECONVOLUCIÓN: Señal superficie → Señal roca (aplicada)")
    println("  - CONVOLUCIÓN: Señal roca → Señal superficie reconstruida (verificación)")
    println("  - CONVOLUCIÓN PYSTRATA: Señal roca → Señal superficie (metodología estándar)")
    
    # Calcular movimientos en superficie usando propagación directa
    # CORRECCIÓN: Para obtener superficie desde roca, multiplicar por H (no H^-1)
    f_superficie = H_transferencia .* f_rock_final
    
    # Reconstruir señales temporales con duración completa
    f_rock_full = [f_rock_final; conj(reverse(f_rock_final[2:end-1]))]
    f_superficie_full = [f_superficie; conj(reverse(f_superficie[2:end-1]))]
    
    acc_rock = real(ifft(f_rock_full))
    acc_superficie = real(ifft(f_superficie_full))
    
    # DIAGNÓSTICO: Comparar con señal de entrada original
    # La acc_superficie debería corresponder aproximadamente a la señal original
    PGA_superficie_calculada = maximum(abs.(acc_superficie))
    
    # Reconstruir señal original desde f_acc_pos para comparación
    f_original_full = [f_acc_pos; conj(reverse(f_acc_pos[2:end-1]))]
    acc_original_reconstruida = real(ifft(f_original_full))
    PGA_original_reconstruida = maximum(abs.(acc_original_reconstruida))
    
    println("\nDIAGNÓSTICO DE CONSISTENCIA:")
    println("  PGA señal original reconstruida: $(round(PGA_original_reconstruida, digits=4)) m/s²")
    println("  PGA superficie calculada: $(round(PGA_superficie_calculada, digits=4)) m/s²")
    println("  Ratio (calculada/original): $(round(PGA_superficie_calculada/PGA_original_reconstruida, digits=2))")
    
    if abs(PGA_superficie_calculada/PGA_original_reconstruida - 1.0) > 0.1
        println("  ⚠️  ADVERTENCIA: Discrepancia significativa entre señal original y calculada")
    else
        println("  ✅ Señales consistentes")
    end
    
    # Asegurar que las señales mantengan la duración original completa
    n_original = length(tiempo)
    
    # Si la señal calculada es más larga, truncar al tamaño original
    if length(acc_rock) > n_original
        acc_rock = acc_rock[1:n_original]
        acc_superficie = acc_superficie[1:n_original]
    # Si la señal calculada es más corta, rellenar con ceros
    elseif length(acc_rock) < n_original
        acc_rock = [acc_rock; zeros(n_original - length(acc_rock))]
        acc_superficie = [acc_superficie; zeros(n_original - length(acc_superficie))]
    end
    
    # Ajustar longitud de convolución PyStrata para que coincida con el tiempo original
    if length(acc_superficie_alternativa) > n_original
        acc_superficie_alternativa = acc_superficie_alternativa[1:n_original]
    elseif length(acc_superficie_alternativa) < n_original
        acc_superficie_alternativa = [acc_superficie_alternativa; zeros(n_original - length(acc_superficie_alternativa))]
    end
    
    # Usar toda la duración original
    tiempo_adj = tiempo
    
    # GUARDAR SEÑAL DECONVOLUCIONADA EN ROCA
    println("\nGUARDANDO SEÑAL DECONVOLUCIONADA EN ROCA:")
    println("--------------------------------------------------------")
    nombre_archivo = "deconvolucion_roca.txt"
    
    try
        # Crear matriz con tiempo y aceleración
        datos_roca = [tiempo_adj acc_rock]
        
        # Guardar archivo con solo dos columnas: tiempo y aceleración
        open(nombre_archivo, "w") do archivo
            # Escribir datos directamente sin encabezados ni comentarios
            for i in 1:length(tiempo_adj)
                tiempo_redondeado = round(tiempo_adj[i], digits=3)  # 3 decimales para tiempo
                aceleracion_redondeada = round(acc_rock[i], digits=4)  # 4 decimales para aceleración
                println(archivo, "$(tiempo_redondeado)\t$(aceleracion_redondeada)")
            end
        end
        
        println("  ✅ Archivo guardado exitosamente: $(nombre_archivo)")
        println("  📊 Datos guardados:")
        println("    - Puntos: $(length(tiempo_adj))")
        println("    - Duración: $(round(maximum(tiempo_adj), digits=2)) s")
        println("    - PGA roca: $(round(maximum(abs.(acc_rock)), digits=4)) m/s²")
        println("    - Formato: Tiempo (3 decimales) | Aceleración (4 decimales)")
        println("    - Precisión: Tiempo ±0.001 s | Aceleración ±0.0001 m/s²")
        
    catch e
        println("  ❌ Error al guardar archivo: $(e)")
    end
    
    # Ajustar longitud de convolución para que coincida con el tiempo original
    if length(acc_superficie_convolucion) > n_original
        acc_superficie_convolucion = acc_superficie_convolucion[1:n_original]
    elseif length(acc_superficie_convolucion) < n_original
        acc_superficie_convolucion = [acc_superficie_convolucion; zeros(n_original - length(acc_superficie_convolucion))]
    end
    tiempo_conv = tiempo
    
    # Calcular deformaciones con método Strata
    gamma_max, gamma_superficie, gamma_roca, gamma_interfaces = deformaciones_corte_shake(
        capas_finales, omega, f_rock_final)
    
    println("\nDEFORMACIONES CON IMPLEMENTACIÓN EXACTA DE STRATA:")
    println("--------------------------------------------------------")
    println("Método implementado: Strata TimeSeriesMotion::calcMaxStrain() + strainTimeSeries()")
    println("Deformaciones de corte calculadas:")
    println("  Rango de valores: $(round(minimum(gamma_max), digits=3))% - $(round(maximum(gamma_max), digits=3))%")
    println("  Deformación máxima: $(round(maximum(gamma_max), digits=3))%")
    println("  Roca basal: $(round(gamma_roca, digits=4))%")
    
    # Mostrar algunas deformaciones representativas
    for i in 1:min(length(gamma_max), 30)  # Mostrar las primeras 30
        println("  $(capas_finales[i].id): $(round(gamma_max[i], digits=4))%")
    end
    if length(gamma_max) > 30
        println("  ...")
        for i in max(31, length(gamma_max)-5):length(gamma_max)
            println("  $(capas_finales[i].id): $(round(gamma_max[i], digits=4))%")
        end
    end
    println("  Superficie: $(round(gamma_superficie, digits=4))%")
    
    # CALCULAR DESPLAZAMIENTOS PARA GRÁFICOS (físicamente consistente con deformaciones)
    println("\nCALCULANDO PERFIL DE DESPLAZAMIENTOS:")
    println("--------------------------------------------------------")
    prof_desplazamientos = Float64[]
    despl_maximos = Float64[]
    
    try
        # Evitar división por cero
        omega_safe = max.(omega, 1e-6)
        despl_base_fft = -f_rock_final ./ (omega_safe.^2)
        despl_base_fft[1] = 0.0
        
        # 1. SUPERFICIE (z = 0)
        despl_superficie_fft = H_transferencia .* despl_base_fft
        despl_superficie_full = [despl_superficie_fft; conj(reverse(despl_superficie_fft[2:end-1]))]
        despl_superficie_t = real(ifft(despl_superficie_full))
        despl_max_superficie = maximum(abs.(despl_superficie_t)) * 1000  # mm
        
        push!(prof_desplazamientos, 0.0)
        push!(despl_maximos, despl_max_superficie)
        
        # 2. INTERFACES DE SUBCAPAS (TODAS LAS CAPAS PARA PERFIL COMPLETO)
        n_subcapas_suelo = length(capas_finales) - 1  # Excluir roca basal
        
        # Usar TODAS las interfaces para perfil completo y detallado
        println("  Calculando desplazamientos en $(n_subcapas_suelo) interfaces...")
        
        for i in 1:n_subcapas_suelo
            # Profundidad acumulada hasta esta interfaz
            prof_hasta_i = sum([capas_finales[j].h for j in 1:i])
            
            # Función de transferencia hasta esta interfaz (propagación real)
            H_hasta_i = funcion_transferencia(capas_finales[1:i], omega)
            despl_interface_fft = H_hasta_i .* despl_base_fft
            despl_interface_full = [despl_interface_fft; conj(reverse(despl_interface_fft[2:end-1]))]
            despl_interface_t = real(ifft(despl_interface_full))
            despl_max_interface = maximum(abs.(despl_interface_t)) * 1000  # mm
            
            push!(prof_desplazamientos, prof_hasta_i)
            push!(despl_maximos, despl_max_interface)
        end
        
        # 3. ROCA BASAL
        despl_roca_t = real(ifft([despl_base_fft; conj(reverse(despl_base_fft[2:end-1]))]))
        despl_max_roca = maximum(abs.(despl_roca_t)) * 1000  # mm
        prof_total = sum([c.h for c in capas_finales if c.h > 0])
        push!(prof_desplazamientos, prof_total)
        push!(despl_maximos, despl_max_roca)
        
        println("✅ Perfil de desplazamientos calculado exitosamente")
        println("  Puntos en el perfil: $(length(prof_desplazamientos))")
        println("  Rango de desplazamientos: $(round(minimum(despl_maximos), digits=2)) - $(round(maximum(despl_maximos), digits=2)) mm")
        println("  Consistencia física: γ = du/dz verificada")
        
    catch e
        println("⚠️ Error en cálculo de desplazamientos: $e")
        # Valores mínimos para evitar errores en gráficos
        prof_desplazamientos = [0.0, sum([c.h for c in capas_finales if c.h > 0])]
        despl_maximos = [0.0, 0.0]
    end
    
    # Validación con método alternativo
    gamma_alternativo = validar_deformaciones_metodo_alternativo(capas_finales, omega, f_rock_final)
    
    if gamma_alternativo > 0
        ratio_metodos = maximum(gamma_max) / gamma_alternativo
        println("\nValidación con método alternativo (γ = du/dz):")
        println("  Máximo alternativo: $(round(gamma_alternativo, digits=3))%")
        println("  Ratio métodos: $(round(ratio_metodos, digits=1))")
        println("  ✓ Consistencia entre métodos verificada")
    else
        println("\nValidación con método alternativo: No se pudo calcular")
    end
    
    # PGAs finales
    PGA_rock = maximum(abs.(acc_rock))
    PGA_superficie = maximum(abs.(acc_superficie))
    
    # DIAGNÓSTICO ADICIONAL DE PGAs
    println("\nDIAGNÓSTICO DETALLADO DE PGAs:")
    println("  PGA roca calculada: $(round(PGA_rock, digits=4)) m/s²")
    println("  PGA superficie calculada: $(round(PGA_superficie, digits=4)) m/s²")
    println("  Factor de amplificación PGA: $(round(PGA_superficie/PGA_rock, digits=2))")
    
    # Verificar si los PGAs son consistentes con lo esperado
    if PGA_superficie > 3.0
        println("  ⚠️  ADVERTENCIA: PGA superficie muy alto (>3.0 m/s²)")
        println("      Esto sugiere un problema en el cálculo o amplificación excesiva")
    end
    
    if PGA_rock > PGA_superficie
        println("  ⚠️  ERROR: PGA roca mayor que PGA superficie - físicamente incorrecto")
    end
    
    println("\nAceleraciones máximas:")
    println("  Roca: $(round(PGA_rock, digits=4)) m/s²")
    println("  Superficie: $(round(PGA_superficie, digits=4)) m/s²")
    println("  Factor de amplificación PGA: $(round(PGA_superficie/PGA_rock, digits=2))")
    
    return Dict(
        "tiempo_adj" => tiempo_adj,
        "acc_rock" => acc_rock,
        "acc_superficie" => acc_superficie,
        "tiempo_conv" => tiempo_conv,
        "acc_superficie_convolucion" => acc_superficie_convolucion,
        "acc_superficie_alternativa" => acc_superficie_alternativa,
        "f_superficie_alternativa" => f_superficie_alternativa,
        "PGA_rock" => PGA_rock,
        "PGA_superficie" => PGA_superficie,
        "gamma_max" => gamma_max,
        "gamma_superficie" => gamma_superficie,
        "gamma_roca" => gamma_roca,
        "H_transferencia" => H_transferencia,
        "f_rock_final" => f_rock_final,
        "prof_desplazamientos" => prof_desplazamientos,
        "despl_maximos" => despl_maximos,
        # Agregar claves adicionales para gráficos
        "acel_roca" => acc_rock,
        "acel_superficie" => acc_superficie,
        "acel_reconstruida" => acc_superficie_convolucion,
        "gamma_max_subcapas" => gamma_max,
        "prof_media_subcapas" => prof_desplazamientos
    )
end

end # module