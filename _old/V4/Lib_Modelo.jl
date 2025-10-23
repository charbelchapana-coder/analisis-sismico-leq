"""
Lib_Modelo.jl - Librería para procesamiento del modelo de suelo
Etapa 2: Definición y procesamiento del modelo de suelo
"""

module Lib_Modelo

using LinearAlgebra
using ..Lib_TransferenciaCore: capa_suelo, parametros_degradacion

export subdividir_capas, periodo_fundamental_sitio

"""
    subdividir_capas(capas::Vector{capa_suelo}, espesor_maximo::Float64)

Función principal para subdividir capas con análisis de espesores optimizado.
Implementa la Etapa 2 del análisis con subdivisión inteligente de capas.
"""
function subdividir_capas(capas::Vector{capa_suelo}, espesor_maximo::Float64)
    println("\n2. PROCESAMIENTO DEL MODELO DE SUELO")
    println("====================================================")
    
    # Función anidada para validación de espesor
    function validar_espesor_subcapa(espesor, espesor_max)
        if espesor <= 0
            throw(ArgumentError("El espesor debe ser positivo"))
        end
        return espesor <= espesor_max
    end
    
    # Función anidada para crear subcapa
    function crear_subcapa(capa_original, indice_subcapa, espesor_subcapa)
        id_subcapa = capa_original.id * "_$(indice_subcapa)"
        return capa_suelo(
            capa_original.ρ,
            capa_original.Vs_inicial,
            capa_original.ζ_inicial,
            espesor_subcapa,
            capa_original.params_degradacion,
            id_subcapa
        )
    end
    
    # Función anidada para análisis de subdivisión
    function analizar_subdivision(capa, espesor_max)
        if capa.h <= 0  # Roca basal
            return 1, [capa.h]
        end
        
        n_subdivisiones = ceil(Int, capa.h / espesor_max)
        espesor_subcapa = capa.h / n_subdivisiones
        espesores = fill(espesor_subcapa, n_subdivisiones)
        
        return n_subdivisiones, espesores
    end
    
    subcapas = capa_suelo[]
    
    for (i, capa) in enumerate(capas)
        if capa.h <= 0  # Roca basal - no se subdivide
            push!(subcapas, capa)
            println("   Capa $(i): $(capa.id) - Sin subdivisión (roca basal)")
        else
            n_subdiv, espesores = analizar_subdivision(capa, espesor_maximo)
            
            if n_subdiv == 1
                push!(subcapas, capa)
                println("   Capa $(i): $(capa.id) - Sin subdivisión (h = $(round(capa.h, digits=1)) m)")
            else
                for j in 1:n_subdiv
                    subcapa = crear_subcapa(capa, j, espesores[j])
                    push!(subcapas, subcapa)
                end
                println("   Capa $(i): $(capa.id) - Subdividida en $(n_subdiv) subcapas de $(round(espesores[1], digits=1)) m")
            end
        end
    end
    
    println("   Total de subcapas: $(length(subcapas))")
    return subcapas
end

"""
    periodo_fundamental_sitio(capas::Vector{capa_suelo})

Calcula el período fundamental del sitio usando el método de cuarto de longitud de onda.
"""
function periodo_fundamental_sitio(capas::Vector{capa_suelo})
    
    # Función anidada para validación de propiedades
    function validar_propiedades(capa)
        if capa.Vs_actual <= 0
            throw(ArgumentError("Velocidad Vs debe ser positiva para $(capa.id)"))
        end
        if capa.h < 0
            throw(ArgumentError("Espesor debe ser no negativo para $(capa.id)"))
        end
    end
    
    # Función anidada para cálculo de tiempo de viaje
    function calcular_tiempo_viaje(capas_validas)
        tiempo_total = 0.0
        for capa in capas_validas
            if capa.h > 0  # Excluir roca basal
                validar_propiedades(capa)
                tiempo_total += 4 * capa.h / capa.Vs_actual
            end
        end
        return tiempo_total
    end
    
    capas_suelo = filter(c -> c.h > 0, capas)  # Excluir roca basal
    
    if isempty(capas_suelo)
        println("   Advertencia: No hay capas de suelo, usando período mínimo")
        return 0.1
    end
    
    tiempo_viaje = calcular_tiempo_viaje(capas_suelo)
    altura_total = sum(c.h for c in capas_suelo)
    
    println("   Altura total del depósito: $(round(altura_total, digits=1)) m")
    println("   Período fundamental: $(round(tiempo_viaje, digits=3)) s")
    
    return tiempo_viaje
end

end # module