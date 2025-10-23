# Ejemplos de Uso - Análisis Sísmico LEQ

Este directorio contiene ejemplos prácticos para diferentes casos de uso del programa.

## 📁 Contenido

### Datos de Ejemplo
- `sismo_ejemplo.txt` - Registro sísmico sintético para pruebas
- `modelo_suelo_basico.jl` - Configuración básica de 3 capas
- `modelo_suelo_complejo.jl` - Configuración avanzada multicapa

### Scripts de Ejemplo
- `ejemplo_basico.jl` - Análisis simple paso a paso
- `ejemplo_parametrico.jl` - Estudio paramétrico automatizado
- `ejemplo_validacion.jl` - Comparación con resultados conocidos

## 🚀 Cómo Usar

1. **Ejecutar ejemplo básico:**
```julia
julia examples/ejemplo_basico.jl
```

2. **Estudio paramétrico:**
```julia
julia examples/ejemplo_parametrico.jl
```

3. **Validación con datos reales:**
```julia
julia examples/ejemplo_validacion.jl
```

## 📊 Resultados Esperados

Cada ejemplo genera:
- Gráficos de análisis en formato PNG
- Archivos de texto con resultados numéricos
- Resumen en consola con métricas clave

## 📝 Personalización

Modifica los parámetros en cada script para adaptar a tu caso específico:
- Propiedades del sismo de entrada
- Configuración del modelo de suelo
- Tolerancias de convergencia
- Archivos de salida