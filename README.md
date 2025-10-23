# 🌍 Análisis Sísmico con Método Lineal Equivalente

[![Julia](https://img.shields.io/badge/Julia-1.6+-blue.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/usuario/analisis-sismico-leq.svg)](https://github.com/usuario/analisis-sismico-leq/stargazers)

Un programa completo en Julia para análisis sísmico de suelos utilizando el método lineal equivalente (LEQ) con arquitectura modular optimizada.

> 🎯 **Ideal para:** Ingenieros geotécnicos, sismólogos, investigadores en dinámica de suelos y estudiantes de posgrado.

## 📋 Descripción

Este proyecto implementa un análisis sísmico completo que permite:

- **Análisis de respuesta sísmica** de perfiles de suelo estratificado
- **Método lineal equivalente** con iteraciones de convergencia
- **Deconvolución sísmica** para obtener la señal en roca basal
- **Cálculo de espectros de respuesta** y funciones de amplificación
- **Análisis de deformaciones de corte** en cada capa
- **Comportamiento no lineal** de suelos blandos y lineal de suelos duros/rocas

## 🏗️ Arquitectura Modular

El proyecto está organizado en **9 librerías especializadas** que implementan las diferentes etapas del análisis:

```
script_LEQ.jl                 # Script principal
├── Lib_TransferenciaCore.jl  # Funciones auxiliares compartidas
├── Lib_Sismico.jl           # Etapa 1: Procesamiento de datos sísmicos
├── Lib_Modelo.jl            # Etapa 2: Modelo de suelo y subdivisión
├── Lib_Preparacion.jl       # Etapa 3: Preparación del análisis
├── Lib_AnalisisLEQ.jl       # Etapa 4: Análisis lineal equivalente
├── Lib_Movimientos.jl       # Etapa 5: Cálculo de movimientos
├── Lib_Espectros.jl         # Etapa 6: Espectros de respuesta
├── Lib_Amplificacion.jl     # Etapa 7: Funciones de amplificación
├── Lib_Graficos.jl          # Etapa 8: Generación de gráficos
└── Lib_Resumen.jl          # Etapa 9: Resumen final
```

## 🛠️ Instalación y Requisitos

### Dependencias de Julia

```julia
using Pkg
Pkg.add(["FFTW", "Plots", "DelimitedFiles", "Statistics", "LinearAlgebra", "Printf"])
```

### Versiones Recomendadas
- Julia ≥ 1.6
- FFTW.jl ≥ 1.4
- Plots.jl ≥ 1.25

## 🚀 Uso Rápido

### 1. Preparar Datos de Entrada

Crea un archivo de texto con el registro sísmico en formato de dos columnas:
```
tiempo[s]    aceleración[m/s²]
0.000        0.1499
0.005        0.1671
0.010        0.1230
...
```

### 2. Configurar el Análisis

Edita el archivo `script_LEQ.jl` en la sección **CONFIGURACIÓN PRINCIPAL**:

```julia
# 1. CONFIGURACIÓN DEL SISMO
archivo_sismo = "tu_sismo.txt"
factor_escala = 1.0

# 2. CONFIGURACIÓN DE DISCRETIZACIÓN
espesor_maximo_subcapas = 1.0  # metros

# 3. CONFIGURACIÓN DEL ANÁLISIS
max_iter = 200
tolerancia = 0.01  # 1%
```

### 3. Definir el Modelo de Suelo

Modifica la sección **DEFINICIÓN DEL MODELO DE SUELO**:

```julia
capas = [
    # Suelo blando (comportamiento no lineal)
    capa_suelo(1400.0, 500.0, 0.05, 50.0, params_diatomaceo, "Arcilla blanda"),
    
    # Suelo duro (comportamiento lineal)
    capa_suelo(2000.0, 850.0, 0.02, 100.0, params_lineal_brecha, "Arena densa"),
    
    # Roca basal (espesor = 0.0)
    capa_suelo(2500.0, 2050.0, 0.02, 0.0, params_lineal_roca, "Roca basal")
]
```

### 4. Ejecutar el Análisis

```julia
julia script_LEQ.jl
```

## 📊 Resultados

El programa genera automáticamente:

### Archivos de Salida
- **`analisis_LEQ.png`**: Gráficos completos del análisis
- **`deconvolucion_roca.txt`**: Señal deconvolucionada en roca basal

### Gráficos Incluidos
1. **Historial de convergencia** del método lineal equivalente
2. **Perfiles de velocidades y amortiguamiento** (inicial vs final)
3. **Series de tiempo** (superficie, profundidades, roca)
4. **Espectros de respuesta** por profundidad
5. **Función de amplificación** del sitio
6. **Deformaciones de corte máximas** por capa

### Información en Consola
- Diagnóstico de la señal sísmica de entrada
- Convergencia del análisis iterativo
- PGA y PSA en diferentes profundidades
- Períodos fundamentales del sitio
- Resumen de deformaciones máximas

## ⚙️ Configuración Avanzada

### Parámetros de Comportamiento No Lineal

Para suelos blandos que muestran degradación con la deformación:

```julia
params_no_lineal = parametros_degradacion(
    1.4,        # a: parámetro de degradación
    0.14,       # γ_ref: deformación de referencia
    -0.15,      # b1: coeficiente cuadrático de amortiguamiento
    0.27,       # b2: coeficiente lineal de amortiguamiento
    0.00773,    # b3: amortiguamiento mínimo
    "equivalente"
)
```

### Parámetros de Comportamiento Lineal

Para suelos duros y rocas:

```julia
params_lineal = parametros_degradacion(
    10000,      # a: sin degradación
    10000,      # γ_ref: sin degradación
    0.0,        # b1: sin variación cuadrática
    0.0,        # b2: sin variación lineal
    0.02,       # b3: amortiguamiento constante (2%)
    "lineal"
)
```

### Formatos de Archivo Soportados

El programa detecta automáticamente separadores:
- **Tabulación** (recomendado)
- **Coma** (CSV)
- **Punto y coma**
- **Espacio**

## 📁 Estructura de Archivos

```
Proyecto/
├── README.md
├── script_LEQ.jl                 # Script principal
├── Lib_*.jl                      # 9 librerías especializadas
├── sismo_*.txt                   # Archivos de entrada
├── deconvolucion_*.txt           # Resultados de deconvolución
├── analisis_*.png                # Gráficos de resultados
└── _old/                         # Versiones anteriores
    ├── V1/, V2/, V3/, V4/        # Historial de versiones
    └── Lineal_Elastico/          # Análisis lineal elástico
```

## 🔬 Metodología

### Método Lineal Equivalente

1. **Análisis elástico inicial** con propiedades iniciales
2. **Iteraciones de convergencia**:
   - Cálculo de deformaciones de corte
   - Actualización de módulos y amortiguamiento
   - Nuevo análisis con propiedades actualizadas
   - Verificación de convergencia
3. **Convergencia** cuando los cambios son < tolerancia especificada

### Deconvolución Sísmica

- **Convolución iterativa** para obtener la señal en roca
- **Función de transferencia** del sitio
- **Validación** mediante reconvolución

### Cálculo de Espectros

- **Espectros de aceleración** para diferentes amortiguamientos
- **Períodos** desde 0.01s hasta 10s
- **Múltiples profundidades** simultáneamente

## 🎯 Casos de Uso

### Ingeniería Sísmica
- Análisis de amplificación sísmica de sitio
- Estudios de microzonificación sísmica
- Evaluación de efectos de sitio

### Geotecnia Sísmica
- Análisis de licuefacción
- Evaluación de estabilidad de taludes
- Diseño de cimentaciones sísmicas

### Investigación
- Caracterización dinámica de suelos
- Validación de modelos constitutivos
- Estudios paramétricos

## ⚠️ Consideraciones Importantes

### Limitaciones del Método
- **Válido para deformaciones moderadas** (< 1-2%)
- **Comportamiento equivalente lineal** (no captura loops de histéresis)
- **Análisis 1D** (propagación vertical de ondas SH)

### Recomendaciones de Uso
- **Discretización**: subcapas ≤ λ/5 (donde λ = longitud de onda mínima)
- **Convergencia**: tolerancia típica 1-5%
- **Iteraciones**: máximo 50-200 según complejidad del modelo

## 🔧 Solución de Problemas

### Problemas Comunes

**Error al leer archivo sísmico**
```
- Verificar formato de columnas
- Comprobar separadores
- Validar que no haya caracteres especiales
```

**No convergencia del análisis**
```
- Aumentar número máximo de iteraciones
- Relajar tolerancia de convergencia
- Revisar parámetros de degradación
```

**Resultados físicamente inconsistentes**
```
- Verificar orden de capas (duras sobre blandas)
- Comprobar valores de Vs y densidad
- Validar parámetros de comportamiento no lineal
```

## 📚 Referencias

- Kramer, S.L. (1996). *Geotechnical Earthquake Engineering*
- Idriss, I.M. & Sun, J.I. (1992). *SHAKE91: A Computer Program for Conducting Equivalent Linear Seismic Response Analyses of Horizontally Layered Soil Deposits*
- Pystrata (https://github.com/arkottke/pystrata)
- Hashash, Y.M.A. et al. (2016). *DEEPSOIL 7.0, User Manual*

## 👨‍💻 Autor

Desarrollado para análisis sísmico geotécnico con Julia.

## 📄 Licencia

Este proyecto es de código abierto para uso académico y profesional.
Este proyecto no se hace responsable de algún uso no debido, ni tampoco puede hacerse responsable por resultados no adecuados para las distintas prácticas presentes.

---

💡 **Sugerencia**: Para obtener mejores resultados, calibra los parámetros de comportamiento no lineal con ensayos de laboratorio específicos del sitio.