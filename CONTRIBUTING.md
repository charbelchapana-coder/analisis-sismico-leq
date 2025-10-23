# Contribuir al Proyecto

¡Gracias por tu interés en contribuir al análisis sísmico LEQ! 🙏

## 🚀 Formas de Contribuir

### 🐛 Reportar Bugs
- Usar [Issues](../../issues) para reportar problemas
- Incluir información del sistema (Julia version, OS)
- Proporcionar datos de entrada que causen el error
- Describir comportamiento esperado vs actual

### ✨ Sugerir Mejoras
- Nuevas funcionalidades en el análisis sísmico
- Mejoras en la documentación
- Optimizaciones de rendimiento
- Nuevos tipos de gráficos o reportes

### 🔧 Enviar Pull Requests

#### Para Cambios Menores
- Correcciones de documentación
- Fixes de bugs pequeños
- Mejoras de código sin cambios funcionales

#### Para Cambios Mayores
1. **Crear Issue primero** para discutir el cambio
2. **Fork** el repositorio
3. **Crear rama**: `git checkout -b feature/nueva-funcionalidad`
4. **Implementar** cambios con tests
5. **Commit**: `git commit -m "Add: nueva funcionalidad para..."`
6. **Push**: `git push origin feature/nueva-funcionalidad`
7. **Crear Pull Request**

## 🧪 Guías de Desarrollo

### Estructura del Código
```
src/
├── Lib_TransferenciaCore.jl  # Funciones compartidas
├── Lib_Sismico.jl           # Procesamiento sísmico
├── Lib_Modelo.jl            # Modelo de suelos
└── ...                      # Otros módulos especializados
```

### Estilo de Código
- **Nombres de funciones**: `snake_case` en español
- **Nombres de variables**: descriptivos y en español
- **Comentarios**: explicar el "por qué", no el "qué"
- **Documentación**: usar docstrings para funciones públicas

### Testing
- Ejecutar tests existentes: `julia test/runtests.jl`
- Agregar tests para nuevas funcionalidades
- Validar con datos conocidos de literatura

### Documentación
- Actualizar README.md si es necesario
- Incluir ejemplos de uso para nuevas funciones
- Mantener consistencia con el estilo existente

## 📋 Checklist para Pull Requests

- [ ] El código sigue las convenciones del proyecto
- [ ] Se agregaron tests para nueva funcionalidad
- [ ] Los tests existentes pasan
- [ ] Se actualizó la documentación
- [ ] Se probó con datos de ejemplo
- [ ] El commit message es descriptivo

## 🏷️ Versionado

Seguimos [Semantic Versioning](https://semver.org/):
- **MAJOR**: cambios incompatibles en API
- **MINOR**: nueva funcionalidad compatible
- **PATCH**: bug fixes compatibles

## 💡 Ideas de Contribución

### Fácil (ideal para empezar)
- [ ] Mejorar mensajes de error
- [ ] Agregar validación de entrada
- [ ] Optimizar gráficos existentes
- [ ] Traducir comentarios a inglés

### Intermedio
- [ ] Implementar nuevos modelos constitutivos
- [ ] Agregar análisis de sensibilidad
- [ ] Mejorar algoritmos de convergencia
- [ ] Crear interfaces gráficas

### Avanzado
- [ ] Análisis 2D/3D
- [ ] Paralelización
- [ ] Integración con otros software
- [ ] Machine learning para calibración

## 🤝 Código de Conducta

- Ser respetuoso y constructivo
- Aceptar feedback y críticas
- Ayudar a otros desarrolladores
- Mantener discusiones técnicas profesionales

## 📞 Contacto

- **Issues**: Para problemas técnicos
- **Discussions**: Para preguntas generales
- **Email**: Para temas sensibles

## 🏆 Reconocimientos

Los contribuidores aparecerán en:
- Lista de AUTHORS.md
- Release notes
- README principal

¡Esperamos tu contribución! 🚀