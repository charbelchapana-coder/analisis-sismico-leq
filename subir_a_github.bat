@echo off
echo ========================================
echo    SUBIR PROYECTO A GITHUB - AUTOMATICO
echo ========================================
echo.

REM Verificar si Git está instalado
git --version >nul 2>&1
if errorlevel 1 (
    echo ❌ Git no está instalado
    echo 💡 Descarga Git desde: https://git-scm.com/download/win
    echo 💡 Después ejecuta este script nuevamente
    pause
    exit /b 1
)

echo ✅ Git detectado correctamente
echo.

REM Cambiar al directorio del proyecto
cd /d "C:\Users\charb\OneDrive\Desktop\Modelos Plaxis\Deconvolucion\Julia"
echo 📁 Directorio actual: %cd%
echo.

REM Verificar si ya es un repositorio Git
if exist .git (
    echo ⚠️  Ya es un repositorio Git
    echo 🔄 Actualizando archivos...
    git add .
    set /p commit_msg="💬 Mensaje del commit: "
    git commit -m "%commit_msg%"
    git push
    echo ✅ Actualización completada
) else (
    echo 🆕 Inicializando nuevo repositorio...
    
    REM Inicializar Git
    git init
    echo ✅ Repositorio inicializado
    
    REM Agregar todos los archivos
    git add .
    echo ✅ Archivos agregados
    
    REM Hacer primer commit
    git commit -m "Primer commit: Análisis sísmico LEQ completo con arquitectura modular"
    echo ✅ Primer commit realizado
    
    REM Cambiar a rama main
    git branch -M main
    echo ✅ Rama principal configurada
    
    echo.
    echo 🌐 CONFIGURAR REPOSITORIO REMOTO
    echo 💡 Ve a https://github.com y crea un nuevo repositorio llamado: analisis-sismico-leq
    echo 💡 NO inicialices con README (ya tienes uno)
    echo.
    set /p github_user="👤 Tu nombre de usuario de GitHub: "
    
    REM Agregar repositorio remoto
    git remote add origin https://github.com/%github_user%/analisis-sismico-leq.git
    echo ✅ Repositorio remoto configurado
    
    REM Subir a GitHub
    echo 🚀 Subiendo archivos a GitHub...
    git push -u origin main
    
    if errorlevel 1 (
        echo ❌ Error al subir archivos
        echo 💡 Posibles soluciones:
        echo    - Verifica tu usuario de GitHub
        echo    - Usa Personal Access Token en lugar de contraseña
        echo    - Verifica que el repositorio existe en GitHub
        pause
        exit /b 1
    )
    
    echo ✅ ¡PROYECTO SUBIDO EXITOSAMENTE!
)

echo.
echo 🎉 TU REPOSITORIO ESTÁ EN:
echo https://github.com/%github_user%/analisis-sismico-leq
echo.
echo 📋 ARCHIVOS INCLUIDOS:
echo ✅ README.md (documentación completa)
echo ✅ Todas las librerías Lib_*.jl
echo ✅ script_LEQ.jl (archivo principal)
echo ✅ Project.toml (dependencias)
echo ✅ examples/ (ejemplos de uso)
echo ✅ LICENSE y .gitignore
echo ✅ Instrucciones de contribución
echo.
echo 💡 PRÓXIMOS PASOS:
echo 1. Ve a tu repositorio en GitHub
echo 2. Agrega una descripción y temas (topics)
echo 3. Considera hacer el repositorio "starred" para promoverlo
echo.
pause