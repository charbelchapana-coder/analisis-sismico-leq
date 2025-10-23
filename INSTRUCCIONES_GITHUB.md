# 🚀 GUÍA PASO A PASO PARA SUBIR A GITHUB

## ANTES DE COMENZAR

### 1. Instalar Git (si no lo tienes)
- Ve a: https://git-scm.com/download/win
- Descarga e instala Git para Windows
- Durante la instalación:
  ✅ Selecciona "Git Bash Here" 
  ✅ Selecciona "Use Git from the Windows Command Prompt"
  ✅ Usa las opciones por defecto para el resto

### 2. Crear cuenta en GitHub (si no la tienes)
- Ve a: https://github.com
- Regístrate con tu email
- Verifica tu email

## PASOS PARA SUBIR TU PROYECTO

### PASO 1: Configurar Git (solo la primera vez)
Abre una terminal (CMD o PowerShell) y ejecuta:

```bash
git config --global user.name "Tu Nombre Completo"
git config --global user.email "tu.email@gmail.com"
```

### PASO 2: Crear repositorio en GitHub
1. Ve a github.com
2. Clic en "New repository" (botón verde)
3. Nombre: `analisis-sismico-leq`
4. Descripción: `Análisis sísmico con método lineal equivalente en Julia`
5. ✅ Marcar "Public"
6. ❌ NO marcar "Initialize with README"
7. Clic "Create repository"

### PASO 3: Subir tu código
Abre terminal en tu carpeta del proyecto y ejecuta ESTOS COMANDOS:

```bash
# Navegar a la carpeta del proyecto
cd "C:\Users\charb\OneDrive\Desktop\Modelos Plaxis\Deconvolucion\Julia"

# Inicializar Git
git init

# Añadir todos los archivos
git add .

# Hacer el primer commit
git commit -m "Primer commit: Análisis sísmico LEQ completo con 9 librerías modulares"

# Cambiar nombre de rama a main
git branch -M main

# Conectar con tu repositorio de GitHub (REEMPLAZA 'tu-usuario' con tu nombre de usuario)
git remote add origin https://github.com/tu-usuario/analisis-sismico-leq.git

# Subir archivos a GitHub
git push -u origin main
```

## COMANDOS FUTUROS PARA ACTUALIZAR

Cuando hagas cambios en el futuro:

```bash
# Ver qué archivos cambiaron
git status

# Añadir cambios
git add .

# Hacer commit con mensaje descriptivo
git commit -m "Descripción de los cambios realizados"

# Subir cambios a GitHub
git push
```

## VERIFICAR QUE TODO ESTÁ BIEN

1. Ve a tu repositorio en GitHub
2. Deberías ver todos tus archivos:
   - ✅ README.md con toda la documentación
   - ✅ Todas las librerías Lib_*.jl
   - ✅ script_LEQ.jl (archivo principal)
   - ✅ Project.toml con dependencias
   - ✅ LICENSE
   - ✅ .gitignore
   - ✅ examples/ con ejemplos
   - ✅ Datos de ejemplo (sismo_*.txt)

## PROBLEMAS COMUNES

### Si Git no está instalado:
- Descargar de git-scm.com
- Reiniciar terminal después de instalar

### Si pide autenticación:
- GitHub eliminó passwords simples
- Usa "Personal Access Token":
  1. GitHub → Settings → Developer settings → Personal access tokens
  2. Generate new token (classic)
  3. Selecciona scopes: repo, workflow, write:packages
  4. Usa este token como password

### Si el repositorio ya existe:
```bash
git remote -v
git remote set-url origin https://github.com/tu-usuario/analisis-sismico-leq.git
```

## 🎉 ¡LISTO!

Tu proyecto estará disponible en:
https://github.com/tu-usuario/analisis-sismico-leq

Otros podrán:
- ⭐ Ver tu código
- 📥 Clonarlo con: `git clone https://github.com/tu-usuario/analisis-sismico-leq.git`
- 🔄 Contribuir con mejoras
- 📚 Leer la documentación completa