@echo off
setlocal

set "CASE_NAME=%~1"
if "%CASE_NAME%"=="" (
    echo Usage: run 1.2.3@1
    exit /b 1
)

for /f "delims=@" %%A in ("%CASE_NAME%") do set "VERSION=%%A"

set "EXE=%~dp0bin\%VERSION%.exe"
if not exist "%EXE%" (
    echo Executable not found: %EXE%
    exit /b 1
)

for /f "usebackq delims=" %%A in (`powershell -NoProfile -Command "(Get-Content '%~dp0config.json' | ConvertFrom-Json).nmlDir"`) do set "NML_DIR=%%A"

if not defined OMP_NUM_THREADS set "OMP_NUM_THREADS=24"
if not defined MKL_NUM_THREADS set "MKL_NUM_THREADS=1"
if not defined OMP_DYNAMIC set "OMP_DYNAMIC=FALSE"
if not defined MKL_DYNAMIC set "MKL_DYNAMIC=FALSE"
if not defined OMP_PROC_BIND set "OMP_PROC_BIND=close"
if not defined OMP_PLACES set "OMP_PLACES=cores"
if not defined KMP_BLOCKTIME set "KMP_BLOCKTIME=0"

"%EXE%" "%NML_DIR%%CASE_NAME%.nml"
