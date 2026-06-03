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

"%EXE%" "%NML_DIR%%CASE_NAME%.nml"
