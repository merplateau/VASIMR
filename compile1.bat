@echo off
setlocal

set "V="
set "DEBUG=0"

rem ---- parse args: -d/--debug enables debug build; the other arg is the version ----
:parse
if "%~1"=="" goto after_parse
if /i "%~1"=="-d" (
    set "DEBUG=1"
    shift
    goto parse
)
if /i "%~1"=="--debug" (
    set "DEBUG=1"
    shift
    goto parse
)
set "V=%~1"
shift
goto parse
:after_parse

if "%DEBUG%"=="1" (
    rem Do not add /fpe-all:0; MKL/PARDISO internals may trigger floating point exceptions.
    rem /check:noarg_temp_created suppresses harmless temporary-array diagnostics.
    set "DBGFLAGS=/Od /debug:full /check:all /check:noarg_temp_created"
    echo [compile1.bat] DEBUG build: %DBGFLAGS%
) else (
    set "DBGFLAGS="
)

set "VS2022INSTALLDIR=C:\Program Files\Microsoft Visual Studio\18\Community"

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

cd /d "%~dp0src"

rem /traceback is always enabled so runtime errors print source file and line information.
ifx /traceback %DBGFLAGS% hypic.f90 fun_b0.f90 fun_record_display.f90 fun_grid.f90 fun_ini.f90 fun_antenna_irf.f90 fun_particles.f90 fun_mcc.f90 fun_fdfd.f90 /Qmkl /object:..\tmp\ /module:..\tmp\ -o ..\bin\%V%.exe /link /STACK:500000000

endlocal
