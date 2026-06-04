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
    rem 注：不加 /fpe-all:0 —— MKL(pardiso)内部的推测式浮点会误触发该陷阱
    set "DBGFLAGS=/Od /debug:full /traceback /check:all"
    echo [compile.bat] DEBUG build: %DBGFLAGS%
) else (
    set "DBGFLAGS="
)

set "VS2022INSTALLDIR=D:\Program Files\Microsoft Visual Studio\2022\Community"

call "D:\Program Files (x86)\Intel\oneAPI\setvars.bat"

cd /d "%~dp0src"

ifx %DBGFLAGS% hypic.f90 fun_b0.f90 fun_record_display.f90 fun_grid.f90 fun_ini.f90 fun_antenna_irf.f90 fun_particles.f90 fun_mcc.f90 fun_fdfd.f90 /Qmkl /object:..\tmp\ /module:..\tmp\ -o ..\bin\%V%.exe /link /STACK:500000000

endlocal
