@echo off
setlocal

set "V=%~1"

set "VS2022INSTALLDIR=D:\Program Files\Microsoft Visual Studio\2022\Community"

call "D:\Program Files (x86)\Intel\oneAPI\setvars.bat"

cd /d "%~dp0src"

ifx hypic.f90 fun_b0.f90 fun_record_display.f90 fun_grid.f90 fun_ini.f90 fun_antenna_irf.f90 fun_particles.f90 fun_mcc.f90 fun_fdfd.f90 /Qmkl /object:..\tmp\ /module:..\tmp\ -o ..\bin\%V%.exe /link /STACK:500000000

endlocal