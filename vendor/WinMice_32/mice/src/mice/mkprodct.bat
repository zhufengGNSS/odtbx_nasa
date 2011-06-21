rem
rem  mkprodct.bat for mice.
rem
rem  Creates mice.dll with MS Visual C++ for use with the MATLAB-CSPICE
rem  interface. 
rem
rem  Version 1.1.0, 16-MAR-2009, (EDW) (BVS)
rem
rem     Changed .dll to .mexw32.
rem
rem  Version 1.0.0, 15-FEB-2005, (EDW) (BVS)
rem

rem
rem  Set the master MATLAB directory.
rem

set MATLAB=C:\Program Files\MATLAB\R2006a

rem
rem  Set library name.
rem 

set NAME=mice

rem
rem  Set file names based on the library name.
rem

set DLLFILE=%NAME%.mexw32
set DEFFILE=%NAME%.def
set LIBFILE=%NAME%.lib
set EXPFILE=%NAME%.exp




rem
rem  Set compile command.
rem

set cl= /c /O2 -I..\..\include -I"%MATLAB%\extern\include" -I"%MATLAB%\simulink\include" -DMATLAB_MEX_FILE -nologo -Zp8

rem
rem  Compile all .c files.
rem

for %%f in (*.c) do cl %%f 

dir /b *.obj > temp.lst

rem
rem  Set MATLAB and SPICE libraries to be linked in.
rem

set LIBMATLAB= libmx.lib libmex.lib libmat.lib
set LIBSPICE= ..\..\lib\cspice.lib

rem
rem  Link DLL.
rem

link /DLL /OUT:%DLLFILE% /DEF:%DEFFILE% /IMPLIB:%LIBFILE% /LIBPATH:"%MATLAB%\extern\lib\win32\microsoft" %LIBMATLAB% %LIBSPICE% @temp.lst 

rem
rem  Delete the intermediate files.
rem
del *.obj
del temp.lst
del %LIBFILE%
del %EXPFILE%

rem
rem  Move the DLL to the lib directory.
rem
move %DLLFILE% ..\..\lib

