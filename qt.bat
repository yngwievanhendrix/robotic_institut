@echo off

rem Select your Visual Studio version
rem set QMAKESPEC=win32-msvc.net
rem set QMAKESPEC=win32-msvc2005
set QMAKESPEC=win32-msvc2008

IF NOT DEFINED BUILDTYPE (SET BUILDTYPE=debug)

SET VTK_DIR=C:\Programme\VTK
SET QTDIR=C:\Qt\4.7.0
SET CMAKE_DIR=C:\Program Files\CMake 2.8\bin

SET PATH=%VTK_DIR%\bin\%BUILDTYPE%;%QTDIR%\bin;%CMAKE_DIR%;%PATH%

qmake -tp vc "CONFIG+=flat" vtk_Reader.pro
pause