@echo off

rem *** Edit according to your configuration: 
rem ***     1. where is your cygwin root directory (CYGWIN)?
rem ***     2. where is your xcrysden root directory (XCRYSDEN_TOPDIR)?

set CYGWIN=C:\cygwin
set XCRYSDEN_TOPDIR=C:\cygwin\home\tone\src\XCrySDen-1.5.18-src

%CYGWIN%\bin\bash -l %XCRYSDEN_TOPDIR%\xcrysden