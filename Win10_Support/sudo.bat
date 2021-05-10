@echo off & setlocal EnableExtensions DisableDelayedExpansion
REM Setup to check and only run elevated; equivaent of running whole script as sudo
REM from https://stackoverflow.com/questions/54658352/passing-quoted-arguments-from-batch-file-to-powershell-start-self-elevation/54701990#54701990
set "Args=%*"
net file 1>nul 2>&1 || (powershell -ex unrestricted -Command Start-Process -Verb RunAs 
-FilePath '%comspec%' --% -ArgumentList '/c %~f0 %Args:"=\""%'
  goto :eof)

REM Code is now running elevated
call %*