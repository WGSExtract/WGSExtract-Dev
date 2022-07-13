@echo off
REM
REM  set_WGSEpath.bat  %1=Path for Installation (will add win10tools\bin)
REM 
REM  Set PATH and WGSE variables in both local environment (immediate) and System registry
REM   Due to System Registry setting (setx), need elevated privileges for those commands
REM   But due to escaping spaces and % variables in PATH, cannot use sudo.bat we created to handle
REM   just elevating a single command.  So elevate this whole script.
REM
REM Copyright (c) 2021, Randy Harr
REM  part of the WGS Extract v3 Win10 Installer release (see that release for license details)
REM

if "%1" == "" (
  echo Usage: Call %~nx0 PATH_ADDITION
  echo   This batch file will prepend the System Registry PATH entry with %%WGSEBIN%%
  echo   and set a System Registry variable named %%WGSEBIN%% with the PATH_ADDITION.
  echo   note: setting the System Registry requires elevated privileges; requested if not set.
  exit /B
) else (
  set "Args=%1"
)

REM Setup to check and run elevated; equivaent of running whole script as sudo in Unix/Linux
REM from https://stackoverflow.com/questions/54658352/passing-quoted-arguments-from-batch-file-to-powershell-start-self-elevation/54701990#54701990
net file 1>nul 2>&1 || (powershell -ex unrestricted -Command Start-Process -Verb RunAs ^
-FilePath '%comspec%' -ArgumentList '/c %~f0 %Args:"=\""%'
  echo Elevated copy of %~nx0 is now running in a separate CMD window
  goto :eof)

REM ========================= Elevated Code Below ============================
cd "%~dp0"
REM Set local variable for path (does not extend beyond this shell)
set "WGSEBIN=%~dp0win10tools\bin"

REM Set System variable to help reduce total System path variable and make future update easier
REM *** requires Elevated privilege ***
setx WGSEBIN "%WGSEBIN%" /m

REM Add WGSEBIN variable to System PATH variable; if not already there (if-else expands the PATH variable)
if ";%PATH:win10tools=%;" == ";%PATH%;" (
  REM Even though PATH has %%WGSEBIN%% set; it gets evaluated here so have to test its value

  REM Save the current PATH just in case User needs to recreate their original one
  echo "%PATH%" > "%USERPROFILE%\WGSE_Win10SavedPath.txt"

  REM Set local path variable -- not really needed as dropped when this script exits
  set "PATH=%WGSEBIN%;%PATH%"

  REM Creating separate System and User Path variables from registry entries; need to set only one
  REM  from https://stackoverflow.com/questions/43934167/how-to-get-only-the-user-path-variable

  setlocal EnableExtensions EnableDelayedExpansion
  set "SystemPath="
  for /F "skip=2 tokens=1,2*" %%N in ('%SystemRoot%\System32\reg.exe query "HKLM\System\CurrentControlSet\Control\Session Manager\Environment" /v "Path" 2^>nul' ) do (
    if /I "%%N" == "Path" (
      set "SystemPath=%%P"
      goto SetPath
    )
  )

  :SetPath
  REM Replace all two consecutive semicolons by a single semicolon. Then remove any trailing semicolon
  set "SystemPath=!SystemPath:;;=;!" 
  if "!SystemPath:~-1!" == ";" ( set "SystemPath=!SystemPath:~0,-1!" )

  REM We prepend the System Path because Windows has some false executables setup to advertise their
  REM  services of the Windows Store -- python, bash, etc.  So we need our copies found first.
  REM We use a variable (WGSEBIN) stored on the PATH so it can be easily updated later.

  REM *** requires Elevated privilege***
  setx PATH "%%WGSEBIN%%;!SystemPath!" /m

  endLocal
)
