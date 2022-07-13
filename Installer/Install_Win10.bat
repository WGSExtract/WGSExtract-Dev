@echo off
REM WGS Extract v3 release
REM Copyright (c) 2021 Randy Harr
REM
REM Now requires the install to download Win10 Bioinformatic tools, stand-alone, pre-built Python, and CygWin environment
REM

REM  Need to be in the installation as CWD could be anywhere if called from the GUI or Symlink
cd %~dp0%

REM Need either python  installed (to do rest of install in python) or BASH/AWK to bypass the Google Drive
REM  "cannot scan for virus" response for large blob files you try to download.  As we removed the Win10 executables
REM  from the general release, we now create a subset CYGWIN tools package that can be downloaded and bootstrap the
REM  process on Win10 systems. Actual file is win10tools-cygwin64.zip
if not exist win10tools\ (
    echo Downloading the CygWin64 environment to bootstrap the rest of the download in a Cygwin Linux BASH shell.
    echo There is often a 30 second pause by Google before the download starts and you are informed of its progress.
REM powershell -Command "Invoke-WebRequest 'https://drive.google.com/uc?&export=download&id=1lYZbFJ3eyDps7e4I_Yu4zVyZz_vjjVcd' -Outfile win10tools.zip"
REM curl -L -o win10tools.zip "https://drive.google.com/uc?&export=download&id=1lYZbFJ3eyDps7e4I_Yu4zVyZz_vjjVcd"
    curl -L -o win10tools.zip "https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjcGO3m9V1V_i-mho/root/content"
    powershell Expand-Archive -LiteralPath "win10tools.zip" -DestinationPath "." -Force
REM unzip -qo win10tools.zip
    del win10tools.zip
    echo Win10 Cygwin bootstrap installed
)

REM Will add win10tools\bin to supplied directory then prepend to PATH; sets variable WGSEBIN -- runs elevated
call set_WGSEpath.bat "%~dp0"

REM Set local PATH path variable for immediate use; registry and PATH not re-read; cannot use variable in PATH
set "WGSEBIN=%~dp0win10tools\bin"
if ";%PATH:WGSEBIN=%;" == ";%PATH%;" ( set "PATH=%WGSEBIN%;%PATH%" )

REM Note: we push the bioinformatic tools, python and similar large downloads to the Upgrade script.
REM       Upgrade is a BASH script that can handled the Google Drive large file downloads
REM       We need to make sure we use our BASH and not the Win10 supplied one or any other
echo Starting Upgrade script to finish install
win10tools\bin\bash.exe Upgrade_v2tov3.sh
if %ERRORLEVEL% NEQ 0 (
  echo Restarting Upgrade script due to a common Cygwin64 on Win10 fork error
  win10tools\bin\bash.exe Upgrade_v2tov3.sh
)

echo:
echo Congratulations!  You finished installing WGS Extract v3!
echo You can start WGS Extract v3 by clicking the Start_Win10.bat file. Make a
echo shortcut, rename it to WGSExtract, and place it in your desktop or Start
echo menu for convenience. You can re-run this installer if desired.
