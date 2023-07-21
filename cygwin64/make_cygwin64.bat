@echo off
REM WGS Extract v5
REM Script to create CygWin64 Base and Development releases for WGS Extract v4.
REM
REM This content was originally in a make_standalone.sh script in v3
REM
REM Copyright (c) 2020-23 Randy Harr as part of the WGS Extract v4 release
REM Released under the xxxxx license

REM Due to the v3 to v4 change from a stripped, minimal cygwin64 release to the (full) base release, we changed
REM the directory from win10tools/ to cygwin64/. We no longer make the stripped, minimal cygwin64 release here.
REM Instead, we capture the versioned Cygwin64 base installer packages to install from and do a local install.

set /a version=4
for /f "delims=" %%# in ('powershell get-date -format "{dd-MMM-yyyy}"') do @set curdate=%%#
set reldir=cygwin64

set "cygwinZIP=%reldir%_v%version%_installer.zip"
set "buildZIP=%reldir%_v%version%_builder.zip"

set "baseURL=https://get.wgse.io/%cygwinZIP%"
set "buildURL=https://get.wgse.io/%buildZIP%"


REM -------------------------------------------------------------------------------------------------------------------
REM We start by creating a cygwin64 base with a few needed libraries and tools. This will be packaged and installed
REM  at each user site.  We find it simpler to package our captured "version" libraries than rely on a CygWin64
REM  site to have the versions we need.  Thus the solution is to capture a static version, downloaded Cygwin64
REM  libraries and installers, redistribute them via the WGSE server, and then locally "install" them at each site.
REM  Pain is we have to redistribute. But can now guarantee a specific version as used to create the bioinformatics
REM  tools. This is much smaller and better to redistribute the uninstalled libraries than redistributing an already
REM  installed system. Get true, local install scripts to run on the local system.

REM Need 7zip installed somewhere on Windows
REM set 7z=C:\Thumb\PortableApps\7-ZipPortable\App\7-Zip64\7zG.exe

cd /d %~dp0%
set cwd=%~dp0%

set "cygbin=%reldir%\bin"
set "cygsetURL=https://www.cygwin.com/setup-x86_64.exe"
set cygsetup=setup-x86_64.exe
set "mirror=https://cygwin.mirror.constant.com/"
set "mirenc=https%%3a%%2f%%2fcygwin.mirror.constant.com%%2f"

set "curlx=%windir%\SysWOW64\curl -kL"

REM Start with a clean slate
rmdir %reldir% %mirenc% /s/q 2>nul
del %cygwinZIP% %buildZIP% 2>nul

echo ================================================================================================================
echo Retrieving cygsetup
echo(
%curlx% -o %cygsetup% %cygsetURL%

echo ================================================================================================================
echo Creating Cygwin64 BASE release with needed libraries for bioinformatic tools (%cygwinZIP%)
echo(
%cygsetup% --root %reldir% --site %mirror% --only-site --quiet-mode --no-shortcuts --arch x86_64 ^
  --no-admin --local-package-dir %cwd% --download --categories base ^
  --packages jq,p7zip,unzip,zip,libbz2-devel,libzip-devel,liblzma-devel,libdeflate-devel,zlib-devel,libncurses-devel,^
libcurl-devel,libssl-devel > setup.txt
copy 00README-cygwin64.txt %reldir% >nul
copy make_cygwin64.bat %reldir% >nul
copy etc.skel.bashrc %reldir% >nul
echo { > %reldir%\%reldir%.json
echo   "%reldir%": { >> %reldir%\%reldir%.json
echo     "version": %version%, >> %reldir%\%reldir%.json
echo     "date": "%curdate%", >> %reldir%\%reldir%.json
echo     "URL": "%baseURL%", >> %reldir%\%reldir%.json
echo     "buildURL": "%buildURL%" >> %reldir%\%reldir%.json
echo   } >> %reldir%\%reldir%.json
echo } >> %reldir%\%reldir%.json
REM echo {"cygwin64.version": %version%^, "cygwin64.date": "%curdate%"^, "cygwin64.baseURL": "%baseURL%"^, "cygwin64.buildURL": "%buildURL%"} > %reldir%\%reldir%.json
rename %mirenc% mirror
powershell Compress-Archive -Path %reldir%,mirror,%cygsetup% -DestinationPath "%cygwinZIP%" -force
rename mirror %mirenc%
REM rmdir %reldir% %mirenc% /s/q
REM del %cygsetup%
REM pause

echo ================================================================================================================
echo Creating Cygwin64 BUILD release with additional libraries for compiling (%buildZIP%)
echo(
REM Previous BASE and libraries with additional packages to compile the bioinformatics software
REM Dependencis bring in Base; so retain from before to avoid respecifying all the additional libraries
REM %curlx% -o %cygsetup% %cygsetURL%
%cygsetup% --root %reldir% --site %mirror% --only-site --quiet-mode --no-shortcuts --arch x86_64 ^
  --no-admin --local-package-dir %cwd% --download --packages mingw64-x86_64-gcc-core,mingw64-x86_64-gcc-g++,gcc-g++,^
autoconf,make,mingw64-x86_64-zlib,mingw64-x86_64-bzip2,mingw64-x86_64-libzip,mingw64-x86_64-zziplib,^
mingw64-x86_64-ncurses,mingw64-x86_64-curl >> setup.txt
REM To bring in cython and tools dependent on it
%cygsetup% --root %reldir% --site %mirror% --only-site --quiet-mode --no-shortcuts --arch x86_64 ^
  --no-admin --local-package-dir %cwd% --download --packages mingw64-x86_64-gcc-fortran,^
python39,pyhton39-cython,python39-numpy,mingw64-x86_64-blas,mingw64-x86_64-cblas,libopenblas  >> setup.txt
rename %mirenc% mirror
powershell Compress-Archive -Path %reldir%,mirror,%cygsetup% -DestinationPath "%buildZIP%" -force
rename mirror %mirenc%
rmdir %reldir% %mirenc% /s/q
del %cygsetup%

echo Finished creating two ZIP files: Cygwin64 BASE installer and Build release.
echo(

echo ================================================================================================================
echo Installing Cygwin64 BASE release with needed libraries for bioinformatic tools (%cygwinZIP%)
echo(
REM %curlx% -o %cygwinZIP% "https://onedrive ......"
powershell Expand-Archive -LiteralPath "%cygwinZIP%" -DestinationPath "." -force
%cygsetup% --root %reldir% --site mirror --only-site --quiet-mode --no-shortcuts --arch x86_64 ^
  --no-admin --local-package-dir %cwd% --local-install --categories base ^
  --packages jq,p7zip,unzip,zip,libbz2-devel,libzip-devel,liblzma-devel,libdeflate-devel,zlib-devel,libncurses-devel,^
libcurl-devel,libssl-devel >> setup.txt
%cygbin%\cygstart.exe %cygbin%\ln.exe -s /cygdrive /mnt
type etc.skel.bashrc >> cygwin64\etc\skel\.bashrc
REM rmdir mirror /s/q
REM del %cygsetup% %cygwinZIP%

echo Finished installing Cygwin64 BASE
echo(

echo ================================================================================================================
echo Installing Cygwin64 BUILD release overlayed onto BASE with compilation tools (%buildZIP%)
echo(
REM %curlx% -o %buildZIP% "https://onedrive ......"
powershell Expand-Archive -LiteralPath "%buildZIP%" -DestinationPath "." -force
%cygsetup% --root %reldir% --site mirror --only-site --quiet-mode --no-shortcuts --arch x86_64^
  --no-admin --local-package-dir %cwd% --local-install --packages mingw64-x86_64-gcc-core,mingw64-x86_64-gcc-g++,^
gcc-g++,autoconf,make,mingw64-x86_64-zlib,mingw64-x86_64-bzip2,mingw64-x86_64-libzip,mingw64-x86_64-zziplib,^
mingw64-x86_64-ncurses,mingw64-x86_64-curl,mingw64-x86_64-gcc-fortran,^
python39,pyhton39-cython,python39-numpy,mingw64-x86_64-blas,mingw64-x86_64-cblas,libopenblas >> setup.txt
rmdir mirror /s/q
del %cygsetup%
REM del %buildZIP%

echo Finished installing Cygwin64 Build. See setup.txt/log for the complete log of commands.
pause
echo(

del setup.log setup.log.full setup.txt

goto:EOF
REM --- exit before getting here ----

REM -------------------------------------------------------------------------------------------------------------------
REM Another option is to create a versioned package list of installed packages (--packages p7zip=15.14.1-1)
REM Note: this is the original packages AND dependencies; a much larger list. So has to be broken up as too long for
REM a command line. Issues are finding repository that retains the versions as they get older.  And to reinstall in the
REM correct order for dependencies. But saves us redistributing the gathered packages.
REM draft code below is BASH version; need to change to cmd.exe version
cygcheck -c -d | sed -e '1,2d' -e 's/ /=/' -e 's/ //g' > /package.lst
split -n l/3 -d /package.lst /package- --additional-suffix=.lst
for file in /package-?.lst ; do paste -sd',' $file > "${file%.*}.str" ; done
7z a -tzip %reldir%_installer_v%version%.zip package-?.str
rm -f package.lst package-?.lst package-?.str

REM http://mirror.isoc.org.il/pub/cygwin/ or maybe use TimeMachine with simple list without versions based on
REM timestamp of installation: http://www.crouchingtigerhiddenfruitbat.org/Cygwin/timemachine.html

REM Back to Windows Batch
for /l %%i in (1,1,3) do (
  set /p pack=<package-%%i.str
  setup-x86_64 --root %reldir% --site https://cygwin.mirror.constant.com --only-site --quiet-mode --no-shortcuts ^
    --no-admin --local-package-dir %cwd% --download ^
    --packages %pack%
)
setup-x86_64.exe --root %reldir% --site https://cygwin.mirror.constant.com --only-site --quiet-mode --no-shortcuts ^
  --no-admin --local_package_dir %cwd% --local-install --categories base ^
  --packages jq,p7zip,unzip,zip,libbz2-devel,libzip-devel,liblzma-devel,libdeflate-devel,zlib-devel,libncurses-devel,libcurl-devel
REM rmdir "https%%3a%%2f%%2fcygwin.mirror.constant.com%%2f" /s/q





