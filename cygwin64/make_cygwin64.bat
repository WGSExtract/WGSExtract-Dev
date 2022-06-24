REM WGS Extract v4
REM Script to create CygWin64 Base and Development releases for WGS Extract v4.
REM
REM This content was originally in a make_standalone.sh script in v3
REM
REM Copyright (c) 2020-22 Randy Harr as part of the WGS Extract v4 release
REM Released under the xxxxx license

REM Due to the v3 to v4 change from a stripped, minimal cygwin64 release to the (full) base release, we changed
REM the directory from win10tools/ to cygwin64/. We no longer make the stripped, minimal cygwin64 release here.
REM Instead, we capture the versioned Cygwin64 base installer packages to install from and do a local install.

set /a version=2
set "baseURL=https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjgTOh_PukuJhL8M0k/root/content"

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

cd %~dp0%
set cwd=%~dp0%

set release=cygwin64
set "cygwin=https://www.cygwin.com/setup-x86_64.exe"
set install=setup-x86_64.exe
set "mirror=https://cygwin.mirror.constant.com/"
set "mirenc=https%%3a%%2f%%2fcygwin.mirror.constant.com%%2f"

set "installer=%release%_installer_v%version%.zip"
set "builder=%release%_build_v%version%.zip"

rmdir %release% %mirenc% /s/q
del %installer%

REM Cygwin64 base release with libraries needed for the bioinformatic tools
%windir%\SysWOW64\curl -kL -o %install% %cygwin%
%install% --root %release% --site %mirror% --only-site --quiet-mode --no-shortcuts ^
  --no-admin --local-package-dir %cwd% --download --categories base ^
  --packages jq,p7zip,unzip,zip,libbz2-devel,libzip-devel,liblzma-devel,libdeflate-devel,zlib-devel,libncurses-devel
copy 00README-cygwin64.txt %release%
copy make_cygwin64.bat %release%
copy etc.skel.bashrc %release%
echo {"cygwin64.version": "%version%"^, "cygwin64.date": "%DATE%"^, "cygwin64.baseURL": "%baseURL%"} > %release%\version.json
powershell Compress-Archive -Path %release%,%mirenc%,%install% -DestinationPath "%installer%" -force
REM rmdir %release% %mirenc% /s/q
REM del %install%
REM pause

REM Previous BASE and libraries with additional packages to compile the bioinformatics software
REM Dependencis bring in Base; so retain from before to avoid respecifying all the additional libraries
REM %windir%\SysWOW64\curl -kL -o %install% %cygwin%
%install% --root %release% --site %mirror% --only-site --quiet-mode --no-shortcuts ^
  --no-admin --local-package-dir %cwd% --download --packages mingw64-x86_64-gcc-core,mingw64-x86_64-gcc-g++,gcc-g++,^
autoconf,make,mingw64-x86_64-zlib,mingw64-x86_64-bzip2,mingw64-x86_64-libzip,mingw64-x86_64-zziplib,mingw64-x86_64--ncurses
powershell Compress-Archive -Path %release%,%mirenc%,%install% -DestinationPath "%builder%" -force
rmdir %release% %mirenc% /s/q
del %install%
pause

REM Installing the packaged Cygwin64 base release from Zip archive
REM %windir%\SysWOW64\curl -kL -o %installer% "https://onedrive ......"
powershell Expand-Archive -LiteralPath "%installer%" -DestinationPath "." -force
%install% --root %release% --site %mirror% --only-site --quiet-mode --no-shortcuts ^
  --no-admin --local-package-dir %cwd% --local-install --categories base ^
  --packages jq,p7zip,unzip,zip,libbz2-devel,libzip-devel,liblzma-devel,libdeflate-devel,zlib-devel,libncurses-devel
cygwin64\bin\cygstart.exe cygwin64\bin\ln.exe -s /cygdrive /mnt
type etc.skel.bashrc >> cygwin64\etc\skel\.bashrc
rmdir %mirenc% /s/q
del %install%
REM del %installer%
pause

REM Installing the packaged Cygwin64 build environment from Zip archive
REM %windir%\SysWOW64\curl -kL -o %builder% "https://onedrive ......"
powershell Expand-Archive -LiteralPath "%builder%" -DestinationPath "." -force
%install% --root %release% --site %mirror% --only-site --quiet-mode --no-shortcuts ^
  --no-admin --local-package-dir %cwd% --local-install --packages mingw64-x86_64-gcc-core,mingw64-x86_64-gcc-g++,gcc-g++,^
autoconf,make,mingw64-x86_64-zlib,mingw64-x86_64-bzip2,mingw64-x86_64-libzip,mingw64-x86_64-zziplib,mingw64-x86_64--ncurses
rmdir %mirenc% /s/q
del %install%
REM del %builder%
pause

del setup.log setup.log.full

exit

REM -------------------------------------------------------------------------------------------------------------------
REM Another option is to create a versioned package list of installed packages (--packages p7zip=15.14.1-1)
REM Note: this is the original packages AND dependencies; a much larger list. So has to be broken up as too long for
REM a command line. Issues are finding repository that retains teh versions as they get older.  And to reinstall in the
REM correct order for dependencies.
cygcheck -c -d | sed -e '1,2d' -e 's/ /=/' -e 's/ //g' > /package.lst
split -n l/3 -d /package.lst /package- --additional-suffix=.lst
for file in /package-?.lst ; do paste -sd',' $file > "${file%.*}.str" ; done
7z a -tzip cygwin64_installer_v%version%.zip package-?.str
rm -f package.lst package-?.lst package-?.str

REM http://mirror.isoc.org.il/pub/cygwin/ or maybe use TimeMachine with simple list without versions based on
REM timestamp of installation: http://www.crouchingtigerhiddenfruitbat.org/Cygwin/timemachine.html

REM Back to Windows Batch
for /l %%i in (1,1,3) do (
  set /p pack=<package-%%i.str
  setup-x86_64 --root cygwin64 --site https://cygwin.mirror.constant.com --only-site --quiet-mode --no-shortcuts ^
    --no-admin --local-package-dir %cwd% --download ^
    --packages %pack%
)
setup-x86_64.exe --root cygwin64 --site https://cygwin.mirror.constant.com --only-site --quiet-mode --no-shortcuts ^
  --no-admin --local_package_dir %cwd% --local-install --categories base ^
  --packages jq,p7zip,unzip,zip,libbz2-devel,libzip-devel,liblzma-devel,libdeflate-devel,zlib-devel,libncurses-devel
REM rmdir "https%%3a%%2f%%2fcygwin.mirror.constant.com%%2f" /s/q





