An image of a cygwin64 Installer files and setup program to run a base (basic) Cygwin64 environment with
 libraries; versioned and matched to the compiled Bioinformatic Tool release.

Copyright (c) 2022. Randy Harr as part of the WGS Extract release.
See https://wgsextract.github.io/ for license details.

Created as part of the WGS Extract (https://wgsextract.github.io) project on the Windows 10/11 platform
 for the support of manipulating Whole-genome sequence files.

See the make_cygwin64.bat file for more information on what is made and how.

Starting in v4, the whole base is recreated and provided.  But still with fixed versions from when the bioinformatic
  tools were compiled.  The Cygwin64 installer is used to locally install the saved, versioned libraries.  The
  bioinformatic tools, instead of being merged / mixed into the /bin directory and other top level locations needed, are
  now put in /usr/local/.  As a result, because we cannot guarantee the PATH is setup correctly, some cygwin64 DLL's
  are copied into /usr/local/bin as well.  To use the bioinformatic tools, you need to add /usr/local/bin to your path.

The Cygwin64 base environment can be (re)created with the following (cmd.exe) commands (if using the same versions):
  \Windows\SysWOW64\curl.exe -kL -o setup-x86_64.exe "https://www.cygwin.com/setup-x86_64.exe"
  setup-x86_64.exe --root . --site https://cygwin.mirror.constant.com --only-site --quiet-mode --no-shortcuts ^
    --no-admin ---categories base --packages jq,p7zip,unzip,zip
  setup-x86_64.exe --root . --site https://cygwin.mirror.constant.com --only-site --quiet-mode --no-shortcuts ^
    --no-admin --packages ^
    libbz2-devel,libzip-devel,liblzma-devel,libdeflate-devel,zlib-devel,libncurses-devel,libcurl-devel,libssl-devel
  rd https%3a%2f%2fcygwin.mirror.constant.com%2f /s/q
  bin\cygstart.exe bin\ln.exe -s /cygdrive /mnt
  type etc.skel.bashrc >> etc\skel\.bashrc
  del setup-x86_64.exe
