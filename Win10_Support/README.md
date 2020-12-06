Scripts and files to support creating binary blobs for Win10 installer.  Similar to making blobs available via a package manager on *Nix systems. 

In the first year of releases, there was a lot of massaging and tuning that had to be done.  Since then, the tools compile directly without modification.

This utilizes a Cygwin64 environment with all the development tools available.

Until this page is expanded with content, rely on reading the main make_standalone.sh file for more information.

For now, we include the processing for three distinct and separate blobs:
(1) Win10Tools blobs (binaries of bioinformatics tools in the CygWin64 stand-alone and redistributable environment
(2) Python standalone release (based on WinPython)
(3) Bioinformatic tools in JAR (Java-based) format but not the portable, standalone Java Runtime release itself

Although Python and Java may seem like they should be similar in form and release, JavaVM's tend to be available and so the installer can check and handle.  Python has the odd issue of being rarely available and Microsoft and Apple have decided to put "fake" executables in their releases which are just advertisements to install their supported "real" releases.  To get around the issues od all this in the installer, and to better assure the environment, we simply carry Python along.  Java JAR distribution is better understood and support. Python "cache" files are not.

If BioConda would extend their releases to include Win10, this would not be necessary.
