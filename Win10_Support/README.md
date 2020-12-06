# Win10 Support Subsystem

Scripts and files to support creating binary blobs for use by the **Win10** installer pf **WGS Extract**.  Similar to making blobs available via a package manager on *Nix systems. **Win10** implies the Microsoft Windows 10 environment. Some have seen success with use in earlier major releases but we have no focus nor support for that here. 

This utilizes a [Cygwin64](https://cygwin.com/) environment with all the development tools available.

Until this page is expanded with content, rely on reading the main make_standalone.sh file for more information.

For now, we include the processing for three distinct and separate blobs:
1. **Win10Tools** blob (binaries of bioinformatics tools in the CygWin64 stand-alone and redistributable environment (native Win10 compiled binaries of C, C++, Nim, etc sourced programs; mainly from [Sanger Institute HTSLib](http://www.htslib.org/) with others)
1. **Python standalone release** blob (based on [WinPython](https://winpython.github.io/)) along with any bioinformatic tool libraries pre-installed
1. **Bioinformatic tools in JAR** blob (Java-based);  but not the portable, standalone Java Runtime release itself (generally the [Broad Institute GATK](https://gatk.broadinstitute.org/hc/en-us) family of tools) 

Although Python and Java may seem like they should be similar in form and release, JavaVM's tend to be available and so the installer can check and handle.  Python has the odd issue of being rarely available and Microsoft and Apple have decided to put "fake" executables in their releases which are just advertisements to install their supported "real" releases.  To get around the issues od all this in the installer, and to better assure the environment, we simply carry Python along.  Java JAR distribution is better understood and support. Python "cache" files are not.

The **Python Pip library** package install may getpulled out of here.  As it is handled in the **Linux** / **MacOSX** installers and should likely be handled, therefore, in the Win10 installer in a similar fashion.

In the first year of releases, there was a lot of massaging and tuning that had to be done.  Since then, the tools compile directly without modification.

If BioConda would extend their releases to include Win10, this would not be necessary.
