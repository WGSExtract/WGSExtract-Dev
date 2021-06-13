# Win10 Support Subsystem

Scripts and files to support creating binary blobs for use by the **Win10** installer of **WGS Extract**.  Similar to making blobs available via a package manager on *Nix systems. **Win10** implies the Microsoft Windows 10 public release environment. Some have had success with earlier major releases but we have no focus nor support for that here. 

This utilizes a [Cygwin64](https://cygwin.com/) environment with all the development tools available to create the releases.

Until this page is expanded with content, rely on reading the main make_standalone.sh file for more information.

For now, we include the processing for three distinct and separate blobs:
1. **Win10Tools** blob (binaries of bioinformatics tools in the CygWin64 stand-alone and redistributable environment (native Win10 compiled binaries of C, C++, Nim, etc sourced programs; mainly from [Sanger Institute HTSLib](http://www.htslib.org/) with others)

The following used to be part of the above binary release but are now in the installer.  Pulled out separately to be built from publicly available releases.
1. **Python standalone release** blob (based on [WinPython 64bit](https://winpython.github.io/)) along with any bioinformatic tool and utility libraries pre-installed
2. **Java JVM (standalone)** blob 
3. **Bioinformatic tools in JAR form** blob (Java-based) (generally the [Broad Institute GATK](https://gatk.broadinstitute.org/hc/en-us) family of tools) 

Python has the odd issue of being rarely available and Microsoft and Apple have decided to put "fake" executables in their releases which are just advertisements to install their supported "real" releases.  The **Python Pip library** package installs are in the installer also; similar to JAR packages of Java based tools.

In the Historical reelase, there was a lot of massaging and tuning that had to be done.  Since then, the tools compile directly without modification (for the most part).

If BioConda would extend their binary releases to include Win10, this whole area would not be necessary.

Cygwin64 has poor support for Linux rmalloc and other libraries.  We had to rewrite the library so BWA would compile and parallelize on Win10 as it does on MacOS and Linux.
