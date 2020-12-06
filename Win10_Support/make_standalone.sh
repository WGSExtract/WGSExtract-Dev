#!/bin/bash

# WGS Extract Beta v3 Win10 Support  ****WORKING DRAFT -- NOT RUNNABLE YET****
#
# Script to create Win10 release of tools needed for WGS Extract install.  Since no package manager is available for Windows, 
#  we need to create pseudo-packages here that we then make available to the installers from our own repository.  The
#  bioinformatic tools are handled in one large package.  And the Python release in another.
#
#  For the bioinformatics tools, we grab the sources and create the executable release structure.  This instead of risking 
#  compiling at each site (which would require a full Cygwin64 development tool environment). Then a simple release zip
#  pack can be downloaded by the installer and made available using local scripts. Just as if using a package manager.
#  We do not install into Win10 like a normal tool release.  Just make the executables available and let the installer
#  worry about how to make them available in the local environment.
#
#  For Java JAR tools, we massage the releases to strip them to just what is needed (Jar files+); we install Java itself as
#   part of the install / release script. Most environments have a Java JVM whereas Python is more unique (see below).
#
#  For Python, we use the more portable WinPython release.  We only need the actual python directory though. And it comes as
#  a zip archive executable.  So we do some minimal processing to trim to only what we need. Win10, as MacOSX, has really messed
#  with the Python environment by putting pseudo executables in that then simply redirect to their software store.  To avoid
#  all the issues of that and getting the tested, correct version, we simply install our own local copy of Python. We simply bring
#  in the necessary packages now as well (PIP libraries; similar to Java Jar files).
#
# This script assumes it is being run in a Cygwin64 BASH shell and the Cygwin development environment is on the exe path
#  The cygwin64 development environment consists of the following main packages (and their dependencies) beyond the basic release:
#    TBD
#
# Copyright (c) 2020 Randy Harr
# Released under the xxxxx license
#  

##############################################################################################################
#  Win10 Implementations of the Bioinformatics Tools (win10tools)
#

# Make our release directory. Copy these scipts simply for documentation to the end user in case they want to remake
# mkdir win10tools; cp $0 win10tools; cp unix.tools win10tools

# Grab the latest releases that we tested with
# curl -L https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 | tar jxf
# curl -L https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz | tar jxf
# curl -L https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 | tar jxf
# curl -L https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar jxf
# curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 | tar jxf -
# curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar jxf -

# First make the various bioinformatic tools and copy the .exe to bin, man page to share/man/man1, etc
cd htslib-1.11   ; configure ; make ; make prefix=../win10tools install ; make clean ; cd ..
cd samtools-1.11 ; configure ; make ; make prefix=../win10tools install ; make clean ; cd ..
cd bcftools-1.11 ; configure ; make ; make prefix=../win10tools install ; make clean ; cd ..
cd bwa      ; make ; cp bwa.exe ../win10tools/bin ; cp bwa.1 ../win10tools/share/man/man1 ; make clean ; cd ..
cd bwa-mem2 ; make ; cp bwa-mem2* ../win10tools/bin ; make clean ; cd ..  # no man page to copy
cd minimap2 ; make ; cp minimap2.exe ../win10tools/bin ; cp minimap2.1 ../win10tools/share/man/man1 ; make clean ; cd ..

# Grab needed and useful Unix utilities from the main CygWin bin
for exe in `cat unix.tools`
do
  cp /bin/$exe win10tools/bin
done

# Check all the *.exe's in bin for dll dependencies; copy all those to the bin as well
for exe in win10tols/bin/*.exe
do
  cygcheck $exe | grep -v Randy | grep -v WINDOWS > dep.tmp
  for dll in `cat dep.tmp`
  do
    cp $dll win10tools/bin
  done
done
cp /bin/curl-config win10tools/bin
rm dep.tmp

# Set loader flag to allow big memory model for certain executables (goes from default 512MB to 2GB max)
# Note: was hoping this would help with low memory utilization in samtools, etc but does not seem to make a difference
for exe in `cat bigmem.tools`
do
  peflags --cygwin-heap=2048 $exe
done

# Done with the source directories
#rm -rf htslib-1.11 samtools-1.11 bcftools-1.11 bwa bwa-mem2 minimap2

# Add a tmp directory for bash (otherwise it complains)
mkdir win10tools/tmp; chmod 777 win10tools/tmp

# Make the package to download (version label it first?)
gzip win10tools ; rm -rf win10tools


##############################################################################################################
#  Java Runtime libraries of needed Bioinformatics Tools (jartools)
#

# Make the release (package) directory
mkdir jartools; cp $0 jartools; cp unix.tools jartools

# These are just JAR files so simply dump in the directory to be found later
wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.11.zip

# Note: GATK4 is more than just a jar file! Python scripts, etc now
# Todo: process GATK4 for install
# Todo: process IGV for install
tar xjf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 
mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar .
rm -f GenomeAnalysisTK-3.8-1-0-gf15c1c3ef GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

# Make the package to download (version label it first?)
gzip jartools ; rm -rf jartools


#############################################################################################################
#  Python release for WGS Extract (not just portable Python; but all the dependent packages as well)
#   (packages are small enough that easier to load here and deliver instead of loading at each install site)

## Note: Not using the regular Python release but instead the WinPython portable release
# Gives us the full release (unlike the embedded one) with the ability for libraries and such.
# WinPython does not have the latest Python release always though ...

## Python.org embedded release we are not using:
##  wget https://www.python.org/ftp/python/3.7.9/python-3.7.9-embed-amd64.zip; mkdir python-3.7.9; mv python-3.7.9-embed.amd64.zip python-3.7.9; cd python-3.7.9; unzip python-3.7.9-embed-amd64.zip; rm python-3.7.9-embed-amd64.zip

# Note: 3.7.7,1 is last release by WinPython although 3.7.9 is available from Python.org. 32 bit vs 64 bit (?)
wget https://github.com/winpython/winpython/releases/download/2.3.20200530/Winpython32-3.7.7.1dot.exe
./winpython32-3.7.7.1.1dot.exe ; cd WPy32-3771 ; cp -r python-3.7.7 ../ ; cd .. ; rm -rf Wpy32-3771

# Need to update PIP; install packages bio, biopython, biosql, numpy, PIL pillow and pyliftover
cd python-3.7.7
./python,exe -m pip update
./python.exe -m pip install bio numpy PIL pillow pyliftover # biopython biosql 
cd ..

gzip python-3.7.7 ; rm -rf python-3.7.7
