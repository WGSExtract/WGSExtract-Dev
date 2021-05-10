#!/bin/bash -x
#
echo WGS Extract Beta v3 Win10 Bioinformatic and Linux Support ****WORKING DRAFT -- NOT EXECUTABLE YET****
exit
#
# Script to create Win10 release of tools needed for WGS Extract install.  Since no package manager is available for
#  Windows, we need to create pseudo-packages here that we then make available to the installers from our own
#  repository.  The bioinformatic tools are handled in one large package.  And the Python release in another.
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
#  For Python, we use the more portable WinPython release.  We only need the actual python directory though. And it comes
#  as a zip archive executable.  So we do some minimal processing to trim to only what we need. Win10, as MacOSX, has
#  really messed with the Python environment by putting pseudo executables in that then simply redirect to their software
#  store.  To avoid all the issues of that and getting the tested, correct version, we simply install our own local copy
#  of Python. We simply bring in the necessary packages now as well (PIP libraries; similar to Java Jar files).
#
# This script assumes it is being run in a Cygwin64 BASH shell and the Cygwin development environment is on the exe path
#  The cygwin64 development environment consists of the following main packages (and their dependencies) beyond the basic
#  release:
#     TBD
#
# Copyright (c) 2020-21 Randy Harr
# Released under the xxxxx license
#

#---------------------------------------------------------------------------------------------------------------------
#  Win10 Cygwin64 release of minimal Unix environment to support the bioinformatic tools and WGS Extract scripts
#     (win10tools-cygwin64)

# Make our release directory. Copy these scipts simply for documentation to the end user in case they want to remake
mkdir -p win10tools/bin
cp make_standalone.sh unix.tools bigmem.tools win10tools

# Grab needed and useful Unix utilities from the main CygWin bin
for exe in $(cat unix.tools); do
  [ -f /bin/$exe ] && cp /bin/$exe win10tools/bin
done

# 7z in /bin is a shell script to run programs out of /usr/lib
mkdir -p win10tools/lib
cp -r /usr/lib/p7zip win10tools/lib
cp /usr/lib/p7zip/* win10tools/bin                   # Let's simply copy to bin as well; as non .exe files are not found

# Special; alias that we cannot store in ZIP file. But needed in Universal Update script
cp -T win10tools/bin/gawk.exe win10tools/bin/awk.exe # simplifies with other OS's that call it awk

# Cygwin does not have sudo; so created our own batch shell script to emulate
cp sudo.bat win10tools/bin

# Need certs for Curl; could never figure out how to get them working. So dropped Curl as Win10 has it

# Check all the *.exe's in bin for dll dependencies; copy all those to the bin as well
files="$(echo win10tools/bin/*exe)"
for exe in $files; do
  cygcheck "$exe" | grep -v Randy | grep -v WINDOWS | grep -v jdk | grep -v exe >dep.tmp
  for dll in $(cat dep.tmp); do
    cp $dll win10tools/bin
  done
done
rm dep.tmp
cp /usr/lib/libbash.dll.a win10tools/lib

mkdir win10tools/tmp  ;  sudo.bat chown "root:root" "win10tools/tmp"  ;  sudo.bat chmod 1777 win10tools/tmp

zip -qr win10tools-cygwin64 win10tools
rm -rf win10tools

#---------------------------------------------------------------------------------------------------------------------
#  Win10 Implementations of the Bioinformatics Tools (win10tools-bio)
mkdir -p win10tools/bin
mkdir -p win10tools/share/man/man1

# Grab the latest releases that we tested with
# curl -L https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 | tar jxf
# curl -L https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz | tar jxf
# curl -L https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 | tar jxf
## BWA is now from our own repository due to a rewrite of the ralloc to get parallel execution
## curl -L https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar jxf
# curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 | tar jxf -
# curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar jxf -
# curl -L https://github.com/OpenGenes/fastp/releases/download/

cd htslib-1.12    ;  configure  ;  make  ;  make prefix=../win10tools install  ;  make clean  ;  cd ..
cd samtools-1.12  ;  configure  ;  make  ;  make prefix=../win10tools install  ;  make clean  ;  cd ..
cd bcftools-1.12  ;  configure  ;  make  ;  make prefix=../win10tools install  ;  make clean  ;  cd ..
cd bwa-konrad     ;  make  ;  cp bwa.exe ../win10tools/bin  ;  cp bwa.1 ../win10tools/share/man/man1  ;  make clean  ;  cd ..
cd bwa-mem2       ;  make  ;  cp bwa-mem2* ../win10tools/bin  ;  make clean  ;  cd .. # no man page to copy
cd minimap2       ;  make  ;  cp minimap2.exe ../win10tools/bin  ;  cp minimap2.1 ../win10tools/share/man/man1  ;  make clean  ;  cd ..
cd fastp          ;  make  ;  cp fastp.exe ../win10tools/bin  ;  make clean  ;  cd ..

# Done with the source directories
#rm -rf htslib-1.12 samtools-1.12 bcftools-1.12 bwa bwa-mem2 minimap2 fastp

# Check all the *.exe's in bin for dll dependencies; copy all those to the bin as well
files="$(echo win10tools/bin/*exe win10tools/bin/bwa-mem2*)"
for exe in $files; do
  cygcheck $exe | grep -v Randy | grep -v WINDOWS | grep -v jdk | grep -v exe >dep.tmp
  for dll in $(cat dep.tmp); do
    cp $dll win10tools/bin
  done
done
rm dep.tmp

# Set loader flag to allow big memory model for certain executables (goes from default 512MB to 2GB max)
# Note: was hoping this would help with low memory utilization in samtools, etc but does not seem to make a difference
for exe in $(cat bigmem.tools); do
  peflags --cygwin-heap=2048 win10tools/bin/$exe
done

# Add a tmp directory for bash; shm for bwa (otherwise they complain)
mkdir win10tools/dev        ;  chmod 770 win10tools/dev
mkdir win10tools/dev/shm    ;  sudo.bat chmod 01777 win10tools/dev/shm
mkdir win10tools/dev/mqueue ;  sudo.bat chmod 01777 win10tools/dev/mqueue
#cp /dev/std* win10tools/dev
#cp /dev/fd win10tools/dev

# Make the package to download (version label it first?)
# Note: have to encrypt now due to Google Drive not liking the exe's in the zip
7z a -tzip win10tools.zip -pWGSEv3 win10tools >/dev/null
mv win10tools.zip win10tools-bioinfo1.12.zip
rm -rf win10tools

#---------------------------------------------------------------------------------------------------------------------
#  Java Runtime libraries of needed Bioinformatics Tools (jartools)
#

# Make the release (package) directory
mkdir jartools
cp $0 jartools

# These are just JAR files so simply dump in the directory to be found later
wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.11.zip
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
wget https://github.com/seppinho/haplogrep-cmd/releases/download/v2.2.9/haplogrep.zip

# Gatk4 install
gunzip gatk-4.1.9.0.zip
mv gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar jartools/GATK4.jar
rm -rf gatk-4.1.9.0.zip gatk-4.1.9.0

# Gatk3 install
tar xjf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar jartools/GATK3.jar
rm -rf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

# Todo: process picard for install
# Todo: process IGV for install

# FastQC Install
gunzip fastqc_v0.11.9.zip
cd fastqc
# FastQC did not do the final step to create a single JAR file with main-entry point manifest specified; so we do it now
echo -e "Class-Path: cisd-jhdf5.jar jbzip2-0.9.jar sam-1.103.jar\nMain-Class: uk.ac.babraham.FastQC.FastQCApplication" >Manifest-Update.txt
jar cfm FastQC.jar Manifest-Update.txt *.jar uk org net LICENS* README* Configuration Template Help *ico
rm -rf cisd-jhdf5.jar jbzip2-0.9.jar sam-1.103.jar uk org net
mv FastQC.jar ../jartools
cd ..
mv fastqc jartools # needs Config, Templates, etc left in there

# Haplogrep Install
gunzip haplogrep.zip
mv haplogrep.jar jartools
rm haplogrep haplogrep.zip

# Make the package to download (version label it first?)
zip -qr jartools jartools
rm -rf jartools

#---------------------------------------------------------------------------------------------------------------------
#  Python release for WGS Extract (not just portable Python; but all the dependent packages as well)
#   (packages are small enough that easier to load here and deliver instead of loading at each install site)

## Note: Not using the regular Python release but instead the WinPython portable release
# Gives us the full release (unlike their embedded one) with the ability to add libraries and such.
# WinPython does not have the latest Python release though ...

## Python.org embedded release we are not using:
##  wget https://www.python.org/ftp/python/3.7.9/python-3.7.9-embed-amd64.zip;
##  mkdir python-3.7.9;
##  mv python-3.7.9-embed.amd64.zip python-3.7.9;
##  cd python-3.7.9;
##  unzip python-3.7.9-embed-amd64.zip;
##  rm python-3.7.9-embed-amd64.zip

# Note: 3.7.7,1 is last release by WinPython although 3.7.9 is available from Python.org. Numpy requires 64bit.
#wget https://github.com/winpython/winpython/releases/download/2.3.20200530/Winpython64-3.7.7.1dot.exe
wget https://github.com/winpython/winpython/releases/download/4.1.20210417/Winpython64-3.8.9.0dot.exe
./winpython64-3.8.9.0dot.exe
cd WPy64-3890  ;  mv python-3.8.9 ..  ;  cd ..  ;  rm -rf Wpy64-3890

# Need to update PIP; install packages pandas pillow pyliftover pyscreenshot (pandas needed by yleaf)
( cd python-3.8.9
  ./python,exe -m pip install --update pip
  ./python.exe -m pip install pillow pyliftover pyscreenshot pandas --no-warn-script-location
   # bio (brings biopython, numpy, requests urllib3), pyfaidx, pyfaidx, pysocks
)

mv python-3.8.9 python
zip -qr python python
rm -rf python
