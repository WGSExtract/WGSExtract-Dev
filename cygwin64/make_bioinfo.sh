#!/bin/bash -x
#
# Script to create Windows release of bioinformatic tools needed for WGS Extract.
# Still have the processing for jartools and similar java toos releases here. Just for documentation.
#
# This content was originally in a make_standalone.sh script in v3
#
# Copyright (c) 2020-22 Randy Harr as part of the WGS Extract release
# Released under the xxxxx license

if [[ $# -eq 1 && "$1" = "download" ]] ; then
  download="True"
else
  download="False"
fi

# For the bioinformatics tools, we grab the sources, make and install into /usr/local of the cygwin64 base installation.
# We do not install into Windows like a normal release.
#
# This script assumes it is being run in a Cygwin64 BASH shell and the Cygwin development environment is on the PATH.
# To create the environment, see the new companion script make_cygwin64.bat

curlx="curl -kLZC - --retry 5"

#----------------------------------------------------------------------------------------------------------------
#  Win10 Implementations of the Bioinformatics Tools (win10tools-bio)
#
BASE=cygwin64-bioinfo
VERS=v2
ARCH=${BASE}_${VERS}.zip
SAMVERS=1.15.1

[ "$download" = "True" ] && \
  (${curlx} https://github.com/samtools/htslib/releases/download/1.15/htslib-${SAMVERS}.tar.bz2 | tar jxf -)
[ -d htslib-${SAMVERS} ] && \
  (cd htslib-${SAMVERS}    ;  ./configure                                 ;  make  ;  make prefix=/usr/local install                 ;  cd ..)

[ "$download" = "True" ] && \
  (${curlx} https://github.com/samtools/samtools/releases/download/1.15/samtools-${SAMVERS}.tar.bz | tar zxf -)
[ -d samtools-${SAMVERS} ] && \
  (cd samtools-${SAMVERS}  ;  ./configure --with-htslib=../htslib-${SAMVERS}  ;  make  ;  make prefix=/usr/local install  ;  make clean  ;  cd ..)

[ "$download" = "True" ] && \
  (${curlx} https://github.com/samtools/bcftools/releases/download/1.15/bcftools-${SAMVERS}.tar.bz2 | tar jxf -)
[ -d bcftools-${SAMVERS} ] && \
  (cd bcftools-${SAMVERS}  ;  ./configure --with-htslib=../htslib-${SAMVERS}  ;  make  ;  make prefix=/usr/local install  ;  make clean  ;  cd ..)

[ -d htslib-${SAMVERS} ] && \
  (cd htslib-${SAMVERS}    ;                                                                                             make clean  ;  cd ..)

[ "$download" = "True" ] && true
  ## BWA is now from our own repository due to a rewrite of the ralloc to get parallel execution
  ## (${curlx} https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar jxf -)
[ -d bwa-konrad ] && \
  (cd bwa-konrad     ;  make  ;  cp bwa.exe /usr/local/bin       ;  cp bwa.1 /usr/local/share/man/man1       ;  make clean  ;  cd ..)

[ "$download" = "True" ] && \
  (${curlx} https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/Source_code_including_submodules.tar.gz | tar zxf -) && \
  vi bwa-mem2-2.2.1/src/utils.h     # see https://github.com/bwa-mem2/bwa-mem2/issues/179 (as of 2022, dup definition)
  #${curlx} https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
  #${curlx} https://github.com/intel/safestringlib/archive/refs/tags/v1.0.0.tar.gz | tar jxf -
  #mv safestringlib-1.0.0 safestringlib; mv -f safestringlib bwa-mem2-2.2.1/ext
  #mkdir bwa-mem2-2.2.1/ext/safestringlib/obj
  # Todo change to sed
[ -d bwa-mem2-2.2.1 ] && \
  (cd bwa-mem2-2.2.1 ;  make  ;  cp bwa-mem2* /usr/local/bin                                                 ;  make clean  ;  cd ..) # no man page to copy

[ "$download" = "True" ] && \
  (${curlx} https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar jxf -)
[ -d minimap2-2.24 ] && \
  (cd minimap2-2.24  ;  make  ;  cp minimap2.exe /usr/local/bin  ;  cp minimap2.1 /usr/local/share/man/man1  ;  make clean  ;  cd ..)

[ "$download" = "True" ] && \
  (${curlx} https://github.com/OpenGene/fastp/archive/refs/tags/v0.21.0.tar.gz) && \
  vi fastp-0.21.0/util.h         # See https://github.com/OpenGene/fastp/issues/318
  # Fastp 0.23 now requies isa-l library from intel
  #  (https://github.com/intel/isa-l) which requires the nasm assembler
  #  (https://www.nasm.us/pub/nasm/releasebuilds/2.15.05/win64/nasm-2.15.05-win64.zip).
  #  See https://github.com/OpenGene/fastp#or-compile-from-source ;
  # ${curlx} https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.2.tar.gz
[ -d fastp-0.21.0 ] && \
  (cd fastp-0.21.0   ;  make  ;  cp fastp.exe /usr/local/bin                                                 ;  make clean  ;  cd ..)

# ToDo HiSat2, PBMM2, bowtie2, bedtools, centrifuge ...
# Make and Install HiSat2
# ${curlx} https://cloud.biohpc.swmed.edu/index.php/s/fE9QCsX3NH4QwBi/download
# Make and Install PBMM2
# ${curlx} https://github.com/PacificBiosciences/pbmm2/archive/refs/tags/v1.7.0.tar.gz

# Done with the source directories
[ "$download" = "True" ] && \
  rm -rf htslib-${SAMVERS} samtools-${SAMVERS} bcftools-${SAMVERS} bwa-mem2-2.2.1 minimap2-2.24 fastp-0.23.2

# Set loader flag to allow big memory model for certain executables (goes from default 512MB to 2GB max)
# Note: was hoping this would help with low memory utilization in samtools, etc but does not seem to make a difference
for exe in $(cat bigmem.tools); do
  peflags --cygwin-heap=2048 /usr/local/bin/$exe
done

# DLLs are not being found in /bin on path from /usr/local binaries; so dup them in /usr/local
files="$(echo /usr/local/bin/*exe)"
rm -f dep.tmp
for exe in $files; do
  cygcheck "$exe" | grep -v Randy | grep -v WINDOWS | grep -v jre | grep -v exe >>dep.tmp
done
sort dep.tmp | uniq > depsu.tmp
for dll in $(cat depsu.tmp); do
  cp $dll /usr/local/bin
done
rm dep.tmp depsu.tmp

# Make the package to download (version label it first?)
# Dropped GoogleDrive as have to encrypt. Google Drive does not like exe's in the zip
# 7z a -tzip win10tools.zip -pWGSEv3 win10tools >/dev/null
cp -r 00README-bioinfo.txt make_bioinfo.sh open_source_licenses bigmem.tools /usr/local
echo ${SAMVERS} > /usr/local/version.txt
rm -f ${ARCH} || true
7z a -tzip ${ARCH} /usr/local >/dev/null
# To retrieve package during install, use 7z x -o/usr ${ARCH} ; rm ${ARCH}

exit

#---------------------------------------------------------------------------------------------------------------------
#  Java Runtime libraries of needed Bioinformatics Tools (jartools)
# NOTE: Only Haplogrep used at this time.  FastQC and IGV are whole directory structures and installed directly.

# Make the release (package) directory
mkdir jartools
cp $0 jartools

# These are just JAR files so simply dump in the directory to be found later
wget https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_2.8.11.zip
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
wget https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip

# Gatk4 install
gunzip gatk-4.1.9.0.zip
mv gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar jartools/GATK4.jar
rm -rf gatk-4.1.9.0.zip gatk-4.1.9.0

# Gatk3 install
tar xjf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar jartools/GATK3.jar
rm -rf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

# Picard install
#  Nothing to do; simply a .jar

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
mv fastqc jartools # needs Config, Templates, etc left in there?
# Multiqc is a python library package

# Haplogrep Install
gunzip haplogrep.zip
mv haplogrep.jar jartools
rm haplogrep haplogrep.zip

# Make the package to download (version label it first?)
zip -qr jartools jartools
rm -rf jartools
