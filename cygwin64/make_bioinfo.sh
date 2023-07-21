#!/bin/bash
# WGS Extract v5
# Script to create Windows release of bioinformatic tools needed for WGS Extract.
# Still have the processing for jartools and similar java toos releases here. Just for documentation.
#
# This content was originally in a make_standalone.sh script in v3
# You need to run cygwin64/cygwin.bat to start a BASH shell in the build environemnt.
# cd /mnt/c/wgse/cygwin64-port
# ./make_bioinfo.sh
#
# Copyright (c) 2020-23 Randy Harr as part of the WGS Extract release
# Released under the xxxxx license

version=5
curdate=$(date +"%d%b%Y")
archive=cygwin64_v${version}_${curdate}_bioinfo.zip
baseURL=https://get.wgse.io/${archive}


if [[ $# -eq 1 && "$1" != "download" ]] ; then
  echo "Usage: $0 [download]"
  echo "  Specify download to download sources (and delete when down)."
  echo "  Else expect them to be already downloadedand available."
  exit 1
fi

curlx="curl -kLZC - --retry 5"

echo "================================================================================================================"
echo "Creating native Windows executables for needed Bioinformatic Tools."
echo

rm -f setup.log || true

# For the bioinformatics tools, we grab the sources, make and install into /usr/local of the cygwin64 base installation.
# We do not install into Windows like a normal release.
#
# This script assumes it is being run in a Cygwin64 BASH shell and the Cygwin development environment is on the PATH.
# To create the environment, see the new companion script make_cygwin64.bat

#----------------------------------------------------------------------------------------------------------------
#  Win10 Implementations of the Bioinformatics Tools (win10tools-bio)
#
htsver=1.17
htsverf=1.17    # Normally a third minor release number like 1.17.1; but not yet created so ...

echo "================================================================================================================"
echo "Building HTSLib"
[ "$1" = "download" ] && \
  (${curlx} https://github.com/samtools/htslib/releases/download/${htsver}/htslib-${htsverf}.tar.bz2 | tar jxf -)
[ -d htslib-${htsverf} ] && \
  (cd htslib-${htsverf}    ;  ./configure --without-libdeflate                ;  make  ;  make prefix=/usr/local install                 ) > setup.log 2>&1

echo
echo "================================================================================================================"
echo "Building Samtools"
[ "$1" = "download" ] && \
  (${curlx} https://github.com/samtools/samtools/releases/download/${htsver}/samtools-${htsverf}.tar.bz | tar zxf -)
[ -d samtools-${htsverf} ] && \
  (cd samtools-${htsverf} || true ;  ./configure --with-htslib=../htslib-${htsverf}  ;  make  ;  make prefix=/usr/local install  ;  make clean  ) >> setup.log 2>&1

echo
echo "================================================================================================================"
echo "Building Bcftools"
[ "$1" = "download" ] && \
  (${curlx} https://github.com/samtools/bcftools/releases/download/${htsver}/bcftools-${htsverf}.tar.bz2 | tar jxf -)
[ -d bcftools-${htsverf} ] && \
  (cd bcftools-${htsverf} || true ;  ./configure --with-htslib=../htslib-${htsverf}  ;  make  ;  make prefix=/usr/local install  ;  make clean  ) >> setup.log 2>&1

[ -d htslib-${htsverf} ] && \
  (cd htslib-${htsverf} || true   ;                                                                                                 make clean  ) >> setup.log 2>&1

echo
echo "================================================================================================================"
echo "Building BWA (Konrads ralloc version)"
[ "$1" = "download" ] && true
  ## BWA is now from our own repository due to a rewrite of the ralloc to get parallel execution
  ## (${curlx} https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar jxf - otherwise)
[ -d bwa-konrad ] && \
  (cd bwa-konrad || true ;              make  ;  cp bwa.exe /usr/local/bin       ;  cp bwa.1 /usr/local/share/man/man1           ;  make clean  ) >> setup.log 2>&1

echo
echo "================================================================================================================"
echo "Building bwa-mem2"
[ "$1" = "download" ] && \
  (${curlx} https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/Source_code_including_submodules.tar.gz | tar zxf -) && \
  vi bwa-mem2-2.2.1/src/utils.h     # see https://github.com/bwa-mem2/bwa-mem2/issues/179 (as of 2022, dup definition)
  #${curlx} https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
  #${curlx} https://github.com/intel/safestringlib/archive/refs/tags/v1.0.0.tar.gz | tar jxf -
  #mv safestringlib-1.0.0 safestringlib; mv -f safestringlib bwa-mem2-2.2.1/ext
  #mkdir bwa-mem2-2.2.1/ext/safestringlib/obj
  # Todo change to sed
[ -d bwa-mem2-2.2.1 ] && \
  (cd bwa-mem2-2.2.1 || true ;              make  ;  cp bwa-mem2* /usr/local/bin                                                 ;  make clean  ) >> setup.log 2>&1

echo
echo "================================================================================================================"
echo "Building minimap2"
[ "$1" = "download" ] && \
  (${curlx} https://github.com/lh3/minimap2/releases/download/v2.25/minimap2-2.25_x64-linux.tar.bz2 | tar jxf -)
[ -d minimap2-2.25 ] && \
  (cd minimap2-2.25  || true ;  make  ;  cp minimap2.exe /usr/local/bin  ;  cp minimap2.1 /usr/local/share/man/man1             ;  make clean  ) >> setup.log 2>&1

echo
echo "================================================================================================================"
echo "Building Fastp"
[ "$1" = "download" ] && \
  (${curlx} https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.4.tar.gz | tar zxf -)
  vi fastp-0.23.4/util.h         # See https://github.com/OpenGene/fastp/issues/318
  # Fastp 0.23 now requies isa-l library from intel
  #  (https://github.com/intel/isa-l) which requires the nasm assembler
  #  (https://www.nasm.us/pub/nasm/releasebuilds/2.15.05/win64/nasm-2.15.05-win64.zip).
  #  See https://github.com/OpenGene/fastp#or-compile-from-source ;
[ -d fastp-0.23.4 ] && \
  (cd fastp-0.23.4  || true ;  make  ;  cp fastp.exe /usr/local/bin                                                            ;  make clean  )  >> setup.log 2>&1

echo
echo "================================================================================================================"
echo "Building Bowtie2"
# [ "$1" = "download" ] && \
#  (${curlx} https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.4.5.tar.gz | tar zxf -)
#[ -d bowtie2-2.4.5 ] && \
#  (cd bowtie2-2.4.5 ; make                                                                                   ;  make clean  )
rm -f bowtie2.zip || true
${curlx} -o bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-mingw-x86_64.zip
unzip -o -qq bowtie2.zip
[ -d bowtie2-2.5.1-mingq-x86_64 ] && (
  cd bowtie2-2.5.1-mingw-x86_64 || true
  mv ./*.exe ./*.bat bowtie2 bowtie2-build bowtie2-inspect /usr/local/bin
  mkdir /usr/local/bowtie2 || true
  cp -r ./* /usr/local/bowtie2
)
rm -rf bowtie2-2.4.5-mingw-x86_64 bowtie2.zip

# Make and Install HiSat2
# ${curlx} https://cloud.biohpc.swmed.edu/index.php/s/fE9QCsX3NH4QwBi/download
# Make and Install PBMM2
# ${curlx} https://github.com/PacificBiosciences/pbmm2/archive/refs/tags/v1.7.0.tar.gz

echo
echo "================================================================================================================"
echo "Cleanup source files; adjust executables"
# Done with the source directories
[ "$1" = "download" ] && \
  rm -rf htslib-${htsverf} samtools-${htsverf} bcftools-${htsverf} bwa-mem2-2.2.1 minimap2-2.24 fastp-0.23.2 bowtie2-2.4.5-mingw-x86_64

# Set loader flag to allow big memory model for certain executables (goes from default 512MB to 2GB max)
# Note: was hoping this would help with low memory utilization in samtools, etc but does not seem to make a difference
while IFS= read -r exe ; do
  peflags --cygwin-heap=2048 "/usr/local/bin/$exe"
done < bigmem.tools

# DLLs are not being found in /bin on path from /usr/local binaries; so dup them in /usr/local
files="$(echo /usr/local/bin/*exe)"
rm -f dep.tmp
for exe in $files; do
  cygcheck "$exe" | grep -v "Randy\|WINDOWS\|jre\|exe" >>dep.tmp
done
while read -r dll ; do
  cp "$dll" /usr/local/bin
done < <(sort dep.tmp | uniq)
rm dep.tmp

echo
echo "================================================================================================================"
echo "Make the Bioinformatic Tools ZIP package"
# Dropped GoogleDrive as have to encrypt. Google Drive does not like exe's in the zip
# 7z a -tzip win10tools.zip -pWGSEv3 win10tools >/dev/null
cp -r 00README-bioinfo.txt make_bioinfo.sh open_source_licenses bigmem.tools /usr/local
make_json() {   # Parameter is top-level name of json file
printf "{ \n"
printf "  \"%s\": { \n"  "$1"
printf "    \"version\": %s, \n"   "${version}"
printf "    \"date\": \"%s\", \n"  "${curdate}"
printf "    \"URL\": \"%s\" \n"    "${baseURL}"
printf " } \n"
printf "} \n"
}
make_json bioinfo > /usr/local/bioinfo.json
#echo "{\"bioinfo.version\": ${version}, \"bioinfo.date\": \"${curdate}\", \"bioinfo.URL\": \"${baseURL}\"}" > /usr/local/version.json
rm -f ${archive} || true
7z a -tzip ${archive} /usr/local >/dev/null
echo "Finished making Bioinformatics tools."
echo "To retrieve package during install, use 7z x -o/usr ${archive}"
echo "Review the file setup.log for warnings and further details of each compile"
echo

exit

#---------------------------------------------------------------------------------------------------------------------
#  Java Runtime libraries of needed Bioinformatics Tools (jartools)
# NOTE: Only Haplogrep used at this time.  FastQC and IGV are whole directory structures and installed directly.

# Make the release (package) directory
mkdir jartools
cp "$0" jartools

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
(
  cd fastqc || true
  # FastQC did not do the final step to create a single JAR file with main-entry point manifest specified; so we do it now
  echo -e "Class-Path: cisd-jhdf5.jar jbzip2-0.9.jar sam-1.103.jar\nMain-Class: uk.ac.babraham.FastQC.FastQCApplication" >Manifest-Update.txt
  jar cfm FastQC.jar Manifest-Update.txt ./*.jar uk org net LICENS* README* Configuration Template Help ./*ico
  rm -rf cisd-jhdf5.jar jbzip2-0.9.jar sam-1.103.jar uk org net
  mv FastQC.jar ../jartools
)
mv fastqc jartools # needs Config, Templates, etc left in there?
# Multiqc is a python library package

# Haplogrep Install
gunzip haplogrep.zip
mv haplogrep.jar jartools
rm haplogrep haplogrep.zip

# Make the package to download (version label it first?)
zip -qr jartools jartools
rm -rf jartools
