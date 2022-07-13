#!/bin/bash
# WGS Extract v2 to v3 in-place Upgrade Script
# Copyright (C) 2021 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
#

# Do any v3 and OS specific environment installs as part of the upgrade
# Grab the main release if only install scripts delivered so far.
# Try to reuse v2 Blobs to avoid downloading again; mostly the reference_genomes directory
# Will delete old v2 release stuff once finished grabbing what we can (if it exists)
# Will fill in any missing reference genomes

# In case called directly, get path to wintools setup first.  The Win10 Update script calls the Install_Win10 which
#  does the initial win10tools grab of the cygwin environment (small enough for Win10 to grab)
case $OSTYPE in
  darwin*|linux*)
    ;;
  msys*|cygwin*)
    if [ ! -d win10tools/ ] ; then
      # We are hosed as somehow this script has been called directly and not from the installation directory
      # We have no access to echo, printf or any error reporting function; but we try anyway just in case
      echo "*** ERROR: Cannot find the Windows 10 Cygwin64 tools directory previously installed."
      exit 1
    fi
    ;;
  *)
    echo "*** Error: unknown OSTYPE of $OSTYPE"
    exit 1
    ;;
esac


# Special function to get large google drive files that have "cannot virus scan" pop-up (>50 MB)
# Inspiration from https://stackoverflow.com/questions/25010369/wget-curl-large-file-from-google-drive
# which has python code version also; for when we switch to that for reference genome management
#  First parameter is Google Drive File ID; Second optional is (final) name of file to download
#  having the 2nd parameter helps when the filename changes or is unknown
get_googledrive_large_file () {
  ggID=$1
  ggURL='https://drive.google.com/uc?export=download'
  filename="$(curl -sc cookie.xml "${ggURL}&id=${ggID}" | grep -o '="uc-name.*</span>' | sed 's/.*">//;s/<.a> .*//')"
  # To get confirm code from initial response; tee above and use below.  Else, can get from cookie as done below
  # code=$(awk '/confirm=[A-Za-z0-9]{4}/ { split($2,arr,"=") ; print arr[2] }' RS="<" FS="&" cookie.html) ;
  code=$(awk '/_warning_/ {print $NF}' cookie.xml)
  code=${code%$'\r'}        # For Win10; cygwin gawk awk adds trailing \r for some reason
  curl -Lkb cookie.xml "${ggURL}&confirm=${code}&id=${ggID}" -o "${filename}"
  if [[ ($# == 2 && ${filename} != "$2") ]]; then
    mv -f "${filename}" "$2"
  fi
  rm cookie.xml
}


# Need to find the installation as CWD could be anywhere if this script called from within GUI (in MacOS)
WGSEDIR=$(dirname "$0")         # Get the directory path of this script
WGSEABS=$(cd "$WGSEDIR" || echo "***Internal ERROR: cd $WGSEABS" ; pwd)   # Extra to get through any symlinks, aliases, etc
WGSEESC=${WGSEABS/ /\\}         # Escape embedded spaces
cd "$WGSEESC" || echo "***Internal ERROR: cd $WGSEESC"


# Do OS Specific Upgrades / additions for v3 from v2
case "$OSTYPE" in
  darwin*)
    echo "*** OSX Specific Upgrades"
    LOCALPORT=/opt/local/bin/port
    if [ ! -f /opt/local/bin/7z ]; then
      sudo -H $LOCALPORT -N install p7zip md5sha1sum # Curl, bunzip2 is already in MacOS basic package
    fi
    sudo -H $LOCALPORT upgrade outdated

    # BWA is not in MacPorts! So grab a compiled version from github user smikkelsendk and install in macports bin
    if [ ! -f /opt/local/bin/bwa ]; then
      arch=$(uname -m)
      case "$arch" in
        x86_64*)        # For old Intel processors historically out there
          curl -Lo bwa https://raw.githubusercontent.com/smikkelsendk/bwa-for-arm/master/bin/bwa0717-mac-x64
          sudo mv -f bwa /opt/local/bin/bwa ;;
        arm*)           # For new M1 processor we start to encounter
          curl -Lo bwa https://raw.githubusercontent.com/smikkelsendk/bwa-for-arm/master/bin/bwa0717-mac-arm64
          sudo mv -f bwa /opt/local/bin/bwa ;;
        *)  echo "*** Error: Unknown MacOS architecture $arch" ;;
      esac
    fi

    PATH="/opt/local/bin:$PATH"
    ;;

  linux*)
    echo "*** LINUX Specific Upgrades"
    # todo check if already installed before requesting reinstall
    sudo apt-get install bwa curl p7zip-full     # Need BWA, Curl bunzip2 and 7zip in v3; bunzip2 already in basic
    ;;

  msys*|cygwin*)
    echo "*** WINDOWS Specific Upgrades / Installs"

    # We check to see if Java command is already available; or if we have locally installed it already
    # Not part of v2 release nor install check; so added to update here
    if ! (command -v java &> /dev/null || [ -f jre/bin/java.exe ]) ; then
      echo Installing Java JRE
      # We do not need standalone Java release but easier to setup that way in batch mode
      openjdk="https://github.com/AdoptOpenJDK/openjdk11-binaries/releases/download/jdk-11.0.10%2B9/"
      curl -Lk -o jre.zip  $openjdk/OpenJDK11U-jre_x64_windows_hotspot_11.0.10_9.zip
      # unzip -qbo jre.zip  &&  rm -f jre.zip
      # 7z x -tzip -y jre.zip > /dev/null  &&  rm -f jre.zip
      powershell Expand-Archive -LiteralPath "jre.zip" -DestinationPath "." -Force  &&  rm -f jre.zip
      mv jdk-11.0.10+9-jre jre
      # Program does two similar checks; either simply available or if Win10, in the jre subdirectory
      echo finished installing Java JRE
    fi

    # WinPython standalone Python release that we pre-setup with specific version and packages (see make_standalone.sh)
    # Was simply redistributed in v2 to everyone no matter what the OS; now only installed with Win10 users
    if [ ! -d python ]; then
      echo Installing Python 3.7.7
      # get_googledrive_large_file "1dfaxQzF_ugBsO7dnZXlkRjSdEd36nuVT" "python.zip"
      curl -L -o python.zip 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjcrtmsj2NG-x08fk/root/content'
      # unzip -qbo python.zip  &&  rm -f python.zip
      # 7z x -tzip -y python.zip > /dev/null  &&  rm -f python.zip
      powershell Expand-Archive -LiteralPath "python.zip" -DestinationPath "." -Force  &&  rm -f python.zip
      echo finished installing Python 3.7.7
      # Already includes PIP library setup; see make_standalone.sh
    fi

    # Grab our own bioinformatics v1.12 release to add to the existing cygwin64 Unix tools already installed
    # Was simply redestributed in v2 to everyone no matter what the OS; now only installed with Win10 users
    if [ ! -f win10tools/bin/samtools.exe ]; then   # Bioinformatics tools (made with CygWin64; merge with Cygwin64)
      echo Installing Bioinformatics tools -- cygwin64 port for Win10
      # Actual file name is win10tools-bioinfo1.12.zip with encrypted files
      # get_googledrive_large_file "1qlBWByCnTKTyeJa3xzeI2qgnyAEDK5BR" "win10tools.zip"
      curl -L -o win10tools.zip 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjdiI35GxxlmrBLU8/root/content'
      # Google decided this zip file, of the many we have, was worthy of 4 level "are your sure" check. So we
      #  encrypted it to bypass their virus check.  Need 7z to unencrypt via command line; no others allow password
      7z x -tzip -y -pWGSEv3 win10tools.zip > /dev/null
      rm -f win10tools.zip
      echo finished installing Bioinformatoc tools.
    fi
    # Sometimes fork error in rm above with Win10 Cygwin64 so try again to make sure
    [ -f win10tools.zip ]   && rm -f win10tools.zip

    # Final cleanup that is not handled by initial Cygin64 install via cmd.exe
    chmod 777 win10tools/tmp
    ;;
  *)
    echo "*** Error: unknown OSTYPE of $OSTYPE" ;;      # Don't have OS specific echo ....
esac


# Get the full WGS Extract v3 release now; just got the Install / Update script(s) during user download
# This is common for all users and includes python code, haplogrep in jartools, yleaf and 400 MB in templates
# Also includes new, basic reference (library) with genome model processing (but no genomes)
if [ ! -d program ]; then
  echo Installing WGS Extract v3 program and its templates ...
  # Actual file: WGSExtractAlphav3j_14Apr2021_Full.zip
  #get_googledrive_large_file "1fbgKbXiVB6ELkRPHV5LAY8SsYJgfJlr5" "WGSExtractv3.zip"
  curl -L -o WGSExtractv3.zip 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjc1OE-wG_vQ1kzdc/root/content'
  case $OSTYPE in
    msys*|cygwin*)
      powershell Expand-Archive -LiteralPath "WGSExtractv3.zip" -DestinationPath "." -Force  &&  rm -f WGSExtractv3.zip
      ;;
    darwin*|linux*)
      # unzip -qbo WGSExtractv3.zip  &&  rm -f WGSExtractv3.zip
      7z x -tzip -y WGSExtractv3.zip > /dev/null  &&  rm -f WGSExtractv3.zip
      ;;
  esac
  [ ! -d open_source_licenses ] && mkdir open_source_licenses
  mv -f WGSExtractv3/open_source_licenses/* open_source_licenses && rmdir WGSextractv3/open_source_licenses
  mv -f WGSExtractv3/* .  &&  rmdir WGSExtractv3
  echo finished installing WGS Extract v3 program and templates.
fi
if [ -d WGSExtractv3 ]; then
  # Win10 Cygwin64 is barfing on a fork call in mv above; so cleanup next time through -- temporary path for now
  mv -f WGSExtractv3/open_source_licenses/* open_source_licenses && rmdir WGSextractv3/open_source_licenses
  mv -f WGSExtractv3/* .  &&  rmdir WGSExtractv3
  echo finished installing WGS Extract v3 program and templates.
fi


#
# Ready to perform actual v2 to v3 transformation of any Blobs we want to transfer, recreate or download
#

# Handle removing old start / install scripts (2b release and patches)

# Avoiding compiled Applescript due to Translocation issues and Apple not allowing distribution outside signed app
[ -d Install_MacOSX.app ]    && rm -rf Install_MacOSX.app
[ -d Start_MacOSX.app ]      && rm -rf Start_MacOSX.app
[ -d Uninstall_MacOSX.app ]  && rm -rf Uninstall_MacOSX.app

# Removed all Apple script methods to start shell files; changed .sh to .command for easy start
[ -f Install_MacOSX.scpt ]   && rm -f Install_MacOSX.scpt
[ -f Start_MacOSX.scpt ]     && rm -f Start_MacOSX.scpt
[ -f Uninstall_MacOSX.scpt ] && rm -f Uninstall_MacOSX.scpt

# renamed to just MacOS due to BigSur 11; renamed .sh to .command for easier click-start
[ -f Install_MacOSX.sh ]     && rm -f Install_MacOSX.sh
[ -f Start_MacOSX.sh ]       && rm -f Start_MacOSX.sh
[ -f Uninstall_MacOSX.sh ]   && rm -f Uninstall_MacOSX.sh

# Changed all OS specific files to Start_xxxx; also renamed MacOS from .sh to .command
[ -f Windows_START.bat ]     && rm -f Windows_START.bat
[ -f MacOS_START.sh ]        && rm -f MacOS_START.sh
[ -f Linux_START.sh ]        && rm -f Linux_START.sh

#
# We simply bundle most of the templates and similar blobs with the new release.
# But if we wanted to transform instead of re-downloading, here is what we would do.
#

# If we wanted to reuse the microarray templates; small enough that simply included with main release
# mkdir program; mkdir program/microarray;
# mv programs/Extract23/raw_file_templates program/microarray
# mv programs/Extract23/All_SNPs* program/microarray
# cd program/microarray
# The GRCh38 version was missing from WGSE v2; create now
# zcat All_SNPs_hg38_ref.tab.gz | grep s/^chr// > ALL_SNPs_GRCh38_ref.tab
# bgzip -i ALL_SNPs_GRCh38_ref.tab
# Need Ploidy override file now; CombinedKit generator assumes all chromosomes are ploidy 2
# echo "* * * F 2\n* * * M 2" > ploidy.txt
# cd ../..

# If we wanted to reuse the yleaf tables; small enough that simply included with main release
# mkdir yleaf
# mv programs/yleaf/Position_files yleaf
# mv programs/yleaf/Hg_Prediction_tables yleaf

# If we wanted to reuse the haplogrep JAR file; small enough that simply included in main release
# mkdir jartools
# mv programs/haplogrep/haplogrep.jar jartools

# Finished with saving files from the old programs/ folder; can now remove directory structure of v2b programs/
[ -d programs ] && rm -rf programs      # Removes old Win10 binary installs as well

#
# In v2 to v3, change from reference_genomes folder to reference/genomes. Reference/genomes should exist already and
# contain the process_reference_genomes.sh file. That is key to creating the index files.
#

# Move v2 Reference Genomes to new v3 area
if [ -d reference_genomes ] && [ -d reference/genomes ]; then
  echo "Saving existing reference genomes from v2 release" ;
  cd reference_genomes || echo "***Internal ERROR: cd reference_genomes" ;
  [ -f hg38.fa.gz ] && echo Moving hg38 && mv -f hg38.fa.gz ../reference/genomes/hg38.fa.gz ;
  [ -f hs37d5.fa.gz ] && echo Moving hs37d5 && mv -f hs37d5.fa.gz ../reference/genomes/hs37d5.fa.gz ;
  [ -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ] && echo Moving hs38 && mv -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ../reference/genomes/hs38.fa.gz ;
  [ -f human_g1k_v37.fasta.gz ] && echo Moving human_g1k_v37 && mv -f human_g1k_v37.fasta.gz ../reference/genomes/human_g1k_v37.fasta.gz ;
  [ -f hg19.fa.gz ] && Moving hg19_wgse && mv -f hg19.fa.gz ../reference/genomes/hg19_wgse.fa.gz ;
  cd ..
  echo "Saved needed reference files. Removing v2 release reference genomes directory" ;
  rm -rf reference_genomes
fi

# Note, a new v3 should simply overlay an existing v3. So a v3 to v3 upgrade is automatic. Even
#  will check if a file exists before trying to create it.

# Download any missing reference genomes; if desired to do now
if [ -d reference/genomes ]; then
  echo ""
  echo "Do you wish to download and process missing reference genomes now?"
  echo " (Could take an hour or more. 4 more of 9 in total. Each is around 1 GB.)"
  select yn in "Yes" "No"; do
    case $yn in
      No )
        echo "OK. You can always rerun the Upgrade / Installer later to download the missing reference genomes."
        echo ""
        exit ;;
      Yes )
        break ;;
    esac
  done
  cd reference/genomes || echo "***Internal ERROR: cd reference/genomes" ;

  echo "Downloading missing reference genomes that were not saved from a previous v2 release."
  [ ! -f hg38.fa.gz ] && echo Grabbing hg38 && curl -o hg38.fa.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz ;
  [ ! -f hs37d5.fa.gz ] && echo Grabbing hs37d5 && curl -k -o hs37d5.fa.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz ;
  [ ! -f hs38.fa.gz ] && echo Grabbing hs38 && curl -o hs38.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ;
  [ ! -f human_g1k_v37.fasta.gz ] && echo Grabbing human_g1k_v37 && curl -o human_g1k_v37.fasta.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;
  # [ ! -f hg19_wgse.fa.gz ] && get_googledrive_large_file "1O-Kb8dpG3wccfbPeClBg0oM33eT7UwYd" "hg19_wgse.fa.gz"
  [ ! -f hg19_wgse.fa.gz ] && echo Grabbing hg19_wgse && curl -L -o hg19_wgse.fa.gz 'https://api.onedrive.com/v1.0/shares/s!AgorjTSMFYpjde1tXQ1bDvCrlFo/root/content'
  echo "Downloading new reference genomes for v3 release that did not exist in v2."
  [ ! -f hg19.fa.gz ] && echo Grabbing hg19 && curl -o hg19.fa.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz ;
  [[ ! -f hg19_yseq.fa.gz && ! -f hg19_yseq.fa.zip ]] && echo Grabbing hg19_yseq && curl -o hg19_yseq.fa.zip http://genomes.yseq.net/WGS/ref/hg19/hg19.zip ;
  [ ! -f Homo_sapiens.GRCh37.dna.toplevel.fa.gz ] && echo Grabbing Homo_sapiens.GRCh37.dna.toplevel && curl -o Homo_sapiens.GRCh37.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz ;
  [ ! -f Homo_sapiens.GRCh38.dna.toplevel.fa.gz ] && echo Grabbing Homo_sapiens.GRCh38.dna.toplevel && curl -o Homo_sapiens.GRCh38.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz ;

  # Process reference genomes for various indices and extracted parameters. Will take an hour more.
  # Critical that the Unix / Linux tools be available on the path for this to work.
  bash process_reference_genomes.sh .  # Create standard indices and support files for each reference model
  cd ../..
fi
