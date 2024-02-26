#!/usr/bin/env bash
# WGS Extract MacOSX Install Script
# Copyright (C) 2018-2020 Marko Bauer
# Copyright (C) 2020-2022 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
#

# Todo Macports has Python and likely Java; may be an easier way to install, maintain and remove later
clear

WGSEDIR=$(dirname "$0")             # Get the script location to determine install directory
WGSEABS=$(cd "$WGSEDIR"; pwd)       # By cd'ing to it, resolve any aliases and symlinks
WGSEFIN="${WGSEABS}"                # Removed escape any embedded spaces; add trailing slash ${WGSEABS/ /\\ }/
cd "${WGSEFIN}"
# echo '******** WGSEFIN:' "${WGSEFIN}"

echo '======================================================================================================'
echo 'WGS Extract v4 Installer for MacOS'
echo
echo 'WGS Extract needs Apple Xcode CLI, Python3, Macports, Java v11+ and Java v8 on your computer.'
echo 'They must install into system directories. So a password is usually required now.'
echo

# Common environment setup for scripts here; sets some variables used later
declare cpu_arch
declare osrelease
source scripts/zcommon.sh dummy

#---------------------------------------------------------------------------------------------------------------
#VERS=`defaults read loginwindow SystemVersionStampAsString`
BVERS=$(sw_vers -buildVersion)
IFS="." read -r -a OVERS <<< "${osrelease}"
if (( OVERS[0] > 10 )) ; then
  echo "MacOSX Version: ${OVERS[0]}.${OVERS[1]}+build${BVERS} on ${cpu_arch}"
else
  echo "MacOSX Version: ${OVERS[0]}.${OVERS[1]}.${OVERS[2]}+build${BVERS} on ${cpu_arch}"
fi
echo

echo '======================================================================================================'
# Install Xcode CLI tools first as some Macports and Pythin PIP tool installs need it
if ! (xcode-select -p &> /dev/null); then
  echo 'Installing Apple Xcode Command Line tools; needed for MacPorts packages and Python modules'
  # sudo xcode-select --install.  # only schedules installer; below touch does the same
  touch /tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress;
  PROD=$(softwareupdate -l | \
      grep "\*.*Command Line" | \
      head -n 1 | awk -F"Command" '{print $2}' | \
      sed -e 's/^ *//' | \
      tr -d '\n')
  softwareupdate -i "Command $PROD" --verbose;
else
  echo 'Apple Xcode CLI tools already installed'
fi
echo

echo '======================================================================================================'
if [ ! -f /usr/local/bin/python3 ]; then
  echo 'Installing Python 3.11.0 (ignoring labotomized Python in MacOS).'
  case "$cpu_arch" in
    x86_64*)  PYTHONF=python-3.11.0-macos11.pkg  ;;   # No longer separate arch installs needed
    arm*)     PYTHONF=python-3.11.0-macos11.pkg  ;;   # Universal installer with native code per arch
    *) echo "*** Error: Unknown MacOS Architecture ${cpu_arch}" ; exit 1 ;;
  esac
  PYTHONP="https://www.python.org/ftp/python/3.11.0/${PYTHONF}"
  # shellcheck disable=SC2154
  ${curlx} -o "$PYTHONF" "$PYTHONP"
  sudo -H installer -package "$PYTHONF" -target /
  rm -f "$PYTHONF"
else
  echo 'Python3 already installed.'
  # Todo Need to check Python3 version and upgrade if needed
fi
echo

echo '======================================================================================================'
# Apple has shadow-execs so /usr/bin/java is fake; need java --version to truly check is working
# But command -v java -v returns true but is not a jvm; instead a dummy executable by Apple
# Need to avoid Oracla Java as it requires interactive license pop-up during install
# Azul is the only site with M1 and x86 code bases; so utilize them.
if ! (ls /Library/Java/JavaVirtualMachines/zulu-17.jre >/dev/null 2>&1 ); then
  echo 'Installing Azul Java JRE 17 LTS ...'
  case "$cpu_arch" in
    x86_64*)  JAVAF="zulu17.32.13-ca-jre17.0.2-macosx_x64"  ;;
    arm*)     JAVAF="zulu17.32.13-ca-jre17.0.2-macosx_aarch64"  ;;
    *) echo "*** Error: Unknown MacOS Architecture ${cpu_arch}" ; exit 1 ;;
  esac
  JAVAP="https://cdn.azul.com/zulu/bin/${JAVAF}.tar.gz"
  ${curlx} -o azul17.tgz "$JAVAP"
  #/opt/local/bin/7z x -tzip -y azul17.tgz >/dev/null && rm -f azul17.tgz
  tar xf azul17.tgz
  sudo -H mv -f "${JAVAF}/zulu-17.jre" /Library/Java/JavaVirtualMachines
  rm -rf "$JAVAF" azul17.tgz
  echo ' ... finished Java v17 install.'
else
  echo 'Java v17 already installed.'
fi

if ! (ls /Library/Java/JavaVirtualMachines/zulu-8.jre >/dev/null 2>&1 ); then
  echo 'Installing Azul Java JRE 8 LTS ...'
  case "$cpu_arch" in
    x86_64*)  JAVAF="zulu8.64.0.15-ca-fx-jre8.0.342-macosx_x64"  ;;
    arm*)     JAVAF="zulu8.64.0.15-ca-fx-jre8.0.342-macosx_aarch64"  ;;
    *) echo "*** Error: Unknown MacOS Architecture ${cpu_arch}" ; exit 1 ;;
  esac
  JAVAP="https://cdn.azul.com/zulu/bin/${JAVAF}.tar.gz"
  ${curlx} -o azul8.tgz "$JAVAP"
  #/opt/local/bin/7z x -tzip -y azul8.tgz >/dev/null && rm -f azul8.tgz
  tar xf azul8.tgz
  sudo -H mv -f "${JAVAF}/zulu-8.jre" /Library/Java/JavaVirtualMachines
  rm -rf "$JAVAF" azul8.tgz
  echo ' ... finished Java v8 install.'
else
  echo 'Java v8 already installed.'
fi
echo


echo '======================================================================================================'
# Now install MacPorts as a way to get additional Unix executables and many of the bioinformatic tools
LOCALPORT=/opt/local/bin/port
if [ ! -f $LOCALPORT ]; then
  ver="2.9.1"   # Required for Sonoma
  echo "Installing MacPorts v${ver} on MacOS version ${OVERS[0]}.${OVERS[1]}."
  case ${OVERS[0]} in
    10)   # Original MacOSX on x86_64 platform
      case ${OVERS[1]} in
        15) MACPORTSF="MacPorts-${ver}-10.15-Catalina.pkg";;
        14) MACPORTSF="MacPorts-${ver}-10.14-Mojave.pkg";;
        13) MACPORTSF="MacPorts-${ver}-10.13-HighSierra.pkg";;
        12) MACPORTSF="MacPorts-${ver}-10.12-Sierra.pkg";;
         *) Echo "***ERROR: Do not understand MacOS version 10.${OVERS[1]} to install MacPorts software"; exit 1;;
      esac;;
    11)  # BigSur is v11.x now.  Minor versions are patches.  @ 11.3 as of May 2021
      MACPORTSF="MacPorts-${ver}-11-BigSur.pkg";;
    12)  # Monterey is v12.x now. Minor versions are patches.  @ 12.4 as of May 2022
      MACPORTSF="MacPorts-${ver}-12-Monterey.pkg";;
    13)  # Ventura is v13.x now.
      MACPORTSF="MacPorts-${ver}-13-Ventura.pkg";;
    14)  # Sonoma is v14.x now.
      MACPORTSF="MacPorts-${ver}-14-Sonoma.pkg";;
    *)
      Echo "***ERROR: Do not undestand MacOS Version: ${OVERS[0]}.${OVERS[1]} to install MacPorts software"; exit 1;;
  esac
  MACPORTSP="https://github.com/macports/macports-base/releases/download/v${ver}/${MACPORTSF}"
  # MACPORTSP="https://distfiles.macports.org/MacPorts/${MACPORTSF}"
  ${curlx} -O "$MACPORTSP"  # MacOS has issues with MacPorts SSL; so curl -k.
  # See https://github.com/WGSExtract/WGSExtract.github.io/discussions/9#discussioncomment-2923858
  sudo -H installer -package "$MACPORTSF" -target /
  rm -f "$MACPORTSF"
else
  echo 'Updating Macports and its package directory ...'
  sudo -H $LOCALPORT selfupdate
  sudo -H $LOCALPORT upgrade outdated
fi
echo

# Filter prolific warnings going to stderr but meant for developers (on MacOS 13 Ventura with new 2.8.0)
# Then filter out detailed stdout messages except Configuring; want to leave unusual / unexpected messages but ...
if [ ! -f /opt/local/bin/7zz ]; then
  echo 'Macports is installing needed Unix utilities and dependencies (can take one hour) ...'
  sudo -H $LOCALPORT -N install bash grep gsed coreutils zip unzip 7zip md5sha1sum jq \
   2> >(grep -v "^Warning: Configuration logfiles" 1>&2) | grep "^--->  Configuring"
  sudo -H ln -s /opt/local/bin/7zz /opt/local/bin/7z    # New 7zip only installs as 7zz
fi
if [ ! -f /opt/local/bin/samtools ]; then
  echo 'MacPorts is installing needed bioinformatics packages and dependencies (can take one hour)...'
  sudo -H $LOCALPORT -N install samtools bcftools htslib \
   2> >(grep -v "^Warning: Configuration logfiles" 1>&2) | grep "^--->  Configuring"
fi
echo 'Finished installing MacPorts base with needed bioinformatic packages.'

# BWA is not in MacPorts! BWA is on homebrew if we decide to switch
# So grab compiled versions from https://github.com/smikkelsendk/bwa-for-arm/tree/master/bin
if [ ! -f /opt/local/bin/bwa ]; then
  echo 'Adding BWA ...'
  case "$cpu_arch" in
    x86_64*)  BWAF=bwa0717-mac-x64  ;;
    arm*)     BWAF=bwa0717-mac-aarch64  ;;
    *)        echo "*** Error: Unknown MacOS Architecture ${cpu_arch}"  ;  exit  ;;
  esac
  BWAP="https://raw.githubusercontent.com/smikkelsendk/bwa-for-arm/master/bin/${BWAF}.tar.gz"
  # BWAP="https://github.com/smikkelsendk/bwa-for-arm/raw/master/bin/${BWAF}.tar.gz"
  ${curlx} -o bwa.tgz "$BWAP"
  if [ -e bwa.tgz ]; then
    tar xf bwa.tgz && rm -f bwa.tgz
    chmod +x $BWAF
    sudo -H mv -f $BWAF /opt/local/bin/bwa    # Not part of MacPorts but put it there for convenience
  else
    echo '***ERROR downloading compiled BWA from smikkelsendk.'
  fi
fi
echo

# Todo bowtie2 -- in Macports but errors on install
# Todo bwa-mem2, minimap2 and fastp install; not in Macports nor Homebrew. Only Linux release on github master.
# Todo HiSat2.2.1 OSX x86 64bit  https://cloud.biohpc.swmed.edu/index.php/s/zMgEtnF6LjnjFrr/download
# Todo pbmm2, a PacBio Minimap2 front-end  https://github.com/PacificBiosciences/pbmm2

echo '======================================================================================================'
echo 'Calling Common script to finish WGS Extract install ...'
echo /opt/local/bin/bash "${WGSEFIN}/scripts/zinstall_common.sh" "$cpu_arch" "$osrelease"
/opt/local/bin/bash "${WGSEFIN}/scripts/zinstall_common.sh" "$cpu_arch" "$osrelease"
status=$?

if [[ $status -eq 0 ]]; then
  echo '======================================================================================================'
  echo
  echo 'Congratulations!  You finished installing WGS Extract v4 on Apple MacOS!'
  echo ' You can start WGS Extract v4 by clicking the WGSExtract.command file. For convenience, make an'
  echo '  alias of it, rename to WGSExtract, and place it in your Application folder or on your desktop.'
  echo
elif [[ $status -eq 10 ]]; then
  exit  # exit silently as was a restart of the Install script
else
  echo '======================================================================================================'
  echo
  echo 'Sorry. Appears there was an error during the WGS Extract v4 install on Apple MacOS.'
  echo ' Please scroll back through the command log and look for any errors.'
  echo
fi
read -n1 -r -p 'Press any key to close this window (after first scrolling up to review for errors) ...'
