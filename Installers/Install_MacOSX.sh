 #!/bin/sh
# WGS Extract MacOSX Install Script
# Copyright (C) 2018-2020 Marko Bauer
# Copyright (C) 2020 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
#
export TERM=xterm       # Needed when run from Applescript
clear
echo "WGSExtract needs Python3, Macports, and Java to be installed on your computer."
echo "It now checks if they and all the software they manage and we need is installed."
echo "Should run this script with 'sudo -H' at the start. If not, will require a password now"
echo

echo "Installing or updating Python3 ... (ignoring native Python3 in Catalina)"
LOCALP3="/usr/local/bin/python3"
if [ ! -f $LOCALP3 ]; then
    echo " ... Python3 is being retrieved and installed."
    PYTHONF="python-3.8.3-macosx10.9.pkg"
    PYTHONP="https://www.python.org/ftp/python/3.8.3/$PYTHONF"
    curl -O $PYTHONP
    sudo -H installer -package $PYTHONF -target /
    rm -f $PYTHONF
else
    echo " ... Python3 already installed"
    # Todo Need to check Python3 version and upgrade if needed
fi
echo

# Sudo not really needed but need to assert -H if not set when running the script itself
echo "Checking for needed Python3 modules ..."
sudo -H $LOCALP3 -m pip install Pillow pyliftover biopython pandas
echo "IGNORE a version warning from PIP in the above. Do NOT upgrade pip as it suggests!"
echo

echo "Checking Xcode CLI tools ..."
xcode-select -p &> /dev/null
if [ $? -ne 0 ]; then
    echo " ... Installing Xcode Command Line tools needed for MacPorts packages"
    # sudo xcode-select --install.  # only schedules installer; below touch does the same
    touch /tmp/.com.apple.dt.CommandLineTools.installondemand.in-progress;
    PROD=$(softwareupdate -l | \
        grep "\*.*Command Line" | \
        head -n 1 | awk -F"Command" '{print $2}' | \
        sed -e 's/^ *//' | \
        tr -d '\n')
    softwareupdate -i "Command $PROD" --verbose;
else
    echo "... Xcode CLI tools OK"
fi

#VERS=`defaults read loginwindow SystemVersionStampAsString`
PVERS=$(sw_vers -productVersion)
BVERS=$(sw_vers -buildVersion)
OVERS=( ${PVERS//./ } )
echo "MacOSX Version: ${OVERS[0]}.${OVERS[1]}.${OVERS[2]}+build${BVERS}"

echo "Installing or updating MacPorts ..."
LOCALPORT=/opt/local/bin/port
if [ ! -f $LOCALPORT ]; then
    MACPORTSP="https://distfiles.macports.org/MacPorts/"
    case ${OVERS[1]} in
      15) MACPORTSF="MacPorts-2.6.2-10.15-Catalina.pkg";;
      14) MACPORTSF="MacPorts-2.6.2-10.14-Mojave.pkg";;
      13) MACPORTSF="MacPorts-2.6.2-10.13-HighSierra.pkg";;
      12) MACPORTSF="MacPorts-2.6.2-10.12-Sierra.pkg";;
       *) Echo "***ERROR: Do not understand MacOSX version ${OVERS[1]} to install MacPorts software"; exit 1;;
    esac
    echo " ... retrieving and installing MacPorts on MacOSX version ${OVERS[1]}."
    curl -O $MACPORTSP/$MACPORTSF
    sudo -H installer -package $MACPORTSF -target /
    rm -f $MACPORTSF
else
    echo "Macports is already found but now upgrading itself and its package directory ..."
    sudo -H $LOCALPORT selfupdate
fi
echo

echo "Macports will now try to install or update packages (samtools, etc) if they are missing."
LOCALSAM=/opt/local/bin/samtools
if [ ! -f $LOCALSAM ]; then
    echo " ... MacPorts is retrieving needed packages and installing."
    sudo -H $LOCALPORT -N install samtools bcftools htslib gsed coreutils dos2unix zip unzip bash grep
else
    echo " ... MacPorts is upgrading already installed packages; as needed"
    sudo -H $LOCALPORT upgrade outdated
fi

# Need to bypass Oracla Java as it requires interactive license pop-up
# No "latest available" link at OpenJDK so whatever was latest version at release time
# Apple has pseudo installs so /usr/bin/java not enough; need java --version to check
echo "Installing Java ..."
LOCALJAVA=/Library/Java/JavaVirtualMachines/adoptopenjdk-11.jre
JAVAP=https://github.com/AdoptOpenJDK/openjdk11-binaries/releases/download/jdk-11.0.8%2B10
JAVAF=OpenJDK11U-jre_x64_mac_hotspot_11.0.8_10.pkg
if [ ! -d $LOCALJAVA ]; then
    curl -L -O $JAVAP/$JAVAF
    sudo -H installer -package $JAVAF -target /
    rm -f $JAVAF
else
    # Todo Need to check Java version and upgrade if needed
    echo "... nothing to install"
fi

# Todo need to grab needed JAR files (GATK, Haplogrep?)

echo "Finished installing or upgrading programs needed by WGSExtract."
Echo
echo " You can now start the program by using ./Start_MacOSX.sh or "
echo " clicking the Start_MacOSX.app. Make an alias of the .app, "
echo " rename it to WGSExtract, and place it in your Application "
echo " folder or on your desktop for convenience."
