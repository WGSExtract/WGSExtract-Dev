#!/bin/sh
# WGS Extract MacOSX Install Script
# Copyright (C) 2018-2020 Marko Bauer
# Copyright (C) 2020 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
#
# Todo Need automatic path to get the latest version of Python and MacPorts; version is hardwired in.
export TERM=xterm	# Needed when run from Applescript
clear
echo "WGSExtract needs Python3 and Macports to be installed on your computer."
echo "It now checks if they and all the software they maanage and we need is installed."
echo "Should run this script with 'sudo -H' at the start. If not, will require a password now"
echo

echo "Checking if Python3 is installed..."
FILE="/usr/local/bin/python3"
if [ ! -f $FILE ]; then
    PYTHON="python-3.8.0-macosx10.9.pkg"
    curl -O https://www.python.org/ftp/python/3.8.0/$PYTHON
    sudo -H installer -package $PYTHON -target /
    rm -f $PYTHON
fi
echo

# Sudo not really needed but need to override -H if not set by running script itself
echo "Testing if any Python3 modules are missing. If that's the case, they will be installed..."
sudo -H /usr/local/bin/python3 -m pip install pillow tk pyliftover pandas
echo "Ignore any version warning from PIP. Do NOT upgrade pip as it suggests!"
echo

echo "Checking if macports is installed..."
FILE=/opt/local/bin/port
if [ ! -f $FILE ]; then
    VERS=`defaults read loginwindow SystemVersionStampAsString`
    case $VERS in
        10.15*) MACPORTS="MacPorts-2.6.2-10.15-Cataline.pkg";;
        10.14*) MACPORTS="MacPorts-2.6.2-10.14-Mojave.pkg";;
        10.13*) MACPORTS="MacPorts-2.6.2-10.13-HighSierra.pkg";;
        10.12*) MACPORTS="MacPorts-2.6.2-10.12-Sierra.pkg";;
    *) Echo "Do not understand MacOSX version to install appropriate MacPorts software"; exit 1;;
    esac
    curl -O https://distfiles.macports.org/MacPorts/$MACPORTS
    sudo -H installer -package $MACPORTS -target /
    rm -f $MACPORTS
fi
echo

echo "Macports will now try to install packages like samtools, etc if they are missing."
sudo -H /opt/local/bin/port selfupdate
sudo -H /opt/local/bin/port install samtools bcftools htslib gsed coreutils dos2unix zip unzip bash grep

echo "Starting WGSExtract ... just to verify everything is OK.  Exit immediately if it starts OK"
WGSEDIR=`dirname $0`
/usr/local/bin/python3 $WGSEDIR/programs/wgsextract/wgsextra
