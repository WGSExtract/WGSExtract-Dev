#!/bin/sh
clear
echo "WGSExtract needs Python3 and Macports to be installed on your computer."
echo "It now does check if all the software that it needs is installed, before it launches."
echo

echo "Checking if Python3 is installed..."
FILE="/usr/local/bin/python3"
if [ ! -f $FILE ]; then
#    echo "Python3 not found in directory /usr/local/bin/"
#    echo "Please download and install this file with standard settings:"
#    echo "https://www.python.org/ftp/python/3.8.0/python-3.8.0-macosx10.9.pkg"
#    echo "After that, start WGSExtract again."
#    exit 1
    curl -O https://www.python.org/ftp/python/3.8.0/python-3.8.0-macosx10.9.pkg
    sudo installer -package python-3.8.0-macosx10.9.pkg -target /
    rm -f python-3.8.0-macosx10.9.pkg
fi
echo "Found Python3 --> OK!"
echo

echo "Testing if any Python3 modules are missing. If that's the case, they will be installed..."
/usr/local/bin/python3 -m pip install pillow tk pyliftover pandas
echo "If pip showed you a warning, then please ignore that. Do NOT upgrade pip as it suggests!"
echo

echo "Checking if macports is installed..."
FILE=/opt/local/bin/port
if [ ! -f $FILE ]; then
#    echo "Macports not found in directory /opt/local/bin/"
#    echo "Please download and install Macports with standard settings:"
#    echo "https://distfiles.macports.org/MacPorts/MacPorts-2.6.2-10.15-Catalina.pkg"
#    echo "After that, start WGSExtract again."
#    exit 2
    curl -O https://distfiles.macports.org/MacPorts/MacPorts-2.6.2-10.14-Mojave.pkg
    sudo installer -package MacPorts-2.6.2-10.14-Mojave.pkg -target /
    rm -f MacPorts-2.6.2-10.14-Mojave.pkg 
#    curl -O https://distfiles.macports.org/MacPorts/MacPorts-2.6.2-10.15-Cataline.pkg
#    sudo installer -package MacPorts-2.6.2-10.15-Catalina.pkg -target /
#    rm -f MacPorts-2.6.2-10.15-Catalina.pkg 
fi
echo "Found Macports --> OK!"
echo

echo "Macports will now try to download packages like samtools, if they are missing."
while [[ ! ( -f /opt/local/bin/samtools &&  -f /opt/local/bin/bcftools &&  -f /opt/local/bin/tabix &&  -f /opt/local/bin/gsed &&  -f /opt/local/bin/gsort  \
&&  -f /opt/local/bin/gcat &&  -f /opt/local/bin/dos2unix &&  -f /opt/local/bin/zip &&  -f /opt/local/bin/unzip &&  -f /opt/local/bin/bgzip \
&&  -f /opt/local/bin/bash &&  -f /opt/local/bin/ggrep &&  -f /opt/local/bin/htsfile ) ]];do   # -f /opt/local/bin/dot &&
	echo "Missing packages detected. Macports will now try to install them. It will need the administrator password for that."
	sudo /opt/local/bin/port install samtools bcftools htslib gsed coreutils dos2unix zip unzip bash grep # graphviz
done
echo "No missing packages."

echo "Starting WGSExtract..."
/usr/local/bin/python3 ./programs/wgsextract/wgsextract.py
