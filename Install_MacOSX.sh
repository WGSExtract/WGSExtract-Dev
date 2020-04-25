#!/bin/sh
export TERM=xterm	# Needed when run from Applescript
clear
echo "WGSExtract needs Python3 and Macports to be installed on your computer."
echo "It now does check if all the software that it needs is installed, before it launches."
echo

echo "Checking if Python3 is installed..."
FILE="/usr/local/bin/python3"
if [ ! -f $FILE ]; then
    curl -O https://www.python.org/ftp/python/3.8.0/python-3.8.0-macosx10.9.pkg
    sudo -H installer -package python-3.8.0-macosx10.9.pkg -target /
    rm -f python-3.8.0-macosx10.9.pkg
fi
echo "Found Python3 --> OK!"
echo

# Sudo not really needed but need to override -H if not set by running script itself
echo "Testing if any Python3 modules are missing. If that's the case, they will be installed..."
sudo -H /usr/local/bin/python3 -m pip install pillow tk pyliftover pandas
echo "If pip showed you a warning, then please ignore that. Do NOT upgrade pip as it suggests!"
echo

# Todo need to check MacOSX version and get appropriate package
# for example, using defaults read loginwindow SystemVersionStampAsString
echo "Checking if macports is installed..."
FILE=/opt/local/bin/port
if [ ! -f $FILE ]; then
    curl -O https://distfiles.macports.org/MacPorts/MacPorts-2.6.2-10.14-Mojave.pkg
    sudo installer -package MacPorts-2.6.2-10.14-Mojave.pkg -target /
    rm -f MacPorts-2.6.2-10.14-Mojave.pkg 
#    curl -O https://distfiles.macports.org/MacPorts/MacPorts-2.6.2-10.15-Cataline.pkg
#    sudo installer -package MacPorts-2.6.2-10.15-Catalina.pkg -target /
#    rm -f MacPorts-2.6.2-10.15-Catalina.pkg 
fi
echo "Found Macports --> OK!"
echo

echo "Macports will now try to download packages like samtools, etc if they are missing."
sudo -H opt/local/bin/port install samtools bcftools htslib gsed coreutils dos2unix zip unzip bash grep # graphviz

echo "Starting WGSExtract..."
WGSEDIR=`dirname $0`
/usr/local/bin/python3 $WGSEDIR/programs/wgsextract/wgsextract.py
