#!/bin/bash

#
# WGS Extract Install Script for Ubuntu Linux (and similar?)
#
if [[ $EUID -ne 0 ]]; then
   echo "This script must be run as root" 
   exit 1
fi

clear
echo "WGSExtract needs some external software to be installed on your computer."
echo "It now does check if all the software that it needs is installed before it launches."
echo "This script has only been tested on Ubuntu. If you are using a different distribution,"
echo "then you might have to adapt this script first."
echo

#while [[ ! ( -f "/usr/bin/samtools" &&  -f "/usr/bin/bcftools" &&  -f "/usr/bin/tabix" &&  -f "/usr/bin/sed" &&  -f "/usr/bin/sort"  &&  -f "/usr/bin/cat" &&  -f "/usr/bin/dos2unix" &&  -f "/usr/bin/zip" &&  -f "/usr/bin/unzip" &&  -f "/usr/bin/bgzip" &&  -f "/bin/bash" &&  -f "/usr/bin/dot" &&  -f "/usr/bin/grep" &&  -f "/usr/bin/htsfile" && -f "/usr/bin/python3" && -f "/usr/bin/pip3" ) ]];do
#	echo "Missing packages detected. apt-get will now try to install them."
#	sudo apt-get install samtools bcftools tabix sed coreutils dos2unix zip unzip bash grep python3 python3-pip python3-pydot python3-pygraphviz python3-tk python3-pil python3-pil.imagetk picard-tools
#done
#echo "No missing packages."
echo "*** Update or add some basic Unix utilities we will need."
apt-get install sed coreutils zip unzip bash grep

echo "*** Now download bioinformatics packages like samtools if missing."
apt-get install samtools bcftools tabix dos2unix

echo "Need to make sure Python3 and its support libraries are there."
apt-get python3 python3-pip python3-pydot python3-pygraphviz python3-tk python3-pil python3-pil.imagetk

# Now assure Java and any Java based applications (note Haplogrep and GATK are JVM based)
apt-get install openjdk-8-jre  #picard-tools

echo "Testing if any Python3 modules are missing. If that's the case, they will be installed..."
pip3 install pillow tk pyliftover pandas
echo

echo "Starting WGSExtract to verify everything went OK ..."
python3 ./programs/wgsextract/wgsextract.py