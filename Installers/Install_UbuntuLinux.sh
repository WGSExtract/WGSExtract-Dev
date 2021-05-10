#!/bin/bash
#
# WGS Extract Install Script for Ubuntu Linux (and similar?)
#


clear
echo
echo "This installer has only been tested on Ubuntu. If you are using a different distribution,"
echo "then you might have to adapt this script first."
echo

echo "*** Update or add or verify some basic Unix utilities we will need."
sudo apt-get install sed coreutils zip unzip bash grep

echo "*** Install needed bioinformatics packages like samtools if missing."
sudo apt-get install samtools bcftools tabix # dos2unix  # bwa, 7z, curl handled in upgrade script below

echo "Need to make sure Python3 and its support libraries are there."
sudo apt-get install python3 python3-pip python3-tk python3-pil python3-pil.imagetk

# Now assure Java and any Java based applications (note Haplogrep and GATK are JVM based)
sudo apt-get install openjdk-11-jre  # picard-tools

echo "Testing if any Python3 modules are missing. If that's the case, they will be installed..."
# Pandas needed by yleaf
pip3 install pyliftover Pillow pyscreenshot pandas --no-warn-script-location --disable-pip-version-check # biopython pandas
echo

echo "Handling WGS Extract specific installs via the Upgrade script"
bash Upgrade_v2tov3.sh        # Handles reference genome library and BWA that is new to v3

echo ""
echo "Congratulations!  You finished installing WGS Extract v3!"
echo " You can start WGS Extract v3 by clicking the Start_UbuntuLinux.sh. Make a softlink, rename it to "
echo " WGSExtract.sh, and place it on your desktop for convenience in starting the program."
