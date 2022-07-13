#!/bin/sh
# WGS Extract MacOSX Uninstall Script
# Copyright (C) 2020 Randolph Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
#
export TERM=xterm       # Needed when run from Applescript
clear

# Really hate using all these "sudo rm -rf" in here.  Prone to catastrophic error

echo "WGS Extract installed Python 3.8, MacPorts, Xcode CLI and openJDK along"
echo "with associated modules like samtools, htslib, pyliftover and such."
echo ""
echo "Should run this script with 'sudo -H' at the start. "
echo "If not, will require a password during the run."
echo

echo "Uninstalling Python 3.8 ... (ignoring native Python3 in Catalina)"
if [ -f /Library/Frameworks/Python.framework/Versions/3.8 -o -a /usr/local/bin/python3 ]; then
    read -p "Do you want to remove Python 3.8 [Y/n]? " -n 1 -r ; echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing Python 3.8 libraries and modules ..."
        sudo -H rm -rf /Library/Frameworks/Python.framework/Versions/3.8
        echo " ... and now the Python 3.8 application is being removed"
        sudo -H rm -rf "/Applications/Python 3.8"
        sudo -H rm -rf /usr/local/bin/python3
    fi
fi
echo


# See https://guide.macports.org/chungked/installing_macports_uninstalling.html
echo "Uninstalling MacPorts ..."
if [ -f /opt/local/bin/port ]; then
    # Find it hard to believe that only MacPorts installs in /opt/local
    #  and we just wipe it out entirely
    read -p "Do you want to remove Macports and all its programs [Y/n]? " -n 1 -r ; echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "MacPorts is deleting packages previously installed ..."
        sudo -H /opt/local/bin/port -fp uninstall installed

        echo " ... and now deleting the main MacPorts installation."
        sudo dscl . -delete /Users/macports
        sudo dscl . -delete /Groups/macports
        sudo -H rm -rf /opt/local \
    	    /Applications/DarwinPorts \
    	    /Applications/MacPorts \
	    /Library/LaunchDaemons/org.macports.* \
	    /Library/Receipts/DarwinPorts*.pkg \
	    /Library/Receipts/MacPorts*.pkg \
	    /Library/StartupItems/DarwinPortsStartup \
	    /Library/Tcl/darwinports1.0 \
	    /Library/Tcl/macports1.0 \
	    ~/.macports
    fi
fi
echo


# if we wipe out the folder; as instructed by Apple to do, then install script errors on trying to reinstall
#echo "Removing Xcode CLI ..."
#if [ -d /Library/Developer/CommandLineTools ]; then
#    read -p "Do you want to remove Xcode CLI [Y/n]? " -n 1 -r ; echo
#    if [[ $REPLY =~ ^[Yy]$ ]]; then
#	sudo -H rm -rf  /Library/Developer/CommandLineTools
#    fi
#fi


echo "Removing Java ... "
if [ -d /Library/Java/JavaVirtualMachines/adoptopenjdk-11.jre ]; then
    read -p "Do you want to remove openJDK Java JRE [Y/n]? " -n 1 -r ; echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
	sudo -H rm -rf /Library/Java/JavaVirtualMachines/adoptopenjdk-11.jre
    fi
fi


echo "Finished uninstalling programs needed by WGSExtract. "
echo "*** Please drag the WGS Extract install folder to the trash yourself. "
echo "If you specified a Reference_Genome and/or Temp directory outside the WGS Extract"
echo "install folder, you will need to move those folders to the trash as well."
