# WGS Extract Betav2b MacOSX Patch File Instructions

This is describing how to overlay the [MacOSX patch file](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/WGSExtract_MacOSX_Patch.zip) into your current [Betav2b release from 18 February 2020](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/docs/README.md). For MacOSX only installations. We rovide a new Install_MacOSX.sh script to do a proper environment install before you run the **WGS Exract** program. And a new, simple Start_MacOSX.sh command to start **WGS Extract** itself. Also include are two APPS that provide GUI window click startup of both these scripts.

1. Click [here to download the MacOSX patch file](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/WGSExtract_MacOSX_Patch.zip).
Your MacOSX will ask if it is OK to download from the github.com site.  Say OK.
1. Go to your Downloads directory and see the already unzipped archive folder there with 4 items inside.
2. Drag those four items to your **WGSExtract** program top-level install directory.  
There should already be 3 **xxx_Start.sh** files there.

Simply double click on **Install_MacOSX.app** to run the installer (this was formerly the **MacOSX_START.sh**)

Simply double click on **Start_MacOSX.app** to run the **WGS Extract** program itself.

You can create an alias of the **Start_MacOSX.app** file (right click, "Make ALias"), rename it (say to **WGSExtract**), 
and move it to wherever you want to have an icon available to simply click and start the program.  Often the Applications
folder in Finder or maybe the desktop).

Since the initial MacOS release back in December, 2019, the **MacOS_start.sh** script has been horribly broken. 
While we documented it and the fixes in the Users Manual, it has been a silent show stopper for
many we fear.  While we fixed it just after the v2B release in February, no one complained so we never went
forward with this patch.  We just assumed there were no MacOSX users. Guess it was more the case of those who
tried and simply gave up without a word. Our aplogies for that.
