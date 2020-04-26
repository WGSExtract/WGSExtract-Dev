# WGS Extract Betav2b MacOSX Patch File Instructions

This is describing how to overlay the [MacOSX patch file](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/WGSExtract_MacOSX_Patch.zip) onto your current [WGS Extract Betav2b release from 18 February 2020](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/docs/README.md). For MacOSX only installations. We rovide a new **Install_MacOSX** app to do a proper environment install before you run the **WGS Exract** program. And a new, simple **Start_MacOSX** command to simply start **WGS Extract** itself (reliably). These are APPS that provide a GUI window click startup.

1. Click [here to download the MacOSX patch file](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/WGSExtract_MacOSX_Patch.zip). Your MacOSX will ask if it is OK to download from the github.com site.  Say OK.
1. Go to your Downloads directory and see the already unzipped archive folder there with 4 items inside.
2. Drag those four items to your **WGSExtract** program top-level install directory.  There should already be 3 **xxx_Start.sh** files there.

Simply double click on **Install_MacOSX.app** to run the installer (this was formerly the **MacOSX_START.sh**)

Simply double click on **Start_MacOSX.app** to run the **WGS Extract** program itself.

You can create an alias of the **Start_MacOSX.app** file (right click, "Make ALias"), rename it (say to **WGSExtract**), 
and move it to wherever you want to have an icon available to simply click and start the program.  Often the Applications
folder in Finder or maybe the desktop).

The first time you run these apps, it will ask for permission to run Terminal inside them.  We do this so you have a command script log of what is happening while they run.  For the Installer, you need to enter your password so the program can have elevated privileges to install programs like Python3 in the system area.  The installer acts as both an installer and updater.  You can run it again in the future to make sure the dependent tools are up to date. (It will not update **WGS Extract** automatically.)

Since the initial MacOS release back in December, 2019, the **MacOS_start.sh** script has been horribly broken. 
While we documented it and the fixes in the Users Manual, it has been a silent show stopper for
many we fear.  While we fixed it just after the v2B release in February, no one complained so we never went
forward with this patch.  We just assumed there were no MacOSX users. Guess it was more the case of those who
tried and simply gave up without a word. Our aplogies for that.
