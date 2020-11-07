# WGS Extract Betav2b MacOSX Patch File Instructions

Overlay the [MacOSX patch file](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/WGSExtract_MacOSX_Patch.zip) onto your current [WGS Extract Betav2b release from 18 February 2020](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/docs/README.md). For MacOSX only installations. This is a new **Install_MacOSX** app for MacOSX to do a proper environment install before you run the **WGS Exract** program. And a new, simple **Start_MacOSX** command to start **WGS Extract** more reliably. These are APPS that provide a GUI window click startup in the MacOS environment.

1. Click [here to download the MacOSX patch file](https://github.com/WGSExtract/WGSExtract-Dev/blob/master/WGSExtract_MacOSX_Patch.zip). 
Your MacOSX will ask if it is OK to download from the github.com site.  Say OK.
2. Go to your Downloads directory and see the already unzipped archive folder there with 4 items inside (if not, click to unzip)
3. Drag those four items to your **WGSExtract** program top-level install directory.  There should be 3 **xxx_START.sh** files there already. For clarity, we assume below you have installed the **WGS Extract program** in your home directory of **~/WGSExtract**.
4. Depending on the unzip program used and method of copy, the executable permission may have been stripped from the files.  If so, via a Terminal command line, execute the following commands:  
   **chmod ugo+x Install_MacOSX.app/Contents/MacOS/applet**  
   **chmod ugo+x Start_MacOSX.app/Contents/MacOS/applet**
5. A reminder (see User Manual) that the MacOSX installation cannot reside in the Documents or Downloads folder due to Gatekeeper.  Best to simply put in your home folder. Read this page completely for other helpful hints.

**`Note`**: as of MacOSX Catalina 15.0.7 update, they have tightened Gatekeeper restrictions to not allow pre-compiled Applescripts to be delivered in unsigned code releases. This functionality of the .app folders is no longer availabe there.

Double click on **Install_MacOSX.app** to run the installer (this was formerly the **MacOSX_START.sh**). If you prefer, you can run the Installer from the Terminal command line via the command: **cd ~/WGSExtract; /bin/sh Install_MacOSX.sh** (you must use this Terminal method for Cataline 15.0.7 and later)

Double click on **Start_MacOSX.app** to run the **WGS Extract** program itself (formerly in **MacOSX_START.sh** also).  If you prefer, you can run the program from the Terminal command line via the command: **cd ~/WGSExtract; /bin/sh Start_MacOSX.sh** (you must use this terminal method for Catalina 15.0.7 and later)

You can create an alias of the **Start_MacOSX.app** file (right click, "Make Alias"), rename it (say to **WGSExtract**), 
and then move it to wherever you want to have an icon available to simply click and start the program.  Often the Applications
folder in Finder or maybe the ~/Desktop. (This no longer works in the MacOSX Cataline 15.0.7 or later update.)

The first time you run these apps, it will ask for permission to run **Terminal** inside them.  We use this so you have a command script log of what is happening.  For the Installer, you need to enter your password so the program can have elevated privileges to install programs like Python3 in the system area.  The installer acts as both an installer and updater.  You can run it again in the future to make sure the dependent tools are up to date. (It will not update **WGS Extract** itself, just the dependent tools.)

For compatability across all the recent MacOSX versions, we have made the underlying shell scripts rely on /bin/sh.  With MacOSX 10.15 Catalina, Apple has chosen to now give a warning when using any shell other than zsh. Please ignore this.  Everything is fine with /bin/sh and this provides compatability across release platforms. As it is, the warning is really about Apple changing the default from BASH to Zsh. It would be nice if Apple presented this warning only at login (or at least only once) and not every time you start a BASH shell script.

Since the initial MacOS release back in December, 2019, the **MacOS_start.sh** script has been broken. 
While we documented the work-around in the Users Manual, it has been a silent show stopper for
many. This patch release finally resolves this while waiting for the next release to include it as part of the base program.
