# Subsystem for WGS Extract Installers

**WGS Extract** supports three main platforms; equally.
1. **Win10**: Microsoft Windows 10 (64 bit only)
1. **Linux**: Flavors of Linux, especially Debian flavors, with Ubuntu Linux 20.04 and 18.04 supported directly
1. **MacOSX**: Apple's Unix environment (although getting tougher)

There are two main scripts for each platform.
1. **Install_XXXXX.sh**
1. **Start_XXXX.sh**

We are working towards having a third, **Uninstall_XXXX.sh** on each as well. Only available for the MacOS at the moment. 

Historically, we had AppleScript wrappers for these Terminal Shell scripts.  The shell scripts are all in Bash for compatibility across the platforms.  But Apple has decided Bash is bad, 
Zsh is good. And reports errors on every instantiation of the Bash script.  Also, Apple has stopped allowing unsigned, compiled AppleScripts to be distributed. So the point-and-click
.app icons are no longer available. We have not figured out how to allow an installer to locally compile the AppleScripts to make them available.  Maybe one day again.

Apple Computer is trying harder everyday to prevent development and release on their platforms except through their licensed Apple Store.
So this platform is becoming harder with each new patch update and release to support and maintain. We may have to bifurcate and treat Apple special and utilize BioConda
for releases.  If BioConda would release on the Win10 platform, it would be a no-brainer and we would simply become a package or docker in that environment. But many who need this tool 
the most are on Win10 platforms.

