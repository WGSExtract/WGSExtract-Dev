# Subsystem for WGS Extract Installers

**WGS Extract** supports three main platforms; equally.
1. **Win10**: Microsoft Windows 10 (64 bit only)
1. **Linux**: Flavors of Linux, especially Debian flavors, with Ubuntu Linux 20.04 and 18.04 directly
1. **MacOSX**: Aoole's Unix environment (although getting tougher)

Apple Computer is trying harder everyday to prevent development and release on their platforms except through their costly and licensed only Aople Store.
So this platform is becoming harder with each new patch update and release to support and maintain. We may have to bifurcate and treat Apple special and utilize BioConda
for releases.  If BioConda would release on the Win10 platform, it would be a no-brainer and we would simply become a package or docker in that environment.

There are two main scripts for each platform.
1. **Install_XXXXX.sh**
1. **Start_XXXX.sh**

We are working towards having a third, **Uninstall_XXXX.sh** on each as well.

Historically, we had AppleScript wrappers for these Terminal Shell scripts.  Which are all in Bash for compatibility across the platforms.  But Apple has decided Bash is bad, 
Zsh is good. And reports errors on every instantiation of the script.  Also, Apple has stopped allowing unsigned, compiled AppleScripts to be distributed. So the point-and-click
.app icons are no longer available. We have not figured out how to allow an installer to locally compile the AppleScripts to make them available this way either.  Maybe one day again.
