# Main WGS Extract Programs

Historically, this was a collection of a number of directories.  Not just the **WGS Extract** code in a subdirectory here. But we are slowly removing the other directories as needs dictate.

## Directories that are still here (or newly added):
1. **yleaf/** -- the modified copies of the yleaf Python code and its needed template files.
1. **wgsextract/** -- the main python code based and related files unique to this program
1. **jartools/** -- Created by all platform installers; a new entry for the expanded JAR file programs (installed and maintained by the installers; not the code base here).
1. **python/** -- Created by the Win10 installer only for a stand-alone Python3 installation (see Win10 Subsystem)
1. **win10tools/** -- Created by the Win10 installeronly for the stand-alone Bioinformatic tools installation (see Win10 subsystem)

## Historic directories in the release not represented here:
1. **python/** -- the Win10 release stand-alone python release that was distributed to all platforms with the tool. See Win10 Subsystem now as it is only installed on Win10 systems by the installer
1. **samtools-cygwin/**, **samtools-mingw/** -- the Win10 release of stand-aline bioinformatic and *nix tools.  Now a singular win10tools/ folder that is only installed by the Win10 installer.
1. **haplogrep/** -- a folder with the single entry haplogrep.jar.  That and other jar programs now redistributed are in the jartools folder here
1. **Extract23/** -- a misnomer as it had the templates of the microarray system developed in v1/v2. It did have an initial Extract23 template script that was relised on but no longer.
The templates are moved into the WGS Extract program directory as a stand-alone module there.
