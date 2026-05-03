# WGS Extract Development
>Note: This repository is woefully out of date.  Use the release packages as they contain a complete copy of all build and release files needed to create a new release. See the latest-release JSON file for the URLs to the current 8 packages. Note that the installer deletes content included in the packages during install.  So you must download and overlay the package contents to get a full development release.  Also, as currently setup, the make_release.sh expects to have the working release folders as setup by the installer. Key being cygwin, msys2 and similar. So better yet is to install the release and then overlay the installer packages on top of that.

This is the main **WGS Extract** Developers repository and support area for coders and others with a deeper understanding of the tool.  Beside the actual code repository that can be git cloned to your development desktop, there are additional documentation and build areas.  See the **WGS Extract** [Developers' Wiki](https://github.com/WGSExtract/WGSExtract-Dev/wiki/) for more information.

See [WGS Extract Users](https://wgsextract.github.io/) for the main entry point for users. That is where to download a release and find documentation related to using the tool. 

See also [WGS Extract Historical](https://github.com/WGSExtract/WGSExtract-Historical "WGSE Historical") for an archive of the original work and tool in 2019 by Marko before others started to get involved and help build on his initial body of work.

Within this area are the source code, build scripts and developer documentation for the tool. Including historical release notes.

# See Also:
* Teemu's AI tool refactoring and porting to pyQT: [WGSExtractv6](https://github.com/teepean/WGSExtractv6)  (stub code in most places; a forked and restart project; not extension of effort here)
* Mac's restructure to allow a true development and relesae out of the GITHub source base: [WGSExtractRepo](https://github.com/theontho/WGSExtractRepo?tab=readme-ov-file) (likely to replace this repository soon)
  * Includes work on AI rewrites into Python of the installer and support shell code (potentially getting to a Python Package release tool)
* Mac's restructure into a CLI first using Microsoft AI tools: [WGSExtract-CLI](https://github.com/theontho/wgsextract-cli)
