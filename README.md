# WGS Extract Development
>Note: This repository is woefully out of date.  Use the release packages as they contain a complete copy of all build and release files needed to create a new release. See the latest-release JSON file for the URLs to the current 8 packages. Note that the installer deletes content included in the packages during install.  So you must download and overlay the package contents to get a full development release.  Also, as currently setup, the make_release.sh expects to have the working release folders as setup by the installer. Key being cygwin, msys2 and similar. So better yet is to install the release and then overlay the installer packages on top of that.

This is the main **WGS Extract** Developers repository and support area for coders and others with a deeper understanding of the tool.  Beside the actual code repository that can be git cloned to your development desktop, there are additional documentation and build areas.  See the **WGS Extract** [Developers' Wiki](https://github.com/WGSExtract/WGSExtract-Dev/wiki/) for more information.

See [WGS Extract Users](https://wgsextract.github.io/) for the main entry point for users. That is where to download a release and find documentation related to using the tool. 

See also [WGS Extract Historical](https://github.com/WGSExtract/WGSExtract-Historical "WGSE Historical") for an archive of the original work and tool in 2019 by Marko before others started to get involved and help build on his initial body of work.

Within this area are the source code, build scripts and developer documentation for the tool. Including historical release notes.

# See Also (later developments; mostly with heavy AI assist):
* Mac's new GUI generator to put the GUI front-end onto the new CLI back-end. Also formal, signed installers. [WGSExtract-GUI](https://theontho.github.io/gui-for-cli/)
* James Kane's rebirth of dna-warehouse with a self-contained, RUST-based, single executable, desktop [Navigator tool](https://jameskane.blog/genomics/2026/07/11/navigator-alpha-run-your-own-genomics-lab.html).
* Mac's refactor into a CLI-first form using AI tools: [WGSExtract-CLI](https://github.com/theontho/wgsextract-cli)
  * Simplifying cross-platform support using Pixi; adding development out of GIT and GIT support for regression testing
* Mac's minor restructure: [WGSExtractRepo](https://github.com/theontho/WGSExtractRepo?tab=readme-ov-file)
  * to allow a true development and relesae out of the GITHub source base
  * Includes work on AI rewrites of the installer and support shell code into Python (potentially getting to a Python Package only release tool)
* Teemu's AI tool refactoring and porting to pyQT: [WGSExtractv6](https://github.com/teepean/WGSExtractv6)  (stub code in many places; a forked and restart project; not extension of effort here)
* Marco's original refactoring (called [WGSE-NG](https://github.com/WGSE-NG/)) but now a rewrite ([double-helix](https://doublehelix.app)). 
  * Lots of work for 6 months but never completed. 
  * See [WGSE-NG/External](https://github.com/wgse-ng/External/) for original work to move msys2 windows binaries creation into GIT YAML script and [attempt to get htslib folks to adopt](https://github.com/samtools/samtools/issues/2064).
