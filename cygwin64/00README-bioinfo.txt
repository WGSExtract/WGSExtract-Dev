An image of a cygwin64 /usr/local folder containing the Bioinformatic tools
 compiled for the x86_64 Windows 10/11 platform. Main components are htslib,
 samtools, bcftools, and bwa.

Copyright (c) 2022. Randy Harr as part of the WGS Extract release.
See https://wgsextract.github.io/ for license details.

Created as part of the WGS Extract (https://wgsextract.github.io) project
 for the support of manipulating Whole-genome sequence files.

See the make_bioinfo.sh file for more information on what is made and how.

Open Source Licenses for the various tools included here are in the
 directory open_source_licences/.

Should be overlayed in a matched cygwin64 base environment with needed libraries.
Should be extracted within a shell in the target cygwin64 environment using the
 BASH command:
  7z x -tzip -o/usr cygwin64-bioinfo_vN.zip



