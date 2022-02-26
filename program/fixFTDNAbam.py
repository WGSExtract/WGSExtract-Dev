# coding: utf8
# Copyright (C) 2022 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""
Standalone script to fix FTDNA BigY BAM files. FTDNA BigY BAM files have a space in the QNAME column
which is not allowed by the SAM spec.  So we change all single spaces found there to an underscore. See
https://gist.github.com/RandyHarr/2cbb8901306891e810d698ab27d212ce for shell script use example to annotate
FTDNA VCF files with yBrowse SNP names and haplogroups.
"""

import sys
import os

# We assume called in a pipeline with samtools -h view used to pass in complete SAM on stdin; filtered SAM on stdout
header = True
for line in sys.stdin:
    if header:                          # Allow spaces in header
        if line[0] == '@':
            os.write(1, str.encode(line))      # print header line; assure \n only even on MS Windows systems
            continue
        else:
            header = False

    # Process every remaining SAM line to remove space from QNAME (first entry)
    # Tabs are used to separate columns.  But space in QNAME field is not allowed.  That is the only white space.
    fields = line.split("\t")     # Split line into column fields using only \t; last field still contains \n
    fields[0] = fields[0].replace(" ",":")
    line = "\t".join(fields)

    os.write(1, str.encode(line))

