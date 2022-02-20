# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.import os

import sys      # for argv[] if called as separate process / program
import os.path  # for os.remove, os.path.exists

from pyliftover import LiftOver     # Main heavyweight lifter

#  Needed WGSE environment as if sub-module; but will call wgse.init() in __main__ if called standalone
import settings as wgse
font = wgse.font
from utilities import DEBUG, nativeOS, universalOS
from commandprocessor import run_bash_script

"""###################################################################################################################
  Wrapper for PyLiftover which is itself based on the UCSC Liftover code and liftover reference file.
    Called by module microarray on BAM extracted CombinedKit that was based on reference model Build38.
    Replaces file it is called with (inplace); often an ~50 MB, 2+million SNP file is extracted from a WGS BAM.
    Although currently called as stand-alone, was recreating and using some settings like WGS Extract (e.g. temp dir)
    so pulled in wgse and similar support modules now to make use of those features.
"""
# Todo need to inline this code with main microarray.py file so does not create new environment


def my_hg38tohg19(combined_file_oFPB, ref="hg38"):
    """ Glue and single call for current use based on HG38 to HG19 liftover of CombinedKit TSV file. """

    # Setup main and intermediate file names (should have been passed the CombinedKit TSV
    combined_txt_oFN        = f'{combined_file_oFPB}.txt'         # Input and eventual output (replaced)
    combined_tmp_oFN        = f'{combined_file_oFPB}.tmp'         # Intermediate; after liftover
    combined_sorted_txt_oFN = f'{combined_file_oFPB}-sorted.txt'  # Intermediate; after liftover then sort

    combined_file_FPB       = universalOS(combined_file_oFPB)     # Prepare for quoted, universalOS versions
    combined_txt_qFN        = f'"{combined_file_FPB}.txt"'
    combined_tmp_qFN        = f'"{combined_file_FPB}.tmp"'
    combined_sorted_txt_qFN = f'"{combined_file_FPB}-sorted.txt"'

    if not os.path.exists(combined_txt_oFN):
        print(f"*** Fatal error. File to liftover does not exist: {combined_txt_oFN}")
        exit()

    print("Converting Microarray Build38 to HG19 positions to maintain compatibility...")
    # Instantiates Liftover class from PyLiftover  # Todo could auto download reference in pyliftover; or reference library?
    lo = LiftOver(wgse.reflib.liftover_chain_oFN)

    # Do liftover with combined_txt as input, combined_tmp as output
    with open(combined_tmp_oFN, "wb") as f_sink:
        with open(combined_txt_oFN, "r") as f_source:
            for source_line in f_source:
                line_tabs = source_line.strip().split("\t")     # strip leading white space; split by tabs
                if '#' != line_tabs[0][0]:      # Ignore comment lines / header (assumes # as first character)
                    line_tabs = source_line.split("\t")
                    oldid = line_tabs[0]
                    oldchrom = ("chr" + line_tabs[1]).replace("chrMT", "chrM")
                    oldpos = line_tabs[2]
                    resultforsnp = line_tabs[3]

                    # Actual call to pyliftover to convert positions
                    newcoord = lo.convert_coordinate(oldchrom, int(oldpos))

                    if newcoord and len(newcoord) > 0:
                        newchrom = newcoord[0][0]
                        newpos = newcoord[0][1]
                        if newchrom in wgse.valid_chromosomes:  # Only liftover the primary 25 sequences
                            newchrom = newchrom.replace("chrM", "chrMT").replace("chr", "")
                            f_sink.write(f'{oldid}\t{newchrom}\t{str(newpos)}\t{resultforsnp}'.encode())
                        else:
                            DEBUG(f'Invalid chromosome from liftover: old: {oldchrom}, {oldpos}; new: {newchrom}, {newpos}')
                    else:
                        DEBUG(f'Cannot liftover HG38 position (ignoring): {oldid}, {oldchrom}, {str(oldpos)}, {resultforsnp}')


    # Going to overwrite the original input file so just delete now
    os.remove(combined_txt_oFN) if os.path.exists(combined_txt_oFN) else ''         # Final file to use

    # Sort the resultant output and pre-pend with a header before making final CombinedKit file
    # Todo Could a more efficient sort be done with BCFTools on the TSV file?

    header_qFN = f'"{wgse.microarray_FP}raw_file_templates/head/23andMe_V3.txt"'

    command  = (
        f"{wgse.sortx_qFN} -t $'\\t' -k2,3 -V {combined_tmp_qFN} > {combined_sorted_txt_qFN}\n"
        f"{wgse.catx_qFN} {header_qFN} {combined_sorted_txt_qFN} > {combined_txt_qFN}\n"
    )

    run_bash_script("LiftoverCleanup", command)

    os.remove(combined_tmp_oFN)
    os.remove(combined_sorted_txt_oFN)

    # Final liftover result is in original source file; replaced the content in combined_txt_oFN


if __name__ == '__main__':
    """ hg38tohg19 liftover MAIN; for when run directly  """
    wgse.init(False)        # Setup local WGSE environment when run stand-alone

    # Process command line arguments
    combined_file_oFPB = nativeOS(sys.argv[1])          # passed as FPB; note in Output and not Temp directory
    source_ref         = sys.argv[2].strip()

    my_hg38tohg19(combined_file_oFPB, source_ref)

#
# Todo how is this liftover working for the Build38 files supplied it?  Should be broken ...
#  Code in original WGSE Beta v2 release always added "chr" to the original name and then stripped it from the final
#  So should only work if passing in GRCh38 as liftover template file is HG38 input
#  And then always generated GRCh37 named chromosomes after doing liftover to Build19
#  The source_ref is ignored in the original; made it more explicit and used here.
#  The only time this is called is to liftover a Build38 file to HG19 naming; but need GRCh37 naming for final?
#  Was working as written? If pass an HG38 original, you would get "chrchr1" as the oldname which would not match
#  Need to investigate more before correcting code here and declaring GRCh38 works as a source file (never tried)
#  Did create and add to code the use of the template for GRCh38 naming for the mpileup/call/annotate function
#