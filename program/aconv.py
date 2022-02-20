# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copytight (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import sys
import os.path
import platform

""""###################################################################################################################
    Module aconv (autosomal microarray file converter / creator)

    Currently setup as a standalone called module not relying on wgse system.
    Assumes have already called microarray to create CombinedKit file of all possible SNPs
    This module is then called for each individual file format you wish to create from the CombinedKit TSV table file
    Todo need to inline this code with microarry.py module; can use WGSE Settings module then
"""

# Globals used in this module; leftover from original implementation Todo should remove need for globals here
target_file_without_suffix_all = None
target_file_without_suffix = None
target_type_suffix = None
target_type_name_all = None
target_type_name = None
f_source_file_called = None     # type:
source_file_called = None
f_targetfile = None
num_target_files = None
last_called_line_read = None
called_chrom  = None
called_pos    = None
called_result = None
templ_chrom = None
templ_pos   = None
templ_id    = None
templates_oFP = None
os_slash = None


def get_template_elements(template_line):
    global target_type_name_all
    global templ_chrom, templ_pos, templ_id
    template_line = template_line.replace('"', '')      # Remove double quotes; do not have WGSE utilities function here
    
    if target_type_name_all in ["FTDNA_V1_Affy", "FTDNA_V2",  "FTDNA_V3",  "MyHeritage_V1", "MyHeritage_V2"]:
        line_elements = template_line.split(",")
        templ_chrom = line_elements[1]
        templ_pos = int(line_elements[2])
        templ_id = line_elements[0]
    elif target_type_name_all == "23andMe_SNPs_API":
        line_elements = template_line.split("\t")
        templ_chrom = line_elements[0]
        templ_chrom = templ_chrom.replace("chr","")
        templ_pos = int(line_elements[1])
        templ_id = line_elements[2].strip()   
    elif target_type_name_all in ["LDNA_V1",  "LDNA_V2"] or \
        "Ancestry" in target_type_name_all or "23andMe" in target_type_name_all:
        line_elements = template_line.split("\t")
        templ_chrom = line_elements[1]
        templ_pos = int(line_elements[2])
        templ_id = line_elements[0]


def chrconv(chrom_to_convert):
        global target_type_name_all
        if target_type_name_all in ["Ancestry_V1", "Ancestry_V2"]:
            if chrom_to_convert not in ["MT", "M", "X", "XY", "Y"]:
                chrom_to_convert = int(chrom_to_convert)
            if chrom_to_convert == 23 or chrom_to_convert == 25:
                chrom_to_convert = "X"
            elif chrom_to_convert == 24:
                chrom_to_convert = "Y"
            elif chrom_to_convert == 26:
                chrom_to_convert = "MT"

        if chrom_to_convert in ["MT", "M"]:
            return 23
        elif chrom_to_convert in ["X", "XY"]:
            return 24
        elif chrom_to_convert in ["Y"]:
            return 25            
        else:
            return int(chrom_to_convert)


def pos1_smaller_than_pos2(chr1,pos1,chr2,pos2):
        if chrconv(chr1) < chrconv(chr2):
            return True
        elif chrconv(chr1) > chrconv(chr2):
            return False
        elif chrconv(chr1) == chrconv(chr2):
            return (pos1 < pos2)


def pos1_equal_to_pos2(chr1,pos1,chr2,pos2):
    if chrconv(chr1) != chrconv(chr2):
        return False
    else:
        return (pos1 == pos2)


def pos1_bigger_than_pos2(chr1,pos1,chr2,pos2):
    if chrconv(chr1) > chrconv(chr2):
        return True
    elif chrconv(chr1) < chrconv(chr2):
        return False
    elif chrconv(chr1) == chrconv(chr2):
        return (pos1 > pos2)


def get_next_called_vars():
    global f_source_file_called, called_chrom, called_pos, called_result, last_called_line_read

    called_line = ""
    for called_line in f_source_file_called:
        if '#' in called_line:
            continue                # Skip comment lines
        #print (called_line)
        break
    if called_line:
        line_columns = called_line.split("\t")
        called_chrom = line_columns[1]
        called_pos = int(line_columns[2])
        called_result = line_columns[3].strip()
    else:
        last_called_line_read = True


def write_line_output_file(output_result):
    global f_targetfile, target_type_name_all
    global templ_chrom, templ_pos, templ_id

    if target_type_name_all in ["FTDNA_V1_Affy", "MyHeritage_V2"]:
        output_line = f'"{templ_id}","{templ_chrom!s}","{templ_pos!s}","{output_result}"'       # quote, comma-sep
    elif target_type_name_all == "FTDNA_V2":
        if templ_id == "rs5939319":
            f_targetfile.write("RSID,CHROMOSOME,POSITION,RESULT\n".encode())
        output_line = f'"{templ_id}","{templ_chrom!s}","{templ_pos!s}","{output_result}"'       # quote, comma-sep
    elif target_type_name_all == "FTDNA_V3":
        output_line = f'{templ_id},{templ_chrom!s},{templ_pos!s},{output_result}'               # comma-sep
    elif target_type_name_all == "23andMe_SNPs_API":
        templ_chrom = templ_chrom.replace("M","MT")
        output_line = f'{templ_id}\t{templ_chrom!s}\t{templ_pos!s}\t{output_result}'            # tab-sep
    elif target_type_name_all in ["LDNA_V1", "LDNA_V2"] or "23andMe" in target_type_name_all:
        output_line = f'{templ_id}\t{templ_chrom!s}\t{templ_pos!s}\t{output_result}'            # tab-sep
    elif target_type_name_all == "MyHeritage_V1":
        if output_result == "CT":
            output_result = "TC"
        elif output_result == "GT":
            output_result = "TG"          
        output_line = f'"{templ_id}","{templ_chrom!s}","{templ_pos!s}","{output_result}"'       # quote, comma-sep
    elif "Ancestry" in target_type_name_all:
        if output_result == "--":
            output_result = "00"
        if output_result == "CT":
            output_result = "TC"
        if output_result == "GT":
            output_result = "TG"
        #if output_result == '.':
        #    output_result = "00"
        output_line = f'{templ_id}\t{templ_chrom!s}\t{templ_pos!s}\t{output_result[0]}\t{output_result[1]}' # tab-sep
    else:
        output_line = ''
    f_targetfile.write(f'{output_line}\n'.encode())


def get_target_type_suffix():
    """ Get the file suffix type for the GLOBAL parameter target_type_name_all. """
    global target_type_name_all
    if target_type_name_all in ["FTDNA_V1_Affy", "FTDNA_V2", "FTDNA_V3", "MyHeritage_V1", "MyHeritage_V2"]:
        return ".csv"
    elif target_type_name_all in ["LDNA_V1", "LDNA_V2"] \
     or "23andMe" in target_type_name_all or "Ancestry" in target_type_name_all:
        return ".txt"


def concat_files():
    """ Create target file by concatinating the header and one or more tail files. Using global paramters. """
    global target_type_name_all, target_file_without_suffix_all, target_type_suffix

    # Read header
    template_head_oFN = f'{templates_oFP}head/{base2}{target_type_suffix}'
    with open(template_head_oFN, "r") as f_template_head:
        head = f_template_head.read()

    # Write output; starting with header
    targetfile_oFN = f'{target_file_without_suffix_all}_{target_type_name_all}{target_type_suffix}'
    with open(targetfile_oFN,"wb") as f_targetfile:
        f_targetfile.write(head.encode())

        # For each target file, open to read, concat into targetfile
        for i in range(num_target_files):
            i_str = str(i+1)
            target_type_name = f'{target_type_name_all}_{i_str}'
            target_file_without_suffix = f'{target_file_without_suffix_all}_{i_str}'
            partfilename_oFN = f'{target_file_without_suffix}_{target_type_name}{target_type_suffix}'
            with open(partfilename_oFN, "r") as f_part:
                part = f_part.read()
            os.remove(partfilename_oFN)
            f_targetfile.write(part.encode())


# Heavily relies on global values passed in and out instead of parameters
def convert_adna():
    global f_source_file_called, f_targetfile
    global called_chrom, called_pos, called_result, last_called_line_read
    global templ_chrom, templ_pos, templ_id
    global source_file_called, templates_oFP, os_slash, num_target_files
    global target_file_without_suffix, target_type_name, target_type_suffix
    global target_type_name_all

    # Parameters in are all globals; reset before each call
    # 4 simultaneous open files; so not using with-as to avoid hideous indentation
    targetfile_oFN = f'{target_file_without_suffix}_{target_type_name}{target_type_suffix}'
    f_targetfile = open(targetfile_oFN, "wb")
    if num_target_files == 1:
        template_head_oFN = f'{templates_oFP}head/{target_type_name}{target_type_suffix}'
        with open(template_head_oFN, "r") as f_template_head:
            head = f_template_head.read()
        f_targetfile.write(head.encode())

    f_source_file_called = open(source_file_called, "r")
    last_called_line_read = False
    get_next_called_vars()

    template_body_oFN = f'{templates_oFP}body/{target_type_name}{target_type_suffix}'
    f_template_body = open(template_body_oFN, "r")
    for template_line in f_template_body:
        get_template_elements(template_line)    # Side effect: sets global templ_chrom, templ_pos, templ_id
        while pos1_smaller_than_pos2(called_chrom, called_pos, templ_chrom, templ_pos) \
              and not last_called_line_read:
            get_next_called_vars()
        if last_called_line_read or templ_chrom == 0 or templ_pos == 0 \
              or pos1_bigger_than_pos2(called_chrom, called_pos, templ_chrom, templ_pos):
            write_line_output_file("00" if "Ancestry" in target_type_name_all else "--")
            continue
        elif pos1_equal_to_pos2(called_chrom, called_pos, templ_chrom, templ_pos):
            write_line_output_file(called_result)
            continue

    f_template_body.close()
    f_source_file_called.close()
    f_targetfile.close()


######################################################################################################################
### Start of Main
#     When called stand-alone as separate program; which is the use currently

if __name__ == '__main__':

    # Designed as stand-alone program and called as such.  No tie in to wgse directly. So must recreate some values.
    os_slash = "\\" if (platform.system() == "Windows") else "/"

    install_oFP = (os.path.dirname(os.path.abspath(__file__))).split(f'program')[0]
    templates_oFP = f'{install_oFP}program{os_slash}microarray{os_slash}raw_file_templates{os_slash}'

    # Simple command line argument capture; only an internal script. File paths are quoted in shell and so still here
    target_type_name_all = sys.argv[1]
    source_file_called = sys.argv[2].strip('"')
    target_file_without_suffix_all = sys.argv[3].strip('"')

    target_type_suffix = get_target_type_suffix()

    if target_type_name_all in ["23andMe_V4", "23andMe_V5"]:
        num_target_files = 2
    elif target_type_name_all == "Ancestry_V1":
        num_target_files = 4
    elif target_type_name_all == "Ancestry_V2":
        num_target_files = 5
    elif target_type_name_all == "FTDNA_V3":
        num_target_files = 3
    else:
        num_target_files = 1

    if num_target_files == 1:
        target_file_without_suffix = target_file_without_suffix_all
        target_type_name = target_type_name_all
        convert_adna()
    else:
        base1 = target_file_without_suffix_all
        base2 = target_type_name_all
        for i in range(num_target_files):
            i_str = str(i + 1)
            target_file_without_suffix = f'{base1}_{i_str}'
            target_type_name = f'{base2}_{i_str}'
            convert_adna()
        concat_files()