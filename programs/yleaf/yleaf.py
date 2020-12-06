###!/usr/bin/env python					# MB Jan2020  Windows
# Copyright (C) 2018-2019 Diego Montiel Gonzalez
# Erasmus Medical Center
# Department of Genetic Identification
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

# YLeaf detection of Y-Haplogroups in Human DNA

import shutil                           # REH 14Mar2020
import time
import subprocess
import string
import random
import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
import collections
import operator
import gc
import platform							# MB Jan2020  Windows

pd.options.mode.chained_assignment = None  # default='warn'

os_slash = "/"							# MB Jan2020  Windows

def get_arguments():
    parser = ArgumentParser()

    parser.add_argument("-fastq", "--Fastq",
            dest="Fastq", required=False,
            help="Use raw FastQ files", metavar="PATH")

    parser.add_argument("-bam", "--Bamfile",
            dest="Bamfile", required=False,
            help="input BAM file", metavar="PATH")

    parser.add_argument("-f", "--fasta-ref",  dest="reference",
            help="fasta reference genome sequence ", metavar="PATH", required=False)

    parser.add_argument("-pos", "--position",  dest="position",
            help="Positions file [hg19.txt or hg38.txt]", metavar="PATH", required=True)

    parser.add_argument("-out", "--output",
            dest="Outputfile", required=True,
            help="Folder name containing outputs", metavar="STRING")

    parser.add_argument("-r", "--Reads_thresh",
            help="The minimum number of reads for each base",
            type=int, required=False,
            default=50)

    parser.add_argument("-q", "--Quality_thresh",
            help="Minimum quality for each read, integer between 10 and 39, inclusive \n [10-40]",
            type=int, required=True)

    parser.add_argument("-b", "--Base_majority",
            help="The minimum percentage of a base result for acceptance \n [50-99]",
            type=int, required=True)

    parser.add_argument("-t", "--Threads", dest="threads",
            help="Set number of additional threads to use during alignment BWA-MEM",
            type=int, default=2)

    # MB Jan2020 Added these two for portability when calling from WGS Extract
    parser.add_argument("-py", "--Python3", dest="python3", help="Path and command to invoke Python3",
            metavar="PATH", required=True)       # MB Jan2020     Windows

    parser.add_argument("-samt", "--Samtools", dest="samtools", help="Path and command to invoke samtools",
            metavar="PATH", required=True)       # MB Jan2020    Windows

    args = parser.parse_args()
    return args


def get_frequency_table(mpileup):
    bases = ["A","T","G","C","+","-"]
    frequency_table = {}

    for i in mpileup.values:
        fastadict = {"A":0,"T":0,"G":0,"C":0}
        sequence = i[9] #actual sequence
        sequence = sequence.upper()
        sequence = trimm_caret(sequence)
        sequence = sequence.replace("$", "")
        indel_pos = find_all_indels(sequence)
        ### Count number of indels
        indels = count_indels(sequence, indel_pos)
        fastadict.update(indels)
        fastadict["-"] += sequence.count("*")
        ### Trimm Indels
        trimm_sequence = trimm_indels(sequence, indel_pos)
        for seq in trimm_sequence:
            if seq in fastadict:
                fastadict[seq] +=1
        frequency_table.update({i[3]:list(fastadict.values())})
    df_frequency_table = pd.DataFrame.from_dict(frequency_table, orient='index')
    df_frequency_table.columns = bases
    return df_frequency_table


def find_all_indels(s):
    find_all = lambda c,s: [x for x in range(c.find(s), len(c)) if c[x] == s]
    list_pos = []
    for i in find_all(s,"-"):
        list_pos.append(i)
    for i in find_all(s,"+"):
        list_pos.append(i)
    return sorted(list_pos)


def count_indels(s, pos):
    dict_indel = {"+":0,"-":0}
    if pos == []:
        return dict_indel
    if len(pos) > 0:
        for i in range(0,len(pos)):
            try: # in case it is not a number but a base pair e.g. A
                dict_indel[s[pos[i]]] += int(s[pos[i]+1])
            except ValueError:
                dict_indel[s[pos[i]]] += 1
                continue
    return dict_indel


def trimm_indels(s, pos):
    ## Receives a sequence and trimms indels
    if pos == []:
        return s
    u_sequence = ""
    start =  pos[0]
    count = (start+1)
    try: # in case it is not a number but a base pair e.g. A
        end = count+int(s[count])+1
    except ValueError:
        end = start+1
    u_sequence = s[:start]
    if len(pos) > 1:
        for i in range(1,len(pos)):
            start = end
            u_sequence += s[start:pos[i]]
            start = pos[i]
            count = (start+1)
            try: # in case it is not a number but a base pair e.g. A
                end = count+int(s[count])+1
            except ValueError:
                end = start+1
            if pos[-1] == pos[i]:
                #print(s[end:])
                u_sequence += s[end:]
    else:
        u_sequence += s[end:]
    return u_sequence


def trimm_caret(s):
    find_all = lambda c,s: [x for x in range(c.find(s), len(c)) if c[x] == s]
    list_pos = []
    for i in find_all(s,"^"):
        list_pos.append(i)
    if list_pos == []:
        return s
    i = 0
    start = 0
    end = 0
    sequence = ""
    while i<len(s):
        if s[i] == "^":
            end = i
            sequence += (s[start:end])
            start = i+1
        elif i >= list_pos[-1]+1:
            sequence += (s[list_pos[-1]+1:])
            break
        i+=1
    return sequence


def chromosome_table(bam_file,bam_folder,file_name):
    output = bam_folder+os_slash+file_name+'.chr'	# MB Jan2020 Windows
    tmp_output = "tmp_bam.txt"

    with open(tmp_output, "w") as f:
        subprocess.run([args.samtools, "idxstats", bam_file], stdout=f)		# MB Jan2020 Windows
        df_chromosome = pd.read_table(tmp_output, header=None)
        total_reads = sum(df_chromosome[2])
        df_chromosome["perc"] = (df_chromosome[2]/total_reads)*100
        df_chromosome = df_chromosome.round(decimals=2)
        df_chromosome['perc'] = df_chromosome['perc'].astype(str) + '%'
        df_chromosome = df_chromosome.drop(columns=[1,3])
        df_chromosome.columns = ['chr','reads','perc']
        df_chromosome.to_csv(output, index=None, sep="\t")

    os.remove(tmp_output)							# MB Jan2020 Windows et al

    if 'Y' in df_chromosome["chr"].values:
        return "Y", total_reads
    elif 'chrY' in df_chromosome["chr"].values:
        return "chrY", total_reads


def check_if_folder(path,ext):
    list_files = []
    if os.path.isdir(path):
        dirpath = os.walk(path)
        for dirpath, dirnames, filenames in dirpath:
            for filename in [f for f in filenames if f.endswith(ext)]:
                files = os.path.join(dirpath, filename)
                list_files.append(files)
        return list_files
    else:
        return [path]

def get_folder_name(path_file):
    #print("pathfile: " +path_file)					# MB Jan2020 debug
    folder      = path_file.split("/")[-1]
    #print("bam file: " +folder)						# MB Jan2020 debug
    folder_name = os.path.splitext(folder)[0]
    #print("bam base name: " +folder_name)				# MB Jan2020 debug
    return folder_name

def create_tmp_dirs(folder):
    print ("Folder in create_tmp_dirs: " + folder)	# MB Jan2020 debug
    if os.path.isdir(folder):
        print("WARNING! Temporary "+folder+" already exists. Removing.")
        try:
            shutil.rmtree(folder)
        except:
            pass        # OS Exception; often because folder in use already if user left something open looking at file
        else:
            os.mkdir(folder)    # Only do mkdir() if removal was successful
    else:
        os.mkdir(folder)							# MB Jan2020 Windows
    return True

def extract_haplogroups(path_Markerfile, Reads_thresh, Base_majority,
                        path_Pileupfile, log_output, fmf_output, Outputfile):
    print("Extracting haplogroups...")
    Markerfile = pd.read_csv(path_Markerfile, header=None, sep="\t")
    Markerfile.columns = ["chr", "marker_name", "haplogroup", "pos", "mutation", "anc", "der"]
    Markerfile = Markerfile.drop_duplicates(subset='pos', keep='first', inplace=False)

    Pileupfile = pd.read_csv(path_Pileupfile, header=None, sep="\t", dtype = {0:str,1:int,2:str,3:int,4:str,5:str})
    Pileupfile.columns = ['chr', 'pos', 'refbase', 'reads', 'align', 'quality']

    log_output_list = []
    log_output_list.append("Total of reads: "+str(len(Pileupfile))) #total of reads

    intersect_pos = np.intersect1d(Pileupfile['pos'], Markerfile['pos'])
    Markerfile = Markerfile.loc[Markerfile['pos'].isin(intersect_pos)]
    Markerfile = Markerfile.sort_values(by=['pos'])
    Pileupfile = Pileupfile.loc[Pileupfile['pos'].isin(intersect_pos)]

    Pileupfile = Pileupfile.drop(['chr'], axis=1)
    df = pd.merge(Markerfile, Pileupfile, on='pos')

    del [[Pileupfile,Markerfile]]
    gc.collect()
    Pileupfile=pd.DataFrame()
    Markerfile=pd.DataFrame()

    log_output_list.append("Valid markers: "+str(len(df))) #valid markers provided

    index_belowzero = df[df["reads"] == 0].index
    df_belowzero = df[df.index.isin(index_belowzero)]
    df_belowzero = df_belowzero.drop(['refbase','align','quality'], axis=1)
    df_belowzero["called_perc"] = "NA"
    df_belowzero["called_base"] = "NA"
    df_belowzero["state"] = "NA"
    df_belowzero["Description"] = "Position with zero reads"

    df = df[~df.index.isin(index_belowzero)]

    df_freq_table = get_frequency_table(df)
    df_freq_table = df_freq_table.drop(['+','-'], axis=1)
    df = df.drop(['refbase','align','quality'], axis=1)

    list_col_indices = np.argmax(df_freq_table.values, axis=1)
    called_base = df_freq_table.columns[list_col_indices]
    total_count_bases = np.sum(df_freq_table.values, axis=1)
    max_count_bases = np.max(df_freq_table, axis=1)
    called_perc = round((max_count_bases/total_count_bases)*100,1)

    bool_anc = np.equal(np.array(called_base), df["anc"].values)
    bool_der = np.equal(np.array(called_base), df["der"].values)

    bool_list_anc = np.where(bool_anc,'A','D')
    bool_list_anc = bool_list_anc.astype('object')
    bool_list_der = np.where(bool_der,'D','A')
    bool_list_der = bool_list_der.astype('object')
    bool_list_state = np.equal(bool_list_anc, bool_list_der)

    df["called_perc"] = np.array(called_perc, dtype=int)
    df["called_base"] = called_base
    df["state"] = bool_list_anc
    df["bool_state"] = bool_list_state

    df_discordantgenotype = df[~bool_list_state]
    df_discordantgenotype = df_discordantgenotype.drop(["bool_state"], axis=1)
    df_discordantgenotype["state"] = "NA"
    df_discordantgenotype["Description"] = "Discordant genotype"

    columns_fmf = df_discordantgenotype.columns

    ## read threshold
    df_readsthreshold = df[df["reads"] < Reads_thresh]
    df_readsthreshold["Description"] = "Below read threshold"

    ## filter by base percentage
    df_basemajority = df[df["called_perc"] < Base_majority]
    df_basemajority["Description"] = "Below base majority"

    index_to_remove = np.concatenate((df_basemajority.index,df_readsthreshold.index,
                                      df_discordantgenotype.index), axis=None)

    df_fmf = pd.concat([df_belowzero,df_readsthreshold, df_basemajority, df_discordantgenotype], axis=0, sort=True)
    df_fmf = df_fmf[columns_fmf]

    df_out = df[~df.index.isin(index_to_remove)]
    df_out = df_out.drop(["bool_state"], axis=1)
    df_out = df_out.sort_values(by=['haplogroup'], ascending=True)

    log_output_list.append("Markers with zero reads: "+str(len(df_belowzero)))
    log_output_list.append("Markers below the read threshold {"+str(Reads_thresh)+"}: "+str(len(df_readsthreshold)))
    log_output_list.append("Markers below the base majority threshold {"+str(Base_majority)+"}: "+str(len(df_basemajority)))
    log_output_list.append("Markers with discordant genotype: "+str(len(df_discordantgenotype)))
    log_output_list.append("Markers without haplogroup information: "+str(len(df_fmf)))
    log_output_list.append("Markers with haplogroup information: "+str(len(df_out)))

    with open(log_output, "a") as log:
        for marker in log_output_list:
            log.write(marker)
            log.write("\n")

    del [[df_basemajority,df_belowzero, df_discordantgenotype, df_readsthreshold, df_freq_table, df]]
    gc.collect()
    df_basemajority=pd.DataFrame()
    df_belowzero=pd.DataFrame()
    df_discordantgenotype=pd.DataFrame()
    df_readsthreshold=pd.DataFrame()
    df_freq_table=pd.DataFrame()
    df = pd.DataFrame()

    df_out = df_out[["chr","pos","marker_name","haplogroup","mutation","anc","der","reads","called_perc","called_base","state"]]
    df_fmf.to_csv(fmf_output, sep="\t", index=False)
    df_out.to_csv(Outputfile, sep="\t", index=False)


def logo():
    print(r"""

           |
          /|\
         /\|/\
        \\\|///
         \\|//
          |||
          |||
          |||

        """)


def run_cmd(cmd,msg,shell=True):
    start_time = time.time()
#    if platform.system() == "Windows":
#        cmd = cmd.replace("/","\\")
    print(cmd)
    subprocess.call(cmd, shell=True)
#    try:
#        result = subprocess.run([cmd], shell=shell)
#    except subprocess.CalledProcessError:
#        if subprocess.CalledProcessError.stdout:
#            print(subprocess.CalledProcessError.stdout)
#        if subprocess.CalledProcessError.stderr:
#            print(subprocess.CalledProcessError.stderr)
#       print(f"--- FAILURE in {msg} ---")
#        exit(subprocess.CalledProcessError.returncode)
#    finally:
#        if result.stdout:
#            print(result.stdout)
    total_time = round(time.time()-start_time)
    print(f"--- {total_time:n} seconds to finish {msg} ---")


# REH 13Mar2020 Execute main code here only if run as a script; not imported
if __name__ == "__main__":

    whole_time = time.time()

    print("\tErasmus MC Department of Genetic Identification \n\n\tYleaf: software tool for human Y-chromosomal \n\tphylogenetic analysis and haplogroup inference v2.1\n")
    print("\t(updated for compatibility with WGS Extract)\n")
    logo()
    args = get_arguments()

    if not (args.Fastq or args.Bamfile):                                      # REH 11 Mar2020 Style
        print("--- Yleaf failed! No FASTQ or BAM file(s) specified... ---")   # REH 11 Mar2020 Style
        exit()
    elif (args.Fastq and args.Bamfile):
        print("--- Yleaf failed! Both FASTQ and BAM file(s) specified... ---")   # REH 11 Mar2020 Style
        exit()

    if args.Fastq and args.reference is None:
        print("--- Yleaf failed! Missing -f reference model for FASTQ processing... ---")   # REH 11 Mar2020 Style
        exit()

    app_folder = os.path.dirname(os.path.realpath(__file__)).replace("\\","/")  # _FP
    app_folder += ("/" if app_folder[-1] != "/" else "")                        # Assure trailing slash

    out_file   = args.Outputfile                                                # _FN
    out_path   = os.path.split(out_file)[0]                                     # _oFP

    cwd        = os.getcwd()                                                    # _oFP
    if os.path.isabs(out_path) or cwd == os_slash:  # MB Jan2020 Windows; style
        out_folder = out_path
    else:
        out_folder = cwd + os_slash + out_path	    # MB Jan2020 Windows
    out_folder += ("/" if out_folder[-1] != "/" else "")                        # Assure trailing slash
    out_folder = out_folder.replace("\\","/")                                   # _FP

    print("App folder: " + app_folder)				# MB Jan2020 debug  REH10Mar2020 uncommented
    print("Out file: " + out_file)					# MB Jan2020 debug  REH10Mar2020 uncommented
    print("Cwd: " + cwd)							# MB Jan2020 debug  REH10Mar2020 uncommented
    print("Out folder: " + out_folder)				# MB Jan2020 debug  REH10Mar2020 uncommented

    if not create_tmp_dirs(out_folder):             # REH 13Mar2020 Style; side-effect
        print("--- Yleaf failed! Cannot create / open Output folder... ---")
        exit()

    print("Starting Yleaf ...")                     # REH 13Mar2020 Style (moved from within conditionals)
    # Algorithm as written assumes all FASTQ and BAM files in tree have different names
    # Algorithm as written cannot handle paired-end FASTQ files
    # Algorithm as originally written only allows FASTQ or BAM processing; not both
    type = 'FastQ' if args.Fastq else ('BAM' if args.Bamfile else '')       # REH 12Mar2020 guaranteed as checked above
    files = check_if_folder(args.Fastq,'.fastq') if args.Fastq else (check_if_folder(args.Bamfile,'.bam') if args.Bamfile else "")

    # Old folders for samtools were samtools-mingw and samtools-cygwin; so need to undo replace of folder name
    samtools = args.samtools
    bcftools = args.samtools.replace("samtools", "bcftools").replace("bcftools-","samtools-")
    bwa      = args.samtools.replace("samtools", "bwa").replace("bwa-","samtools-")

    for file in files:
        folder_name = get_folder_name (file)
        folder = os.path.join(app_folder,out_folder,folder_name,"")    # REH 11Mar2020 ugly side effect; force trailing slash
        folder = folder.replace("\\","/")       # Force into Universal OS slash form
        if not create_tmp_dirs(folder):         # REH 11 Mar2020 ugly side effect
            print(f"--- Yleaf failed! Cannot create Output folder {folder} for {type} file: {file}. Skipping... ---")
            continue

        bam_file_sorted = ""
        if args.Fastq:
            fastq_file = file
            sam_file = folder + folder_name + ".sam"       # MB Jan2020 Windows
            cmd = f"{bwa} mem -t {args.threads} {args.reference} {fastq_file} > {sam_file}"
            run_cmd(cmd, 'aligning reads to reference')

            bam_file = folder + folder_name + ".bam"        # MB Jan2020 Windows
            cmd = f"{samtools} view -@ {args.threads} -bS {sam_file} | samtools sort -@ {args.threads} -m 2G -o {bam_file}"
            run_cmd(cmd, 'converting SAM to BAM')

            os.remove(sam_file)                         # REH 13Mar2020 Windows

            cmd = f"{samtools} index -@ {args.threads} {bam_file}"
            run_cmd(cmd, 'indexing BAM file')
        else: # args.Bamfile:    # Must be args.bamfile (guaranteed)
            bam_file = file
            if not os.path.exists(bam_file + '.bai'):   # REH 13Mar2020 added back in but guaranteed not needed in WGSExtract
                bam_file_sorted = folder + folder_name + ".sorted.bam"
                cmd = f"{samtools} sort -m 2G -@ {args.threads} {bam_file} > {bam_file_sorted}"
                run_cmd(cmd, 'sorting BAM file')

                cmd = f"{samtools} index -@ {args.threads} {bam_file_sorted}"
                run_cmd(cmd, "indexing BAM file")
                bam_file = bam_file_sorted

        Outputfile = folder + folder_name + ".out"     # MB Jan2020 Windows
        log_output = folder + folder_name + ".log"		# MB Jan2020 Windows
        fmf_output = folder + folder_name + ".fmf"		# MB Jan2020 Windows
        pileupfile = folder + folder_name + ".pu" 		# MB Jan2020 Windows

        header,total_reads = chromosome_table(bam_file,folder,folder_name)
        # "samtools mpileup" changes to "bcftools mpileup --no-reference"; but formal VCF and not TSV file then also
        cmd = f"{samtools} mpileup -AQ{args.Quality_thresh} -r {header} {bam_file} > {pileupfile}"  # MB Jan2020
        run_cmd(cmd, 'executing Pileup')

        if os.path.exists(bam_file_sorted):
            os.remove(bam_file_sorted)
            os.remove(bam_file_sorted+".bai")

        start_time = time.time()
        extract_haplogroups(args.position, args.Reads_thresh, args.Base_majority,
                                    pileupfile, log_output, fmf_output, Outputfile)
#        os.remove(pileupfile)
        print("--- %.2f seconds in extracting haplogroups --- " % (time.time() - start_time) )

        # REH 13Mar2020 Ugly; need to replace with import of the file; not sub-process call
        cmd = f"{args.python3} {app_folder}predict_haplogroup.py -input {out_folder} -out {out_file}" # MB Jan2020 Windows
        run_cmd(cmd, "identifying haplogroup")

    print("--- %.2f seconds to complete Yleaf  ---" % (time.time() - whole_time))
