# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
 General Settings module support for the WGS Extract system

  This module is imported as wgse everywhere.  For convenience, we additionally use a local name to reduce name length
  in one instance (font).  Therefore, the code to use this module in other modules is:
    import settings as wgse
    font = wgse.font

  Because Global variables and constants are defined here, this file is included by all others.
  We must simply import; variable names cannot be selectively included with a "from" in Python.
  Imports of any needed System library functions here are buried in function definitions. (None
  are called often so there is not much penalty for doing so.)

  As Python has no structures, and global variable values are read when the file is imported, and redirected naming
   can get unwieldly, we do not use a Class that we later instantiate. It adds yet another layer of naming.
   So we insdtead define a single function init() here that is similar to a Class __init__ and that is expected
   to "instantiate" this "Class" and its global variables.  Then all the variables are simply available. This is all
   simply to reduce an extra level of naming.  All variables can simple be gotten by wgse.name where name is the
   variable name.

  See description at end of this file for naming conventions on variables of files / paths.

 Class(es): none
 Globals: see long list of primary level names below
 Function(s): init(), load(), save()
"""

# -----------------------------------------------------------------------------------------------------------------
# Program global constants
#
# There are no real constants in Python. But we treat these as such.  Read-only, set here once.

__version__ = "Beta ver3 (15 Jun 2021)"
manual_url = "https://wgsextract.github.io/"
# manual_url = "https://bit.ly/35IziTY"       # Version 3 Alpha manual (not publicly released yet)
ftdna_pubytree_url  = "https://www.familytreedna.com/public/y-dna-haplotree/"
yfull_searchsnp_url = "https://yfull.com/search-snp-in-tree/"
isogg_tree_url      = "https://isogg.org/tree/"

# Key window sizes (fixed size to force word wrap; target min is macbook pro 1366x768)
mainwindow_size   = '700x750'
microarr_winsize  = '615x750'
bamstats_winsize  = '1030x768'  # Need every vertical pixel allowed ...
simResult_winsize = '700x700'
simResult_maxw    = 700
simResult_maxh    = 700
yHgResult_winsize = '850x750'
yHgResult_maxw    = 850
yHgResult_maxh    = 750

font = {  # Boy do we have a hand-coded window GUI.  This at least gets us a little more generality ...
    '12':  ("Times New Roman", 12),         # For version/date/time in Title and various result pages (below min size)
    '13':  ("Times New Roman", 13),         # Stats by-chr table; Title manual / exit buttons (below min size)
    '14':  ("Times New Roman", 14),         # Minimum size per Marko (for laptops / Macbook with 1366x768 screen)
    '14b': ("Times New Roman", 14, "bold"),
    '14o': ("Times New Roman", 14, "overstrike"),
    '14u': ("Times New Roman", 14, "underline"),
    '16':  ("Times New Roman", 16),
    '16b': ("Times New Roman", 16, "bold"),
    '18':  ("Times New Roman", 18),
    '18b': ("Times New Roman", 18, "bold"),
    '20b': ("Times New Roman", 20, "bold"),
    '28b': ("Times New Roman", 28, "bold"),
    '32':  ("Arial Black", 32)              # Main Window Program Title
}

# Expected run times of various commands in seconds; for Please Wait window (used in module commandprocessor)
# Number below based on Randy's 40x, 57GB BAM file on his 2 core, AMD A10-5700 processor using CygWin htslib 1.10
expected_time = {    # Minutes * 60 Seconds
    'GetBAMHeader':               2,  # ## samtools view -H (less than a second usually)
    'ButtonAlignBAM':   8 * 60 * 60,  # ## bwa align of FASTQ to BAM
    'ButtonBAMStats':            15,  # ## samtools idxstats (2 secs if index file there; sometimes a little longer)
    'ButtonBAMStats2':       3 * 60,  # ## Nanopore long read; read more lines to get read length
    'ButtonBAMNoIndex':     35 * 60,  # ## samtools idxstats (35 min if file not indexed on 30x WGS)
    'ButtonBAMNoSort':      70 * 60,  # ## samtools idxstats (70+ min if file not sorted nor indexed on 30x WGS)
    'GenSortedBAM':         30 * 60,  # ## samtools sort
    'GenBAMIndex':          30 * 60,  # ## samtools index
    'ButtonMicroarrayDNA':   5 * 60,  # ## 1 to 12 aconv() calls; quick
    'ButtonCombinedKit':    50 * 60,  # ## Samtools mpileup, bcftools call on WGS; an hour or more for single processor
    'ButtonYHaplo':         10 * 60,  # ## yleaf haplo
    'ButtonYHaplo2':         5 * 60,  # ## yleaf lookup
    'ButtonMTHaplo':        20 * 60,  # ## haplogrep (java)
    'ButtonUnalignBAM':     75 * 60,  # ## samtools fastq on name sorted bam
    'ButtonUnmappedReads':  45 * 60,  # ## samtools view
    'ButtonMitoFASTA':       3 * 60,  # ## ButtonMitoVCF; samtools fasta
    'ButtonMitoBAM':         3 * 60,  # ## samtools view
    'ButtonMitoVCF':         3 * 60,  # ## bcftools mpileup, call
    'ButtonYandMT':          3 * 60,  # ## samtools view
    'ButtonYonly':           3 * 60,  # ## samtools view
    'CRAMtoBAM':            90 * 60,  # ## samtools cram to bam (1hr); samtools index (30min)
    'BAMtoCRAM':            30 * 60,  # ## samtools bam to cram (1hr); samtools index (30min)
    'CoverageStats':        45 * 60,  # ## samtools coverage to get per chromosome coverage
    'CoverageStatsWES':     60 * 60,  # ## samtools depth w/ WES BED file to get Coverage and Avg Read Depth
    'CreateAlignIndices':   60 * 60,  # ## bwa index on fasta file
    'AlignCleanup':        115 * 60,  # ## Fixmate, Sort, Markdup, Index
    'LiftoverCleanup':            5,  # ## Sort and Compress of CombinedKit file
    'AnnotatedVCF-yOnly':   10 * 60,  # ## Extract Y-only VCF from BAM and annotate
    'UnsortBAM':            10 * 60   # ## samtools reheader (to change coord sorted to unknown) (DEBUG_MODE only)
}

# Names used in calls to commandprocessor "run_bash_script"; called from modules mainwindow, bamfiles, microarray,
#     and hg38tohg19; keys also defined as keys in the langstrs i8n dictionary keys (languages.xlsx)
#   Times are adjusted relative to a BAM File Size of 60 GB.  So partial BAMs may get weighted down significantly.

# Used in module bamfiles (specifically, calc_stats()
valid_autos = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
               "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]
valid_other = ["M", "X", "Y"]  # REH 20Mar2020 valid split to sort autosomes numerically; REH 25Feb2021 removed '*'

# Used exclusively in hg38tohg19 Liftover module (microarray); just brought here for clarity
valid_chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                     "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                     "chr20", "chr21", "chr22", "chrM", "chrX", "chrY"]

# Currently, we just fork off 25 jobs as if 32 processors.  Can we more smartly manage when less processors?
# This is a load balance grouping based on the max available processors power of two groups
parallelBAM = {
    32: ["chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8", "chr9",
         "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
         "chr20", "chr21", "chr22", "chrM", "chrX", "chrY"],
    16: ["chr1", "chr2", "chr3,chr22", "chr4,chr21", "chr5,chr20", "chr6,chr19", "chr7,chr18", "chr8,chr17",
         "chr9,chr16", "chr10,chr15", "chr11,chr14", "chr12,chr13", "chrM,chrX,chrY"],
     8: ["chr1,chrM,chrX,chrY", "chr2,chr8,chr17", "chr3,chr9,chr16,chr22", "chr4,chr6,chr19,chr21",
         "chr5,chr12,chr13,chr20", "chr7,chr10,chr15,chr18", "chr11,chr14"],
     4: ["chr1,chr11,chr14,chrM,chrX,chrY", "chr2,chr7,chr8,chr10,chr15,chr17,chr18",
         "chr3,chr5,chr9,chr12,chr13,chr16,chr20,chr22", "chr4,chr6,chr19,chr21"],
     2: ["chr1,chr2,chr7,chr8,chr10,chr11,chr14,chr15,chr17,chr18,chrM,chrX,chrY",
         "chr3,chr5,chr9,chr12,chr13,chr16,chr20,chr22,chr4,chr6,chr19,chr21"],
     1: ["."]
}

# Same order as Valid Chromosomes; numbers by James Kane from GATK CallableLoci
# https://docs.google.com/spreadsheets/d/1S3a69mxHeiwRb2l3s_uyFK5KGjHhL9Qv1odObJBgOdQ/view#gid=477598766
# Note: override M which actually has 1. Causes breadth of coverage to become 100.01%
nadjust = {
    38: {"1": 18475410, "2": 1645301, "3": 195424, "4": 461888, "5": 2555066, "6": 727457, "7": 375842, "8": 370500,
         "9": 16604167, "10": 534460, "11": 552880, "12": 137493, "13": 16381203, "14": 18525617, "15": 17349864,
         "16": 8532402, "17": 337237, "18": 283680, "19": 2851219, "20": 499910, "21": 8671412, "22": 13746488,
         "X": 1147866, "Y": 33591060, "M": 0},
    37: {"1": 23970000, "2": 4994855, "3": 3225295, "4": 3492600, "5": 3220000, "6": 3720001, "7": 3785000,
         "8": 3475100, "9": 21070000, "10": 4220009, "11": 3877000, "12": 3370502, "13": 19580000, "14": 19060000,
         "15": 20836626, "16": 11470000, "17": 3400000, "18": 3420019, "19": 3320000, "20": 3520000, "21": 13023253,
         "22": 16410021, "X": 4170000, "Y": 36389037, "M": 0},
    19: {"1": 23970000, "2": 4994855, "3": 3225295, "4": 3492600, "5": 3220000, "6": 3720001, "7": 3785000,
         "8": 3475100, "9": 21070000, "10": 4220009, "11": 3877000, "12": 3370502, "13": 19580000, "14": 19060000,
         "15": 20836626, "16": 11470000, "17": 3400000, "18": 3420019, "19": 3320000, "20": 3520000, "21": 13023253,
         "22": 16410021, "X": 4170000, "Y": 36389037, "M": 0}
}

''' Variables contain File Name Path, Base, and Suffix (FPBS) naming conventions explained at the end of this file. '''
# --------------------------------------------------------------------------------------------------------------------
# Program global variables
#

# Variables initially set here but used globally in the WGSE system. Accessed via wgse.<varname>

DEBUG_MODE    = None  # Global DEBUG mode setting (similar to __debug__ removed by python -O)
wsl_bwa_patch = None  # To bypass Win10 BWA which is single processorlanguage

# Class object instantiation points
outdir = None  # Output Directory subsystem
reflib = None  # Reference Library subsystem
tempf  = None  # TemporaryFile subsystem (user override of default location possible)
lang   = None  # Internationalization Language subsystem (fixed location)
window = None  # main Tk() window (fixed once created unless language change asked)
BAM    = None  # main BAM Class instance with all its stats, attributes, etc


# OS / HW Platform specific variables
os_plat     = None  # type: [str]  # OS Platform; to avoid multiple calls to system.platform
os_slash    = None  # type: [str]  # OS Platform file path slash (forward or backward)
os_batch_FS = None  # type: [str]  # OS Platform batch file suffix (.bat or .sh)
os_threads  = 1     #              # OS Platform number of virtual processors; will set later but default is 1 until set
os_pid      = 1000  #              # OS Process ID of this run; just set to 1000 till read to get the type

# from Typing import Optional
# Users home area (determined from OS) and setting files we will look for and manipulate there
User_oFP     = None  # Users home area directory path (~, users/name)
debugset_oFN = None  # type: [str]  # File name for DEBUG mode toggle
wgseset_oFN  = None  # type: [str]  # File name for global settings
wslbwa_oFN   = None  # type: [str]  # File name for WSL BWA Patch toggle

# Key paths all determined from where this settings file is located.
install_FP   = None  # type: [str]  # Installation path for WGS Extract; program/ below is where this file is
install_oFP  = None  # type: [str]  # Installation path for WGS Extract; program/ below is where this file is
prog_oFP     = None  # type: [str]  # File path of this Python file (location of WGS Extract program)
prog_FP      = None  # type: [str]  # File path of this Python file (location of WGS Extract program)
language_oFN = None  # type: [str]  # File path and name of Internationalization language settings
image_oFP    = None  # type: [str]  # Icon image in program header file pointer
dnaImage     = None  #              # Placeholder for banner image (used in MainWindow; need gloval variable)
icon_oFP     = None  # type: [str]  # For Windows Icon image file pointer

# Key external program path, names and suffixes (OS dependent)
python3_FP    = None  # type: [str]  # Installation path to Python3 installation
python3x_FN   = None  # type: [str]  # Actual file name w/ path of Python3 executable

yleaf_FP       = None  # type: [str]  # Installation path for our version of yLeaf
microarray_FP  = None  # type: [str]  # Installation path for our microarray extract tool (currently in wgesextract)
microarray_oFP = None  # type: [str]  # Installation path for our microarray extract tool (currently in wgesextract)
# liftover_chain_oFN = None  # Fully qualified name of liftover chain file in reference_library

javax_FN      = None  # type: [str]  # Java runtime path executable file name
javargs       = None  # type: [str]  # "-Xmx2g --module-path=xxxx -jar"
javax_FNp     = None  # type: [str]  # Java runtime executable with standard parameters (javaargs)

jartools_FP   = None  # type: [str]  # Installation path for jars
haplogrepx_FN = None  # type: [str]  # Haplogrep command (Java Jar file)
picardx_FN    = None  # type: [str]  # Broad Institute Picard (part of GATK3)
gatk3x_FN     = None  # type: [str]  # Broad Institute GATK3
gatk4x_FN     = None  # type: [str]  # Broad Institute GATK4 (maybe setup subfolder for release)
igvx_FN       = None  # type: [str]  # IGV Genomics Viewer (maybe setup subfolder for release)

bashx_oFN     = None  # type: [str]  # Bash shell (needed in Windows only; to run a shell)
headx_qFN     = None  # type: [str]  # Unix head command
awkx_qFN      = None  # type: [str]  # Unix awk command
grepx_qFN     = None  # type: [str]  # Unix grep command
sortx_qFN     = None  # type: [str]  # Unix sort command
catx_qFN      = None  # type: [str]  # Unix cat/type command
sedx_qFN      = None  # type: [str]  # Unix sed command
zipx_qFN      = None  # type: [str]  # Unix zip command
unzipx_qFN    = None  # type: [str]  # Unix unzip command

samtools_oFP  = None  # type: [str]  # Installation path to HTSLib Bioinformatic tools (samtools, bcftools, etc)
samtools_FP   = None  # type: [str]  # Installation path to HTSLib Bioinformatic tools (samtools, bcftools, etc)
samtoolsx_qFN = None  # type: [str]  # HTSlib samtools bioinformatics command (WGSE v3, samtools 1.10+)
bcftoolsx_qFN = None  # type: [str]  # HTSLib bcftools bioinformatics command
tabixx_qFN    = None  # type: [str]  # HTSLib tabix bioinformatics command
bgzipx_qFN    = None  # type: [str]  # HTSlib bgzip bioinformatics command
bwax_qFN      = None  # type: [str]  # BWA MEM Alignment command
bwamem2x_qFN  = None  # type: [str]  # BWA MEM v2 Alignment command
fastqcx_FN    = None  # type: [str]  # Fastqc analysis tool


# Todo make table of all external execs needed (i.e. a registry) with name, status if available, version, etc.
#   fill table at start; then check JIT when needed for availability (no error till needed).  That way,
#   can simply return with NOOP if necessary programs are not available.  Make sure to account for PIP/PythonLib
#   and Java Jar availability also; not just executables?  Confirm actually used in scripts here (or likely needed)
#   AND in the install environment

# -----------------------------------------------------------------------------------------------------------------
# Program global routine(s)
#

def init(top_level=False):
    """
    Main routine, like a Class init, to setup all the global variables. We do local imports here as this module / file
    is imported everywhere.
    """
    import os  # os.path, etc
    import locale  # For setting locale for number format (thousands separator, etc)
    import platform  # platform.system
    # import multiprocessing

    # Local imports as this module / file is imported everywhere else; no .h/,c separation like in C
    from utilities import DEBUG, is_legal_path, nativeOS, universalOS, wgse_message
    from utilities import TemporaryFiles, LanguageStrings, OutputDirectory       # subsystem classes
    from commandprocessor import is_command_available
    from referencelibrary import ReferenceLibrary           # subsystem class
    from mainwindow import init_mainWindow

    # Globals we want to access from in here
    global DEBUG_MODE, wsl_bwa_patch                        # Some universal settings
    global tempf, lang, outdir, reflib, window, BAM         # Some universal classes
    global os_plat, os_slash, os_batch_FS, os_threads, os_pid
    global User_oFP, debugset_oFN, wgseset_oFN, wslbwa_oFN  # , langset_oFN
    global prog_oFP, prog_FP, language_oFN, image_oFP, dnaImage, icon_oFP
    global install_FP, install_oFP
    global python3_FP, python3x_FN, yleaf_FP, microarray_FP, microarray_oFP     # , liftover_chain_oFN
    global javax_FN, javargs, javax_FNp, jartools_FP, haplogrepx_FN, picardx_FN, gatk3x_FN, gatk4x_FN, igvx_FN
    global bashx_oFN, headx_qFN, awkx_qFN, grepx_qFN, sortx_qFN, catx_qFN, sedx_qFN, zipx_qFN, unzipx_qFN
    global samtools_oFP, samtools_FP, samtoolsx_qFN, bcftoolsx_qFN, tabixx_qFN, bgzipx_qFN, bwax_qFN, bwamem2x_qFN
    global fastqcx_FN

    # OS / current host platform Specific values
    os_plat = platform.system()     # Windows, Darwin, ...
    os_slash, os_batch_FS = ('\\', '.bat') if os_plat == "Windows" else ('/', '.sh')
    os_threads = os.cpu_count()     # For parallelizing commands that can use multiple processors
    os_pid = os.getpid()            # Will use to create a unique area in the Temp directory for our run

    #
    # Stored settings set per user / run
    #

    # Users home directory for settings (Unix hidden style with dot (.))
    User_oFP = os.path.expanduser("~")
    if User_oFP[-1] != os_slash:  # In case returns as root (/) or (C:\\)
        User_oFP += os_slash  # Assure trailing os_slash
    debugset_oFN  = f'{User_oFP}.wgsedebug'   # No content; just existence of file turns on debugging
    wslbwaset_oFN = f'{User_oFP}.wgsewslbwa'  # No content; just existence turns on WSL BWA patch on Win10 systems
    wgseset_oFN   = f'{User_oFP}.wgsextract'  # General settings save / restore

    # Start global debug messages if requested (utilities.py); start after TemporaryFiles so it can clean directory
    if os.path.exists(debugset_oFN) and os.path.isfile(debugset_oFN):
        DEBUG_MODE = True
        DEBUG("***** Debug Mode Turned On *****") if top_level else ''

    # Special for Win10 patch to use WSL BWA as CygWIn64 is not running multiproc
    if os.path.exists(wslbwaset_oFN) and os.path.isfile(wslbwaset_oFN):
        wsl_bwa_patch = True
        DEBUG("***** WSL BWA Patch Mode Turned On *****") if top_level else ''

    # Need to set default locale; for decimal number (:n) formatting;
    #   override to POSIX when needed for consistent sort order (BAM header signature)
    # On Unix-based systems, sets language in OS messagebox dialogs also (Win10 does not allow altering)
    DEBUG(f'Locale: {locale.getlocale()[0]}')
    DEBUG(f'Default Locale: {locale.getdefaultlocale()[0]}')
    locale.setlocale(locale.LC_ALL, '')

    #
    # Global, hardwired locations of key files, areas based on where this file is located
    #

    # Setup all the needed paths and filenames; based mostly on where this program (Python script) resides
    prog_oFP = os.path.dirname(os.path.abspath(__file__))       # Where does this Python program reside?
    # Todo should call is_legal_path() here at first use; but language system not setup yet to report error
    if prog_oFP[-1] != os_slash[0]:  # In case returns as root (/) or (C:\\); never should as buried in installation ...
        prog_oFP += os_slash  # Assure trailing os_slash
    prog_FP = universalOS(prog_oFP)

    # Some key static files we need to read in and use
    image_oFP = f"{prog_oFP}img{os_slash}dna.png"
    icon_oFP  = f"{prog_oFP}img{os_slash}favicon.ico"
    language_oFN = f'{prog_oFP}language.xlsx'  # Moved from subdir of files in v2 to single file here; .txt to .csv

    # Find program installation by assuming this file is in program/ below the installation root
    install_FP  = prog_FP.split('program')[0]  # Has trailing os_slash
    install_oFP = prog_oFP.split('program')[0]

    tempfiles_oFN = f'{install_FP}temp/'        # Default installation location
    reflib_oFN = f'{install_FP}reference/'      # Default installation location

    #
    # Setup bioinformatic, OS and other tools we need to access too
    #

    # Three guaranteed folders in install: program (prog_FP), yleaf, jartools and temp
    yleaf_FP       = f'{install_FP}yleaf/'
    jartools_FP    = f'{install_FP}jartools/'
    microarray_FP  = f'{prog_FP}microarray/'     # Note: below the program/ folder; not from install root
    microarray_oFP = nativeOS(microarray_FP)

    if is_command_available("java", "--version", True, 0):
        javax_FN = "java"       # Most installations, if installed, on the call path and so just mentioned
    elif is_command_available(f'{install_FP}jre/bin/java', "--version", True, 0):
        javax_FN = f'{install_FP}jre/bin/java'      # If Win10 installer had to install, then local to WGSE
    else:
        javax_FN = None
    # javargs = f'-Xmx2g --module-path={jartools_FP} -jar'
    javargs = f'-Xmx2g -jar'        # Not sophisticated enough yet to have modules in jartools
    javax_FNp = f'{javax_FN} {javargs}'

    # Set all the OS platform specific path / program names; see file variable naming spec at end; 'x' is executable
    if os_plat == "Windows":  # .exe extension works in .sh and requried in .bat commands for Windows
        python3_FP  = f'{install_FP}python/'
        python3x_FN = f'{python3_FP}python.exe'

        samtools_oFP = f'{install_oFP}win10tools\\bin\\'
        samtools_FP  = f'{install_FP}win10tools/bin/'

        bashx_oFN = f'{samtools_oFP}bash.exe'       # Exception, native-OS and no quote
        headx_qFN = f'"{samtools_FP}head.exe"'
        awkx_qFN  = f'"{samtools_FP}gawk.exe"'
        grepx_qFN = f'"{samtools_FP}grep.exe"'
        sortx_qFN = f'"{samtools_FP}sort.exe"'      # Only in hg39tohg19 liftover module
        catx_qFN  = f'"{samtools_FP}cat.exe"'       # Only in hg39tohg19 liftover module
        sedx_qFN  = f'"{samtools_FP}sed.exe"'       # Only in microarray
        zipx_qFN  = f'"{samtools_FP}zip.exe"'       # Only in microarray
        unzipx_qFN = f'"{samtools_FP}unzip.exe"'     # Only in microarray

        samtoolsx_qFN = f'"{samtools_FP}samtools.exe"'
        bcftoolsx_qFN = f'"{samtools_FP}bcftools.exe"'
        tabixx_qFN    = f'"{samtools_FP}tabix.exe"'
        bgzipx_qFN    = f'"{samtools_FP}bgzip.exe"'  # Only in export_unmapped_reads
        bwax_qFN      = f'"{samtools_FP}bwa.exe"'
        bwamem2x_qFN  = f'"{samtools_FP}bwamem2.exe"'
    elif os_plat == "Darwin":  # Apple MacOS
        python3_FP   = '/usr/local/bin/'
        python3x_FN  = f'{python3_FP}python3'

        bashx_oFN = ''  # only used in Windows
        headx_qFN = '"/usr/bin/head"'
        awkx_qFN  = '"/usr/bin/awk"'
        grepx_qFN = '"/usr/bin/grep"'
        sortx_qFN = '"/usr/bin/sort"'
        catx_qFN  = '"/bin/cat"'
        sedx_qFN  = '"sed"'
        zipx_qFN  = '"zip"'
        unzipx_qFN = '"unzip"'

        samtools_oFP  = '/opt/local/bin/'
        samtools_FP   = '/opt/local/bin/'
        samtoolsx_qFN = f'"{samtools_FP}samtools"'
        bcftoolsx_qFN = f'"{samtools_FP}bcftools"'
        tabixx_qFN    = f'"{samtools_FP}tabix"'
        bgzipx_qFN    = f'"{samtools_FP}bgzip"'  # Only used in export_unmapped_reads at the current time
        bwax_qFN      = f'"{samtools_FP}bwa"'
        bwamem2x_qFN  = f'"{samtools_FP}bwamem2"'
    else:  # All others (Unix, Linux); most simply found in path
        python3_FP   = ''
        python3x_FN  = 'python3'

        bashx_oFN = ''  # only used in Windows
        headx_qFN = '"head"'
        awkx_qFN  = '"awk"'
        grepx_qFN = '"grep"'
        sortx_qFN = '"sort"'
        catx_qFN  = '"cat"'
        sedx_qFN  = '"sed"'
        zipx_qFN  = '"zip"'
        unzipx_qFN = '"unzip"'

        samtools_oFP  = ''
        samtools_FP   = ''
        samtoolsx_qFN = '"samtools"'
        bcftoolsx_qFN = '"bcftools"'
        tabixx_qFN    = '"tabix"'
        bgzipx_qFN    = '"bgzip"'  # Only used in export_unmapped_reads at the current time
        bwax_qFN      = '"bwa"'
        bwamem2x_qFN  = '"bwamem2"'

    #
    # Start up key subsystems:  Language Translation, Temporary Files directory, Reference Library directory
    #

    # window not null says init_mainWindow has been run, dnaImage not null says setup_mainWindow has been run
    window = init_mainWindow() if top_level else None   # pseudo __init__ call for non-class; creates root window

    # Start i18n Language subsystem (utilities.py)  (start first so warnings / error messages can be translated)
    lang = LanguageStrings(language_oFN)  # Start language subsystem (utilities.py)
    # lang.read_language_setting(langset_oFN)  # Try stored lang setting; otherwise ask user if missing / error

    # Start Temporary Files subsystem (utilities.py) with default location
    tempf = TemporaryFiles(tempfiles_oFN, top_level)  # Initiate TemporaryFiles @ default, do clean of directory

    # Start Reference Library subsystem (referencelibrary.py) with default location
    reflib = ReferenceLibrary(reflib_oFN)  # Initiate ReferenceLibrary @ default location

    # Start Output Directory subsystem (utilities.py) with none (no default).
    outdir = OutputDirectory()      # No default value

    BAM = None        # No init for BAMFile class; just keep it unset for now

    load_settings(top_level)  # Grab default setting overrides from users directory; finish subsystem setup

    # TODO Is this still needed now that we are quoting all file paths and names?
    if not is_legal_path(nativeOS(prog_FP)):
        # No wgse.window yet; but need wgse.lang.i18n[] for the translated error message
        wgse_message("error", 'InvalidPathWindowTitle', True,
                     lang.i18n['InvalidPathErrorMessage'].replace("{{PATH}}", prog_FP))
        exit()


# noinspection PyUnresolvedReferences
def load_settings(top_level):
    """
    To read / load saved settings on program startup. For now, just a select few values restored.
    """
    import os
    import json
    from utilities import DEBUG, wgse_message
    from bamfiles import BAMFile, BAMContentError, BAMContentWarning

    global outdir, reflib, tempf, lang, BAM
    global wgseset_oFN

    settings_to_restore = {}
    if os.path.exists(wgseset_oFN):
        try:
            with open(wgseset_oFN, "r") as f:
                settings_to_restore = json.load(f)
            DEBUG(f'Settings to restore: {settings_to_restore}')
        except:  # Exception here if file processing error on existing file; so delete
            DEBUG(f'*** Error processing settings file; deleting')
            os.remove(wgseset_oFN)
            return
        print("Processing saved settings file ...\n")

    # Restore value IF set in file and class subsystem setup to receive the value; else retain old value
    # We have avoided setting up each subsystem until this call as a stored setting may override default
    #  (or in the case of language, the stored setting or pop-up to user selects the language)

    if not lang:
        DEBUG("*** FATAL ERROR: Language subsystem not available!")
        raise
    lang.change_language(settings_to_restore.get('lang.language', lang.language))

    if not reflib:
        DEBUG("*** FATAL ERROR: Reference Library subsystem not available!")
        raise
    reflib.change(settings_to_restore.get('reflib.FP', reflib.FP))

    if not tempf:
        DEBUG("*** FATAL ERROR: Temporary Files subsystem not available!")
        raise
    tempf.change(settings_to_restore.get('tempf.FP', tempf.FP))

    if not (tempf.oFP and os.path.isdir(tempf.oFP)):
        # Cannot proceed with setting up BAM file from settings if there is no valid temporary files area
        # Either the default directory (temp/ in the installation) or a stored setting must be valid
        DEBUG(f'*** FATAL ERROR - no valid Temporary File Path: {tempf.oFP}')
        wgse_message("error", 'InvalidTempDirTitle', False, 'errTempDirPath')
        exit()

    if not outdir:
        DEBUG("*** FATAL ERROR: Output Directory subsystem not available!")
        raise
    outdir.change(settings_to_restore.get('outdir.FP', outdir.FP))

    if top_level:
        saved_BAM_FN = settings_to_restore.get('BAM.file_FN', "")
        if saved_BAM_FN:
            try:
                BAM = BAMFile(saved_BAM_FN)
            except BAMContentError as err:
                BAM = None      # Ignore errors from trying restore; simply ignore setting
            except BAMContentWarning as warn:
                pass
            else:
                pass    # Ignore issuing warnings about CRAM, unindexed or unsorted file
            if outdir.FP and BAM:  # Sets part of output File variables that includes BAM base name
                outdir.oFPB = outdir.oFP + BAM.file_FB
                outdir.FPB  = outdir.FP  + BAM.file_FB


def save_settings():
    """
    To store "saved" settings on program exit or whenever changed.  For now, just a select few values:
    DEBUG_MODE, language, output directory, reference library directory, temporary directory, wsl_bwa_patch
    """
    import json
    from utilities import DEBUG

    global outdir, BAM, reflib, tempf, lang
    global wgseset_oFN

    # DEBUG(f'WGSE JSON Setting file: {wgseset_oFN}')
    settings_to_save = {}

    # Only save value if set and not default (ones with default have "set" variable)
    if lang and lang.language:
        settings_to_save['lang.language'] = lang.language
    if reflib and reflib.set:
        settings_to_save['reflib.FP'] = reflib.FP
    if tempf  and tempf.set:
        settings_to_save['tempf.FP'] = tempf.FP
    if outdir and outdir.FP:
        settings_to_save['outdir.FP'] = outdir.FP
    if BAM and BAM.file_FN:
        settings_to_save['BAM.file_FN'] = BAM.file_FN
    # DEBUG(f'WGSE JSON Settings to save: {settings_to_save}')

    # Do not save if no values set; language should always beset though
    if settings_to_save:
        try:
            with open(wgseset_oFN, 'w') as f:
                json.dump(settings_to_save, f)
        except:
            DEBUG("*** Error writing settings file; skipping.")
            pass    # We try. If cannot write file then just move on silently


"""      ******************  FILE VARIABLE NOMENCLATURE  ***********************
 We do a lot of file name construction and manipulation. So we need clear conventions for naming the variables.
 File Variable Naming Acronyms (_suffix):
  FN    -- fully qualified File Name (path, base, suffix; i.e FPBS) (absolute or relative) (DEFAULT is uFN)
  FP    -- File Path name only (including trailing slash; head in split + /)
  FB    -- File Base name only (NO path, NO suffix) (simple root; splitext[0].split[1])
  FS    -- File Suffix only (including dot prefix; also called extension in splitext)
  FBS   -- File Base and Suffix (traditional file name; tail in split, basename())
  FPB   -- File Path and Base name (no suffix) (called root in splitext)
  qFN   -- Quoted fully qualified File Name (would never quote anything else) (with possible spaces in name)
  oFN   -- Native OS (slash) File Name (what you pass to open(); also oFP if directory)
  wFN   -- Windows (back slash) path File Name (or FP, etc)
  uFN   -- Unix/Universal (forward slash) path File Name (or FP, etc) (DEFAULT, if not specified)
            uFN is DEFAULT for ALL OS's in BASH and Win10 PowerShell; C:/users/name/.wgselang is OK
  x_qFN -- indicates actual executable program file name; usually always quoted but not in Native OS format

  In our terminology here:
   os.path.split elements:    [ oFP , FBS ]
   os.path.splitext elements: [ oFPB , FS ]
   os.path.basename returns FBS (same as split[1])
  But remember os.path is OS specific.  Pathlib is the OS independent representation.

  As most of the manipulation is done for BASH scripts; the default stored is Unix paths (forward slash)
  Only when needed in a native / root OS form for Windows, will it be converted (i.e. to Windows back slash)
  wgse.os_slash will contain the appropriate slash for file paths in the current OS
  outFPB is the historic outfws + '/' (fws was file withOUT suffix)
  
  Python "open" command needs oFN without quotes (native OS, uses appropriate slash).
  Bash scripts need qFN (in Unix slash format); even on Windows machines when path includes drive letter format.
  Python "process" command is a list of strings so no quoting of individual strings is necessary; not shell parsed.
"""

# Todo
#  Beginning setup here for a Check-Generate/Gather subsystem to allow button commands to simply specify what they need.
#  Allowing the Check-Generate/Gather to try and find AND make available the necessary files. Can maybe automatically
#  invoke major functions like alignment, variant calling, etc. if the button command does not care
#  about the parameters used to get from a to b. Also allows test vendor supplied files (or previously user generated
#  ones) to be looked for and utilized directly.  To be moved up into the main section once settled and implemented.
#  Thinking of using bit masks for efficiency over dictionaries.  But may need more parameters attached so maybe not?
'''
# Mask Categories for Required Files
mask_UserFiles
mask_RefGenome
mask_Liftover
mask_Microarray
mask_CombinedKit

# Masks for Required User Input Files
mask_FASTQ              = 0x0001
mask_FASTQ_PE           = 0x0002
mask_BAM                = 0x0010
mask_BAM_Sorted         = 0x0020
mask_BAM_Index          = 0x0040
mask_BAM_Unaligned      = 0x0080
mask_CRAM               = 0x0100
mask_CRAN_Index         = 0x0200
mask_VCF_RAW
mask_VCF_RAW_Index
mask_VCF_Called
mask_VCF_Called_Index
mask_CombinedKit

# Masks for Required Microarray Template files
mask_23andMe_API
mask_23andMe_v2
mask_23andMe_v3
mask_23andMe_v4
mask_23andMe_v5
mask_Ancestry_v1
mask_Ancestry_v2
mask_FTDNA_v1
mask_FTDNA_v2
mask_FTDNA_v3
mask_LivDNA_v1
mask_LivDNA_v2
mask_MyHeri_v1
mask_MyHeri_v2
mask_NGG_v2
mask_NGG_v2NG
'''