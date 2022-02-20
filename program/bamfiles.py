# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
    The first major module / subsystem that is a proper object Class.  BAMFile (class; bamfiles module) is the central
    handler for all things BAM.  Processing, determining stats, manipulating (sorting, indexing; eventually CRAM/SAM
    processing)  The idea is a class instance for each BAM file specified.  Eventually allowing for multiple to be
    known at the same time.  Once a BAM file is selected in the OS interface in the main Window module, the BAM
    object is instantiated and processing begins. Is the main driver of the referencegenome module.
"""

import os       # for path, stat
#from tkinter.ttk import Button, Label
from tkinter import Toplevel, Radiobutton, Button, Label, StringVar

import settings as wgse
font = wgse.font
from utilities import DEBUG, is_legal_path, nativeOS, universalOS, unquote, Error, Warning, wgse_message
from commandprocessor import run_bash_script


######################################################################################################################
# BAM File Processing Error Classes
class BAMContentError(Error):
    def __init__(self, reason):
        self.reason = wgse.lang.i18n[reason]

class BAMContentWarning(Warning):
    def __init__(self, reason):
        self.reason = wgse.lang.i18n[reason]

######################################################################################################################
# BAM File Processing class

class BAMFile:
    """
        BAM File class object; including all stats and processing thereof
        In init, simply assert if error and let outer caller worry about deleting the object
        (i.e. we assume within invoked from within a try...except). Errors will be of class BAMContentError
        with an additional message parameter. The message parameter must be a LanguageStrings dictionary key.
    """
    Refgenome = None
    askrefgenomewindow = None

    def __init__(self, BAM_FN):        # Formerly known as process_BAM() before becoming Class.__init__
        """ Check BAM File name / path is legal; fill in basic stats from BAM file if so. Assert on error. """

        BAM_oFN = nativeOS(BAM_FN)
        if not is_legal_path(BAM_oFN):
            raise BAMContentError('errBAMFileSpecialChars')
        if not os.path.isfile(BAM_oFN) or os.path.getsize(BAM_oFN) < 10000:    # file is bad
            raise BAMContentError('errBAMContent')

        # OK, committed to this new file.  Let's start processing
        self.set_default_BAM_file_settings()  # Just to make sure we start with a clean slate

        # set_BAM_file_settings ... just doing file names. Most content set later while processing file.
        self.file_oFN  = BAM_oFN
        self.file_FN   = BAM_FN
        self.file_qFN  = f'"{self.file_FN}"'

        self.file_FP   = os.path.dirname(BAM_FN)
        self.file_FP  += '/' if self.file_FP[-1] != '/' else ''  # Assure trailing slash; universalOS so forward slash

        self.file_oFPB = os.path.splitext(BAM_oFN)[0]
        self.file_FPB  = universalOS(self.file_oFPB)

        self.file_FBS  = os.path.basename(BAM_FN)
        self.file_FB, self.file_FS = os.path.splitext(self.file_FBS)
        self.disp_FBS  = self.file_FBS if len(self.file_FBS) < 31 else f'{self.file_FBS[:12]} ... {self.file_FBS[-13:]}'

        self.file_type = self.file_FS[1:].upper()    # strip leading dot and upcase
        if not self.file_type in ["BAM", "CRAM", "SAM"]:    # Internal call may pass an incorrect file
            raise BAMContentError('errBAMFileExtension')

        self.file_stats = os.stat(self.file_oFN)
        if self.file_stats.st_size == 0:
            raise BAMContentError('errBAMFileEmpty')
        self.relfsize = max(0.01, self.file_stats.st_size / (45 * 10.0 ** 9))    # Relative size to 45GB (time adjust)

        # Do a bunch of quick things, custom here, that lets us characterize the BAM a bit. May pop-out witha Raise ...
        wgse.BAM = self                 # May need within the below calls.
        self.process_bam_header()       # Quick but may have error exceptions
        self.process_bam_body()         # Little bit longer (10 seconds max) but key ones like read length, etc.

        # We have three longer running samtools stats commands we can run. But IDXStats is only a
        # second if a BAM file with an index.  So try as it enables the rest of the buttons in the GUI
        # (change_status_of_actiona_buttons).  May "raise" a warning.
        self.get_samtools_stats(button_directly=False)

        #  There are three get**stats calls. Each which will call process**stats if the file exists (or it creates it
        #  due to the button_directly=True or the special condition given above).  Will call all three get**stats calls
        #  when clicking on the Stats button (with button_directly=True on samtools_stats only). Idea is if any stats
        #  files exist from a previous run, will immediately read and process them.

        #return self                    # really the default for a class __init__ but just reminding ourselves


    def clear_BAM_file_settings(self):
        """ Clear any BAM file settings (simply delete class object; commit suicide) """
        del self
        # moved creating default, initial values to self.set_default_BAM_file_settings()


    def set_default_BAM_file_settings(self):
        """ Sets default values for all global BAM file settings; just as a safety to set all at the creation time. """

        self.file_oFN  = None   # Various renditions of the BAM/CRAM file name, path, etc
        self.file_oFPB = None
        self.file_FN   = None
        self.file_qFN  = None
        self.file_FBS  = None
        self.file_FPB  = None
        self.file_FP   = None
        self.file_FB   = None
        self.file_FS   = None

        self.file_type = None    # for "SAM", "BAM", or "CRAM" (file_FS without the leading dot)

        self.Header    = ""      # Will save complete BAM Header in text form (samtools view -H)
        self.file_stats = None   # File Stat result (for file size and other items)
        self.relfsize  = None    # Relative file size (to 45 GB) (to scale time values)

        self.Refgenome = None    # {hg,GRCh}{19/37,38}, hs37d5, hs38DH (not all combinations valid)
        self.RefgenomeNew = None # Following Reference Model study: hg, 1kg, ebi, ncbi
        self.Refgenome_qFN = None  # File name of reference genome used for BAM
        self.RefMito   = None    # rCRS, Yoruba, RSRS
        self.SNTypeC   = None    # "Chr", "Num" or "Acc" (NCBI Accession #)
        self.SNTypeM   = None    # M or MT
        self.SNCount   = 0       # Count of '@SQ:\tSN' entries in Header
        self.Build     = None    # Major Build (19?, 36, 37, 38)
        self.ReadType  = None    # Read type: Paired-end or Single-end

        self.gender    = None

        self.Stats    = False       # Will delay calculating basic stats if file not indexed
        self.Sorted   = False
        self.Indexed  = False

        # todo The following values are immutable once set; could be a tuple or fixed-key dictionary instead
        self.raw_gigabases              = 0
        self.raw_avg_read_depth_WES     = 0  # WES (Exome) region ARD (full WES region)
        self.raw_avg_read_depth_full    = 0  # Traditional, Full (old style) includes N's in model
        self.raw_avg_read_depth_NoN     = 0  # No N's included (new)
        self.raw_segs_read              = 0  # Total number of segments read (in BAM, or FASTQs)
        self.raw_avg_read_length        = 0
        self.mapped_gbases              = 0
        self.mapped_avg_read_depth_WES  = 0  # WES (Exome) region ARD (only tested portion)
        self.mapped_avg_read_depth_full = 0  # Traditional, Full (old style) includes N's in model
        self.mapped_avg_read_depth_NoN  = 0  # No N's included (new)
        self.mapped_segs_read           = 0  # Total # segments mapped (localized and unlocalized)
        self.mapped_reads_percent       = 0
        self.coverage                   = None  # Breadth of Coverage for whole BAM (valid value may be zero)
        self.coverage_WES               = None  # Breadth of Coverage for WES (valid value may be zero)
        self.low_coverage               = False
        self.long_read                  = False

        self.stats_chroms = []  # Will save IDXStats based summary table for display with "Show Stats" button
        self.chrom_types = {"A": 0, "X": 0, "Y": 0.0001, "M": 0, "*": 0, "O": 0}  # Stores # read segments for each
                            # Keeping Y >0 to avoid div0 error


    def process_bam_header(self):
        """ Reads and stores BAM header. Sets BAM.CoordSorted also. Does not need BAM Index file (.bai). """
        global askRefgenWindow

        # Rewritten to process more of header; and process header in Python and not BASH AWK/GREP scripts
        samtools = wgse.samtoolsx_qFN
        bamhead_qFN  = f'"{wgse.tempf.FP}bamheader.tmp"'
        bamfile  = self.file_qFN

        commands = f'{samtools} view -H --no-PG {bamfile} > {bamhead_qFN} \n'

        run_bash_script("GetBAMHeader", commands)

        bamhead_oFN = unquote(nativeOS(bamhead_qFN))
        if not os.path.exists(bamhead_oFN) or os.path.getsize(bamhead_oFN) == 0:  # Command success check
            raise BAMContentError('errBAMHeader')
        try:
            with open(bamhead_oFN, "r") as f:
                self.Header = f.read()
        except:
            raise BAMContentError('errBAMHeader')


        if "SO:coordinate" in self.Header:
            self.Sorted = True
        elif "SO:unsorted" in self.Header:
            self.Sorted = False
        else:
            DEBUG(f"ERROR: BAM coord sorted state in Header?: {self.Header[0]}")
        DEBUG(f"Bam Coord Sorted? {self.Sorted}")
        
        self.Indexed = self.check_for_bam_index()

        self.determine_reference_genome()       # Only need header available to determine reference_genome


    def process_bam_body(self):
        """
        Processes the BAM body; but mostly  via idxstats.  Operate on first 100K sequence entries.
        Called internally from the "Show Stats" button as well.
        Will not have run idxstats yet if detected unsorted or unindexed when initially set.
        """

        # First do some of our own custom processing; then lastly call samtools idxstats to finish up -- all very quick
        samtools = wgse.samtoolsx_qFN
        bamfile  = self.file_qFN
        rdlfile_qFN  = f'"{wgse.tempf.FP}readlength.tmp"'     # Sum of sequence field lengths divided by # of seq
        flagfile_qFN = f'"{wgse.tempf.FP}flags.tmp"'          # Flags field to get paired end
        cmd_head = wgse.headx_qFN
        cmd_awk  = wgse.awkx_qFN
        cmd_cat  = wgse.catx_qFN

        cram_opt = f'-T {self.Refgenome_qFN}' if self.file_type == "CRAM" else ""
        awkopt = "-F \"\\\"*\\t\\\"*\" '{sum += length($10); count++} END {if (count+0 != 0) print sum/count}'"  # REH 11Mar2020 FTDNA BAM Bug fix

        # Could make one line with "tee" but "tee" is not distributed with Win10 release tools
        commands = (
            f'{samtools} view {cram_opt} {bamfile} | {cmd_head} -100000 > {flagfile_qFN} \n'
            f'{cmd_cat} {flagfile_qFN} | {cmd_awk} {awkopt} > {rdlfile_qFN} \n'
        )

        run_bash_script('ButtonBAMStats2', commands)  # not running stats again; so simply do it

        # Process Flags file (currently for 0x1 mask; 1 to indicate paired-end or 0 for single-end reads
        flagfile_oFN = nativeOS(unquote(flagfile_qFN))
        if not os.path.exists(flagfile_oFN):
            raise BAMContentError('errBAMNoFlagsFile')
        with open(flagfile_oFN, "r") as flags_file:
            scale = 0  # We do not know which
            for flag_line in flags_file:
                line_columns = flag_line.split("\t")
                scale += 1 if int(line_columns[1]) & 0x1 else -1  # 2nd column is flags integer; if bit 0 is 1, is paired
               # Try to make sure is clearly one or the other

        self.ReadType = "Paired" if scale > 10000 else "Single" if scale < -10000 else "Unknown" # Schroedingers paradox
        DEBUG(f'Read Segment End Type: {self.ReadType} (scale:{scale})')
        if self.ReadType == "Unknown":
            raise BAMContentWarning('warnBAMBadReadType')

        # Process Read Length File (rdlfile)
        rdlfile_oFN = nativeOS(unquote(rdlfile_qFN))
        if not os.path.exists(rdlfile_oFN) or os.path.getsize(rdlfile_oFN) == 0:    # Catch script failures
            raise BAMContentError('errBAMNoReadLengthFile')
        with open(rdlfile_oFN, "r") as f:
            raw_avg_read_length = round(float(f.read(1000).strip().replace(",", ".")))  # REH 11Mar2020 style
        if raw_avg_read_length > 410:       # Increased to account for ySeq 400 bp read WGS now available
            # Nanopore long read; need to gather more values as shorter ones are packed first: 2 million enough ?
            DEBUG("Long Read BAM ... reprocessing to get the read length using more read samples")
            commands = f'{samtools} view {cram_opt} {bamfile} | {cmd_head} -2000000 | {cmd_awk} {awkopt} > {rdlfile_qFN} \n'

            run_bash_script('ButtonBAMStatsLong', commands)    #

            if not os.path.exists(rdlfile_oFN) or os.path.getsize(rdlfile_oFN) == 0:
                raise BAMContentError('errBAMNoReadLengthFile')
            with open(rdlfile_oFN, "r") as f:
                raw_avg_read_length = round(float(f.read(1000).strip().replace(",", ".")))  # REH 11Mar2020 style
        DEBUG(f"Read Length: {raw_avg_read_length!s}")
        self.raw_avg_read_length = raw_avg_read_length


    def get_samtools_stats(self, button_directly=False):
        """
            Called immediately after (re)setting a BAM file to update the main window summary results display;
            samtools idxstats results are saved for the detailed per-chromosome table from the stats button
            While it is much faster with an index file available, it is not required.  The file must be sorted though.
        """

        if self.Stats:      # Set by a successful process_samtools_stats run at the end here
            return          # Already ran idxstats; no need to run again. Just display from stored values.

        # May be called before output directory is set; so use temporary files area if so (that must be set)
        path = wgse.tempf.FP if not wgse.outdir and wgse.outdir.FP is None else wgse.outdir.FP
        idxfile_qFN = f'"{path}{self.file_FB}_idxstats.csv"'
        idxfile_oFN = nativeOS(unquote(idxfile_qFN))

        # If idxstats file does not exist, then must create it. Could be 30+ minutes if CRAM or unindexed BAM but only
        #  do that if button_directly was True. Cannot use mod time as could be stats from before CRAM / BAM conversion
        if not (os.path.isfile(idxfile_oFN) and os.path.getsize(idxfile_oFN) > 500):
           #    and os.path.getmtime(idxfile_oFN) > os.path.getmtime(self.file_oFN)):

            # Unless user hit button directly, then stop this now as it may take awhile ...
            if not button_directly:
                return      # Moved all warnings out to BAMFile (init) and in here so allows to proceed if csv exists
                # if not self.Sorted or not self.Indexed:
                    # Warn about not collecting stats if not sorted and/or indexed.
                #    raise BAMContentWarning('warnBAMNoStatsNoIndex')
                # elif self.file_type == "CRAM":  # Do not need to warn about both; one is enough
                    # Warn without collecting stats on CRAM file.  Even if indexed, samtools idxstats on a CRAM reads the
                    # whole file.  Better to convert to BAM than do that.
                #    raise BAMContentWarning('warnCRAMNoStats')

            # So either Stats button was hit directly (and not run yet) or in a new, indexed BAM so will run it as if
            #  button hit directly (quick run) to enable full table stats so can set rest of GUI buttons

            # Create and run Samtools idxstats command on the BAM file
            samtools = wgse.samtoolsx_qFN
            bamfile = self.file_qFN

            commands = f'{samtools} idxstats {bamfile} > {idxfile_qFN} \n'

            # Seconds to minutes to hours ....
            title = "ButtonBAMStats" if self.Sorted and self.Indexed and self.file_type == "BAM" else \
                    "ButtonBAMNoIndex" if self.Sorted else \
                    "ButtonBAMNoSort"
            run_bash_script(title, commands)

            # If still not there then report an error as could not create it when wanted to
            if not (os.path.isfile(idxfile_oFN) and os.path.getsize(idxfile_oFN) > 500):
                raise BAMContentError('errBAMNoIDXStatsFile')

        self.process_samtools_stats(idxfile_oFN)


    def process_samtools_stats(self, idxfile_oFN):
        """
        Only called from one place and only if the file exists and likely good. Passed the samtools idxstats result
        file to process and store in various settings.  One of three (possibly long) stats run commands and files.
        """

        # Only capture primary chromosomes; will append Alt Contigs later and save Unmapped in special summary variables
        #valid_autos = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
        #               "15", "16", "17", "18", "19", "20", "21", "22"]     # Moved to settings module
        #valid_other = ["X", "Y", "M"]  # moved to settings module although exclusively used here; removed '*'

        # These are the arrays to save the IDXstats data and additional processed values to display
        # Eight columns: keySN, storedSN, Modelen, Model N Count, Loc & Unlocalized Seg Counts,  (five initially)
        #                Mapped Gigabases (calc), Mapped Average Read Depth (calc),
        #                Breadth of Coverage (added later; just add 0 now)
        # SNCount rows: SNs but not including the special IDXstats row of unmapped ('*')
        stats_autos = []
        stats_other = []  # REH 20Mar2020 pulled out to easily numeric sort autosomes: X, Y, MT
        stats_altcont = ["O", wgse.lang.i18n['OtherChr'], 0, 0, 0]  # Alt Contigs row collecting all others
        # stats_unmap = ["*", "*", 0, 0, 0]    # Special stats_unmap row to capture special IDXStats unmapped value
        stats_total = ["T", wgse.lang.i18n['Total'], 0, 0, 0]  # Total row at bottom; before appending new columns
        stats_unmap = 0

        # To store different types of chromosome totals: Autosome, X, Y, Mitochondrial, Unmapped, and Other (alt contigs)*
        chrom_types = {'A': 0, 'X': 0, 'Y': 0.0001, 'M': 0, 'O': 0, '*': 0}  # non-zero to handle divide by zero later

        # Read the idxstats CSV table generated earlier; keep autosomes and others separate
        with open(idxfile_oFN, "r") as stats_file:
            for stats_line in stats_file:
                line_columns = stats_line.split("\t")
                chromosome = line_columns[0].strip()  # Save as found in file for printout, etc
                chromnum = chromosome.upper().replace("CHR", "").replace("MT", "M") # normalize to numeric only, M only
                len_in_model_bp = int(line_columns[1])
                map_seg_read = int(line_columns[2]) + int(line_columns[3].strip())  # Sum localized and unlocalized

                # Only process and store those with mapped segment values (and "*"); due to non-WGS BAMs (e.g. y-only)
                # Stats arrays are [nominal Chr SN, actual Chr SN in file, len in model, # map seg]
                # The nominal is consistent for computer processing (uppercase, no "chr" prefix, MT -> M, etc)
                if map_seg_read > 0:
                    if chromnum in wgse.valid_autos:
                        Ncount = wgse.nadjust[self.Build][chromnum]
                        stats_autos.append([chromnum, chromosome, len_in_model_bp, Ncount, map_seg_read])
                    elif chromnum in wgse.valid_other:
                        Ncount = wgse.nadjust[self.Build][chromnum] if chromnum != "*" else 0
                        stats_other.append([chromnum, chromosome, len_in_model_bp, Ncount, map_seg_read])
                    elif chromnum == '*':
                        stats_unmap = map_seg_read   # unlocalized column 3 has unmapped value total; save special
                    else:  # Alt Contigs, etc that are mapped (could even be non-std named autosomes)
                        if not chromosome in ["EBV", "NC_007605"]:  # Ignore EBV in 1KG as throws off stats; not mapped
                            stats_altcont[2] += len_in_model_bp
                            stats_altcont[3] = 0                # We have not loaded N count for contigs; only prinary
                            stats_altcont[4] += map_seg_read
                # Todo Need pythonic indices; not numeric position. For both read idxstats file and values table built

        DEBUG(f"Other_values: {str(stats_altcont)}")

        # If BAM is empty of mapped human genome SN's or unmapped-only (i.e. only alt-contigs),
        #   then raise / report error of BAM content as we do not know what we have
        if not (stats_autos or stats_other or stats_unmap > 0):
            DEBUG("ERROR: BAM has no human genome model key elements, what to do?")
            self.gender = None
            self.chrom_types = chrom_types  # chrom_types already zero'ed out
            raise BAMContentError('errBAMNonHumanGenome')

        # Numerically sort the Autosomes by chromosome number before building full table    # REH 20Mar2020
        # Most files are sorted alphabetically and not numerically. So chr10 follows chr1.
        stats_autos.sort(key=lambda chr: int(chr[0]))  # Sort entries based on chromosome #

        # Build table
        stats_chroms = stats_autos  # start with autosomes
        stats_chroms.extend(stats_other)  # Add somal, MT, other (alt_contigs, etc) but not '*' anymore
        if (stats_altcont[4] > 0):   # REH 09Mar2020 so only MAPped alt contigs included
            stats_chroms.append(stats_altcont)
        #if (stats_unmap[4] > 0):     # Removed unmap row from internal table (and display)
        #    stats_chroms.append(stats_unmap)

        # Develop column totals of table; column order first; append as last row to end of table (now before adding other)
        for row in range(len(stats_chroms)):
            #if stats_chroms[row][0] == '*':        # Skip unmapped row in totals
            #    break                              # Removed unmapped row from internal table (and display)
            for col in range(2, 5):  # Skip first 2 columns of chromosome labels
                stats_total[col] += int(stats_chroms[row][col])  # stats_total initialized to zero at start
        stats_chroms.append(stats_total)       # Add MAPPED, chromosome & mt Total row to end of table

        # Calculate summary columns (Mapped gbases, mapped avg read depth, coverage); append to each row
        for i in range(len(stats_chroms)):  # Fall back to C / Assembly, use i for row index :)
            len_in_model_bp = stats_chroms[i][2] + 0.0001   # make non-zero to help with divide by zero; especially '*'
            nlen_bp = stats_chroms[i][3]
            map_seg_read = stats_chroms[i][4]
            temp_mapped_gbases = float(map_seg_read * self.raw_avg_read_length)
            stats_chroms[i].append(round(temp_mapped_gbases / (10**9), 2))      # gbases column -- col 5
            stats_chroms[i].append(round(temp_mapped_gbases / (len_in_model_bp - nlen_bp)))  # ARD (no N) -- col 6
            stats_chroms[i].append(0)       # Save a zero'ed out coverage row to be filled in later -- col 7

            # Calculate the total segments for each type; to understand if WGS, Y-only, MT-only, unmapped only, etc
            if stats_chroms[i][0] != "T":       # Skip Total row
                # We lump all autosomes into A; assume always all or none. Otherwise, leave others (X, Y, M, O)
                #   alt contigs already lumped in O for Other
                name = "A" if stats_chroms[i][0] in wgse.valid_autos else stats_chroms[i][0]
                chrom_types[name] += stats_chroms[i][4]  # Total the mapped of each type

        chrom_types['*'] = stats_unmap   # Capture unmapped reads here (no longer kept in main array)
        DEBUG(f"Chrom_types: {str(chrom_types)}")

        # Give nice names to use in formulas below
        #  note: have to recalc Avg Read Depth as already made into string with added 'x'
        total_len_in_model = stats_chroms[-1][2] + 0.0001  # for divide by zero
        total_nlen_in_model = stats_chroms[-1][3]
        total_map_seg_read = stats_chroms[-1][4]    # Note: will not include '*' already
        total_unmap_seg_read = stats_unmap
        total_map_gbases = stats_chroms[-1][5] if total_map_seg_read > 0 else 0.0
        total_segs_read = total_map_seg_read + total_unmap_seg_read + 0.0001  # for divide by zero when unmapped BAM
        total_bases = total_segs_read * self.raw_avg_read_length

        # Set global stats that will be retained, displayed and used for decision making
        self.raw_gigabases = total_bases / (10**9)
        self.raw_avg_read_depth_NoN = 0 if total_len_in_model < 0.001 else \
            round(total_bases / (total_len_in_model - total_nlen_in_model))
        self.raw_avg_read_depth_full = 0 if total_len_in_model < 0.001 else \
            round(total_bases / total_len_in_model)     # Old style; may stop showing
        self.raw_segs_read = total_segs_read
        # self.total_reads_percent = total_segs_read / total_segs_read  #  Will always be 1 or 100%
        # self.raw_avg_read_length = raw_avg_read_length    # Calculated before this call and stored then
        self.mapped_gbases = total_map_gbases
        self.mapped_avg_read_depth_NoN = \
            round(total_map_seg_read * self.raw_avg_read_length / (total_len_in_model - total_nlen_in_model))
        self.mapped_avg_read_depth_full = \
            round(total_map_seg_read * self.raw_avg_read_length / total_len_in_model)   # Old style; may stop showing
        self.mapped_segs_read = total_map_seg_read
        self.mapped_reads_percent = total_map_seg_read / total_segs_read
        # self.mapped_avg_read_length       # is the same as raw

        self.stats_chroms = stats_chroms  # Save for detailed stats button display in mainWindow routines
        self.chrom_types = chrom_types
        if chrom_types["Y"] < 1 and chrom_types["X"] < 1:  # No X and Y reads
            self.gender = 'Unknown'
        elif (chrom_types["X"] / chrom_types["Y"]) > 20:
            self.gender = 'Female'
        else:   # Males: num Y reads ~x4 num X reads, Females: X > x20 Y reads
            self.gender = 'Male'

        # Process and post a few warnings if appropriate
        if self.mapped_avg_read_depth_NoN < 10 and \
                (self.chrom_types["A"] > 1 or self.chrom_types["Y"] > 1 or self.chrom_types["M"] > 1):
            """ Assume X is same as A for this check. If unmapped ('*') only, then mapped avg read depth is 0. """
            wgse_message("warning", 'LowCoverageWindowTitle', False, 'LowCoverageWarning')
            self.low_coverage = True

        if self.raw_avg_read_length > 410:
            wgse_message("warning", 'LongReadSequenceTitle', False, 'LongReadSequenceWarning')
            self.long_read = True

        self.Stats = True


    def get_coverage_stats(self, window, button_directly=False):
        """
            Called from a button on the BAM Stats display page to fill in the Breadth of Coverage column by running
            samtools coverage command.  Takes another 30 minutes so only give it as a second step option. Is only
            internally called so simply need to add to the self.stats_chroms array and return.  mainWindow will
            redisplay appropriately once returned.

            Coverage command file header:
            #rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
            Stats are based on full length; not less N coverage. So do not use coverage, meandepth
            Instead use covbases divided by (Model Len - N Bases; in stats_chroms)
            Using meanbaseq and meanmapq is OK
        """

        covfile_qFN = f'"{wgse.outdir.FPB}_coverage.csv"'
        covfile_oFN = nativeOS(unquote(covfile_qFN))

        # Long enough operation; if have file from before then reuse
        if button_directly and \
           not (os.path.isfile(covfile_oFN) and os.path.getsize(covfile_oFN) > 100 and
                os.path.getmtime(covfile_oFN) > os.path.getmtime(self.file_oFN)):

            # Run samtools coverage
            samtools = wgse.samtoolsx_qFN
            bamfile = self.file_qFN
            cramopts = f'--reference {self.Refgenome_qFN}' if self.file_type == "CRAM" else ""

            commands = f'{samtools} coverage {cramopts} -o {covfile_qFN} {bamfile} \n'
            run_bash_script("CoverageStats", commands, parent=window)

        # Check if coverage file generated; even if generated, skip if unmapped-only BAM
        if not os.path.isfile(covfile_oFN) or os.path.getsize(covfile_oFN) < 100:
            if button_directly:
                # Error in generating Coverage Stats file
                pass    # Todo maybe should report via pop-up before returning if button_directly
            return
        elif len(self.stats_chroms) == 1:     # Is unmapped only BAM so coverage is 0 and pre-filled
            # Pretend we processed it because is empty (all zero's) anyway which has been defaulted already
            self.coverage = self.stats_chroms[-1][7]        # Setting coverage summary value (even 0) says we processed
            return

        self.process_coverage_stats(covfile_oFN)


    def process_coverage_stats(self, covfile_oFN):
        """
        Single place call to actuall process the coverage stats file (we either found or created)
        """

        # Note: Total row at stats_chroms[-1] must always exist; Other row at stats_chroms[-2] may exist
        with open(covfile_oFN, "r") as stats_file:
            stats_file.readline()        # Skip header line that we purposely left in for users perusing file
            for stats_line in stats_file:
                line_columns = stats_line.split("\t")
                found = False
                for i in range(len(self.stats_chroms)):
                    if str(line_columns[0]) == self.stats_chroms[i][1]:
                        self.stats_chroms[i][7] = int(line_columns[4]) / \
                                                  (self.stats_chroms[i][2] - self.stats_chroms[i][3])
                        self.stats_chroms[-1][7] += int(line_columns[4])
                        found = True
                        break
                if not found and str(line_columns[0]) not in ["*", "EBV", "NC_007605"] and \
                        self.stats_chroms[-2][0] == 'O':    # Must be Alt Contig lumped in Other that we count
                    self.stats_chroms[-2][7] += int(line_columns[4])
                    self.stats_chroms[-1][7] += int(line_columns[4])

            if self.stats_chroms[-2][0] == 'O':     # If Other row exists
                self.stats_chroms[-2][7] = self.stats_chroms[-2][7] / self.stats_chroms[-2][2]
            # Final finish up on Total row
            self.stats_chroms[-1][7] = self.stats_chroms[-1][7] / (self.stats_chroms[-1][2] - self.stats_chroms[-1][3])

            # Set summary value from total row at bottom; indicates this has been finished successfully
            self.coverage = self.stats_chroms[-1][7]


    def get_WES_stats(self, window, button_directly=False):

        covfile_qFN = f'"{wgse.outdir.FPB}_wescvg.csv"'
        covfile_oFN = nativeOS(unquote(covfile_qFN))

        # Long enough operation; if we have file from before then simply reuse
        if button_directly and \
           not (os.path.isfile(covfile_oFN) and os.path.getsize(covfile_oFN) > 50 and
                os.path.getmtime(covfile_oFN) > os.path.getmtime(self.file_oFN)):

            # Run samtools depth with WES Bed file to get avg read depth and coverage
            # The normal samtools coverage does not accept a bed file; samtools bedcov does not work well
            samtools = wgse.samtoolsx_qFN
            awk = wgse.awkx_qFN
            bamfile = self.file_qFN
            bedfile_FN = f'"{wgse.reflib.get_wes_bed_file_FN(self.Build, self.SNTypeC)}"'
            cramopts = f'--reference {self.Refgenome_qFN}' if self.file_type == "CRAM" else ""

            # Because the depth file is so large; we process it via a pipe and awk; save per chromosome for user perusal
            script = """'{ names[$1]=$1 ; if($3==0){zero[$1]++} else {nz[$1]++ ; sumnz[$1]+=$3} } END {
              printf("%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n","chr","zero","nonzero","sum nz","fract nz","avg nz","avg all");
              for (x in names) { totalbc = zero[x]+nz[x]+1 ; printf("%s\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",
              x,zero[x],nz[x],sumnz[x],nz[x]/totalbc,sumnz[x]/nz[x],sumnz[x]/totalbc) } }'"""

            commands = f'{samtools} depth -a -b {bedfile_FN} {cramopts} {bamfile} | {awk} {script} > {covfile_qFN}'

            run_bash_script("CoverageStatsWES", commands, parent=window)

        # Check if coverage file generated; even if generated, skip if unmapped-only BAM (only headers; 48 bytes)
        if len(self.stats_chroms) == 1:       # Is unmapped only BAM (total only row) so coverage is 0
            self.raw_avg_read_depth_WES = 0     # these were not defaulted to zero
            self.mapped_avg_read_depth_WES = 0  #
            self.coverage_WES = 0               # Setting coverage summary value, even if 0, says we processed the data
        elif not os.path.isfile(covfile_oFN) or os.path.getsize(covfile_oFN) < 50:     # header only is 48 bytes
            if button_directly:
                pass        # Todo maybe should report via pop-up before returning
                # Error in generating WES Coverage Stats file
        else:
            self.process_WES_stats(covfile_oFN)


    def process_WES_stats(self, covfile_oFN):
        """
        Process our custom WES stats file created from a run of samtools depth and a custom awk script post-processor.
        Even though the file has per sequence values, we only want the summary values saved and displayed.
        Have gotten here from single internal call that has already assured the file exists.  Just filling in the data
        table from the file content.
        """
        # For WES Stats, only save total average read depth (zero count and non-zero); as well as Coverage. Not per chr
        zero_bases_cnt = nonzero_bases_cnt = sumnz_bases = 0
        with open(covfile_oFN, "r") as stats_file:
            next(stats_file)        # Skip header line that we added in for user persuing file
            for stats_line in stats_file:
                columns = stats_line.split("\t")
                zero_bases_cnt    += int(columns[1])
                nonzero_bases_cnt += int(columns[2])
                sumnz_bases       += int(columns[3])

            total_bases = zero_bases_cnt + nonzero_bases_cnt
            self.raw_avg_read_depth_WES    = round(sumnz_bases / total_bases)
            self.mapped_avg_read_depth_WES = round(sumnz_bases / nonzero_bases_cnt)
            self.coverage_WES              = nonzero_bases_cnt / total_bases


    def determine_reference_genome(self):
        """
        Try to determine the reference genome from just the header information.
        Has gotten more complicated since we accept any BAM; including subsets we produce (Y, MT, unmapped only).

        The BAM Header is set from the reference genome by the alignment tool creating the BAM.  It is generally
        carried through in full when subsets of the BAM are made (e.g. unmapped BAM).  So you may get an answer here
        even if only an unmapped BAM body content.  Header SQ records should be based only on the reference model and
        not alignment content.

        Original, generic, major model naming comes from:
         hg if old style "chr22/chrM" sequence names, GRCh if numeric-only
         19 if old style hg with Yoruba, 37 if build 37 model lengths, 38 if Build 38 model lengths
         Number of Sequence Names in the header is a good indicator of class and other model characteristics
        """

        # Determine SN names as given as opposed to relying on reference model to determine
        self.SNTypeC = "Chr" if "@SQ\tSN:chr" in self.Header else "Acc" if "@SQ\tSN:CM0006" in self.Header else "Num"
        self.SNTypeM = "MT" if "@SQ\tSN:MT" in self.Header or "@SQ\tSN:chrMT" in self.Header else "M"
        self.SNCount = self.Header.count('@SQ\tSN:')   # (Reliable) analysis model (file) by count of SN entries?
        DEBUG(f"SN: {self.SNTypeC}, {self.SNTypeM}; SN#:{self.SNCount}")

        # Possible header has been shortened if subset of BAM file; so use an Autosome (chr1), Y and X to check length
        # If only MT or Unmapped (*) in header, then cannot tell RefGenome
        Build37 = any(x in self.Header for x in ["LN:249250621", "LN:59373566", "LN:155270560"])
        Build38 = any(x in self.Header for x in ["LN:248956422", "LN:57227415", "LN:156040895"])
        Build36 = None          # Todo Should add Build36 here as still see them infrequently
        self.Build = 37 if Build37 else 38 if Build38 else 36 if Build36 else None
        DEBUG(f"Build: {self.Build:3d}")

        # Look for Mitochondria model to determine its type. This model is independent of the Major / Minor ref model
        self.RefMito = "Yoruba" if any(x in self.Header for x in ["M\tLN:16571", "MT\tLN:16571"]) else\
                       "rCRS"   if any(x in self.Header for x in ["M\tLN:16569", "MT\tLN:16569"]) else\
                       None     # Todo handling RSRS model -- look for spacers?
        if self.RefMito == "Yoruba":
            self.Build = 19 if Build37 else 18 if Build36 else 20   # Use old style HG build #'s if Yoruba
            # All Build 36 are Yoruba though so ....

        DEBUG(f"Ref Genome Mito: {self.RefMito}")

        # New form is based on Major / Class mechanism in Reference Genome study: https://bit.ly/34CO0vj
        #  Used to rely on SNCount.  But oddball, unrecognized ref genomes in BAMs may have something close set by user.
        #  So, rely more on key items found in special ref genomes.
        if "SN:NC_007605" in self.Header and self.SNCount in [84, 85, 86]:
            # hs37 class all have SN:NC007605 (EBV in Numeric naming)
            # human_1kg model is 84 and one less SN than hs37.fa.gz; hs37d5 is 86 and one more than hs37 (hs37d5 decoy)
            self.Refgenome = "hs37d5" if "@SQ\tSN:hs37d5" in self.Header else "hs37"
            self.RefgenomeNew = "1k37g"
        elif "SN:chrEBV" in self.Header and self.SNCount in [195, 456, 3366]:
            # hs38DH is 3366 SN count and has SN:HLA- unique; hs38 is 195 and hs38a is xxx; all uniquely have chrEBV
            self.Refgenome = "hs38DH" if "SN:hs38dh" in self.Header else "hs38"         # SN:HLA- also unique to hs38dh
            self.RefgenomeNew = "1k38h"
        elif self.SNCount in [25, 93, 297, 455, 639]:    # SNCounts of 6 hg/ebi models; 2 duplicated
            # Handle the hg (UCSC) and GRCh (EBI) models here
            # Note: new technique not used here only gives 19 to Build37 with Yoruba Mito model
            self.Refgenome  = "hg" if self.SNTypeC == "Chr" else "GRCh"
            self.Refgenome += str(self.Build)   # Build already takes into account mito for 19/37 differentiation
            self.RefgenomeNew = self.Refgenome.replace("GRCh","EBI") + ("g" if self.SNTypeC == "Chr" else "h")
        DEBUG(f'Ref Genome: {self.Refgenome}, by New nomenclature: {self.RefgenomeNew}')

        # Using the MD5 on the BAM header, see if we know the refgenome more refined (actual model file); for CRAM later
        # Note: MD5 method requires full WGS header; so the first above using chr1 would have worked at minimum

        pass  # Todo implement code to lookup MD5Sum on BAM Header (J Rhys)

        # We give up if still not set, simply ask the user
        if not self.Refgenome:  # not set yet: and (self.chrom_types["A"] > 1 or self.chrom_types["Y"] > 1):
            self.ask_reference_genome(in_align=False)
            DEBUG(f"Ref Genome (User): {self.Refgenome}")
            # todo ask user to send BAM header so we can get its md5 signature for future

        self.Refgenome_qFN = wgse.reflib.get_reference_genome_qFN(self.Refgenome, self.SNCount)
        tempFN = unquote(self.Refgenome_qFN)
        DEBUG(f'Ref Genome File: {tempFN}')


    # Todo should this be in mainwindow? or referencelibrary?  Window creation in bamfiles? As it is, they will likely
    #  have less a clue than we can determine (when not in the alignment command)
    def ask_reference_genome(self, in_align=False):
        """ If we cannot figure the reference genome, ask the user.  Gives them option to say "unknown". """

        # mainWindow.withdraw()
        rbRefgenome = StringVar(wgse.window)

        # Technically, may not be mainwindow that we are forked from; could be anywhere. But have no current window.
        askRefgenWindow = Toplevel(wgse.window)
        askRefgenWindow.transient()
        askRefgenWindow.protocol("WM_DELETE_WINDOW", 0)
        askRefgenWindow.title(wgse.lang.i18n['SelectReferenceGenome'])
        # askRefgenWindow.geometry("")
        askRefgenWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
        askRefgenWindow.columnconfigure(0, weight=1)
        askRefgenWindow.rowconfigure(0, weight=1)

        textToDisplay = ""
        if not in_align:        # Failed trying to determine Reference Genome of BAM; so ask user
            textToDisplay = wgse.lang.i18n['CouldntDetermineReferenceGenome'].replace("{{SN}}", str(self.SNCount))
            textToDisplay += wgse.lang.i18n['CautionSelectTheCorrectRefGenome']
        textToDisplay += '\n' + wgse.lang.i18n['PleaseSelectReferenceGenome']

        Label(askRefgenWindow, text=textToDisplay, font=font['14'],
              anchor="w", justify="left").grid(column=0, row=0, columnspan=2, padx=5, pady=3)
        Label(askRefgenWindow, text=self.refgenome_str(), font=font['14'],
              anchor="w", justify="center").grid(column=0, row=1, columnspan=2, padx=5, pady=3)

        # Setup series of radio buttons for available reference genomes; disable buttons that clearly cannot
        #  be (incorrect Build or Chr Name type)
        # Note: hs37 and Human_g1k_v37 are slightly different; by one SN value.  But we do not have a separate internal
        # value for human_g1k yet so we use hs37 to represent human_g1k
        rowcnt = 2
        rbhs37d5 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHS37D5'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="hs37d5")
        rbhs37d5.grid(column=0, row=rowcnt, padx=5, pady=3)
        if not (self.Build == 37 and self.SNTypeC == "Num"):
            rbhs37d5.configure(state='disabled')
        rbhs37   = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHumanG1k'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="hs37")
        rbhs37.grid(column=1, row=rowcnt, padx=5, pady=3); rowcnt += 1
        if not (self.Build == 37 and self.SNTypeC == "Num"):
            rbhs37.configure(state='disabled')
        rbhs38   = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHS38'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="hs38")
        rbhs38.grid(column=0, row=rowcnt, padx=5, pady=3)
        if not (self.Build == 38 and self.SNTypeC == "Chr"):
            rbhs38.configure(state='disabled')
        rbhs38dh = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHS38DH'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="hs38DH")
        rbhs38dh.grid(column=1, row=rowcnt, padx=5, pady=3); rowcnt += 1
        if not (self.Build == 38 and self.SNTypeC == "Chr"):
            rbhs38dh.configure(state='disabled')
        rbhg37   = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHG37'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="hg37")
        rbhg37.grid(column=0, row=rowcnt, padx=5, pady=3)
        if not (self.Build == 37 and self.SNTypeC == "Chr"):
            rbhg37.configure(state='disabled')
        rbhg19   = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHG19'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="hg19")
        rbhg19.grid(column=1, row=rowcnt, padx=5, pady=3); rowcnt += 1
        if not (self.Build == 19 and self.SNTypeC == "Chr"):
            rbhg19.configure(state='disabled')
        rbhg38   = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHG38'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="hg38")
        rbhg38.grid(column=0, row=rowcnt, padx=5, pady=3)
        if not (self.Build == 38 and self.SNTypeC == "Chr"):
            rbhg38.configure(state='disabled')
        rbhg37wg = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsHG37wgse'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="hg37wg")
        rbhg37wg.grid(column=1, row=rowcnt, padx=5, pady=3); rowcnt += 1
        if not (self.Build == 37 and self.SNTypeC == "Chr"):
            rbhg37wg.configure(state='disabled')
        rbgrch37 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsGRCh37'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="GRCh37")
        rbgrch37.grid(column=0, row=rowcnt, padx=5, pady=3)
        if not (self.Build == 37 and self.SNTypeC == "Num"):
            rbgrch37.configure(state='disabled')
        rbgrch38 = Radiobutton(askRefgenWindow, text=wgse.lang.i18n['RefGenomIsGRCh38'], font=font['14'],
                               indicatoron=0, variable=rbRefgenome, value="GRCh38")
        rbgrch38.grid(column=1, row=rowcnt, padx=5, pady=3); rowcnt += 1
        if not (self.Build == 38 and self.SNTypeC == "Num"):
            rbgrch38.configure(state='disabled')
        Radiobutton(askRefgenWindow, text=wgse.lang.i18n['Unknown'], font=font['14'],
                    indicatoron=0, variable=rbRefgenome, value="Unknown"
                    ).grid(column=0, row=rowcnt, columnspan=2, padx=5, pady=3); rowcnt += 1
        Button(askRefgenWindow, text=wgse.lang.i18n['Done'], font=font['14'],
               command=lambda win=askRefgenWindow:self.ask_reference_genome_finish(win, rbRefgenome.get())
               ).grid(column=0, row=rowcnt, columnspan=2, padx=5, pady=3); rowcnt += 1

        askRefgenWindow.update()
        askRefgenWindow.grab_set()          # ala mainloop for Toplevel, transient windows
        askRefgenWindow.wait_window()


    def ask_reference_genome_finish(self, window, Refgenome):
        # from mainwindow import resume_main_window
        global askRefgenWindow

        self.Refgenome = Refgenome
        if self.Refgenome == "Unknown":
            self.Refgenome = None
        if not self.RefMito:
            self.RefMito = "Yoruba" if self.Refgenome == "hg19" else\
                           "rCRS"   if self.Refgenome else\
                           None
        self.RefgenomeNew = None    # Todo How to set up new refgenome after the fact?

        window.destroy()
        # resume_main_window()


    def chrom_types_str(self):
        """
        Create a string of comma-separated BAM SN components: Auto(somes), X, Y, Mito, etc.
        Typical 30x WGS male BAM file would return: "Auto, X, Y, Mito, Unmap, Other" (dependent on language)
        """
        count = 0
        chrom_types_str = ""
        for key, value in self.chrom_types.items():
            if (value > 1 and not (key == 'Y' and self.gender == 'Female')):  # Do not show Y on Female
                count += 1
                chrom_types_str += (", " if chrom_types_str != "" else "") + wgse.lang.i18n[key]
                ''' Key values are one of: ("A", "X", "Y", "M", "O", "*"). Single character lookup in language file. '''
        if count == 1:
            chrom_types_str += f' {wgse.lang.i18n["Only"]}'
        return chrom_types_str


    def refgenome_str(self):
        """ Create string of ref genome, mito genome, and SN Count  """

        result  = self.Refgenome if self.Refgenome else wgse.lang.i18n["Unknown"]
        if self.SNTypeC:
            result += f' ({self.SNTypeC})'
        if self.RefMito:
            result += ", " if result else ""
            result += self.RefMito
        if self.SNCount:
            result += ", " if result else ""
            result += f'{self.SNCount} {wgse.lang.i18n["SNs"]}'
        return result


    def filestatus_str(self):
        """ Create string of sorted, indexed and file size status """
        result = ""
        #if self.Sorted:
        result += ", " if result else ""
        result += wgse.lang.i18n["Sorted"] if self.Sorted else wgse.lang.i18n["Unsorted"]
        #if self.Indexed:
        result += ", " if result else ""
        result += wgse.lang.i18n["Indexed"] if self.Indexed else wgse.lang.i18n["Unindexed"]
        # File size
        result += ", " if result else ""
        result += f'{round(self.file_stats.st_size / (10.0 ** 9), 1)} {wgse.lang.i18n["GBs"]}'
        return result


    def get_chr_name(self,chrtype):
        """ Return chromosome name of "type" (MT or Y) specific to the BAM reference model. """
        if self.SNTypeC == "Chr":
            chrM = "chrM" if self.SNTypeM == "M" else "chrMT"
            chrY = "chrY"
        elif self.SNTypeC == "Num":
            chrM = "MT" if wgse.BAM.SNTypeM == "MT" else "M"
            chrY = "Y"
        else:
            chrM = chrY = ""
        return chrM if chrtype == "M" else chrY


    def check_for_bam_index(self):
        """ Marko's original BAM Index file exists check.  """
        # todo need to add checks in OUT directory, if specified, as well. And add option to use in samtools calls if there
        self.Indexed = \
            (os.path.isfile(self.file_oFN + ".bai") or os.path.isfile(self.file_oFN + ".BAI")) if self.file_type == "BAM" else \
            (os.path.isfile(self.file_oFN + ".crai") or os.path.isfile(self.file_oFN + ".CRAI")) if self.file_type == "CRAM" else None
        # todo seems case sensitive check is naive for case sensitive file systems. how does samtools handle?
        return self.Indexed


    def find_FASTQs(self):
        """
        A special to search for commonly-named FASTQs from vendor.  As a precursor to creating them from the BAM.
        In here as if they are found, will set class variable so they can be used and associated.  Likely will do
        similar for VCF.  Eventually, add "get_" routine like for reference_genome to ask the user. Saves time and
        disk space to reuse. Helpful especially when we start saving and restore BAM library stats.
        """

        # First look in output directory for ones we create (our naming convention)
        if wgse.outdir.FPB:
            r1fastq_FN = f'{wgse.outdir.FPB}_R1.fastq.gz'
            r2fastq_FN = f'{wgse.outdir.FPB}_R2.fastq.gz'
            if os.path.exists(nativeOS(r1fastq_FN)) and os.path.exists(nativeOS(r2fastq_FN)):
                self.R1fastq_FN = r1fastq_FN
                self.R2fastq_FN = r2fastq_FN
                return True  # todo checking dates and sizes

        # Now look in BAM directory for common test company names: Nebula Genomics
        r1fastq_FN = f'{self.file_FPB}_1.fq.gz'
        r2fastq_FN = f'{self.file_FPB}_2.fq.gz'
        if os.path.exists(nativeOS(r1fastq_FN)) and os.path.exists(nativeOS(r2fastq_FN)):
            self.R1fastq_FN = r1fastq_FN
            self.R2fastq_FN = r2fastq_FN
            return True  # todo checking dates and sizes

        # Now look in BAM directory for common test company names: Dante Labs
        r1fastq_FN = f'{self.file_FPB}_SA_L001_R1_001.fastq.gz'
        r2fastq_FN = f'{self.file_FPB}_SA_L001_R2_001.fastq.gz'
        if os.path.exists(nativeOS(r1fastq_FN)) and os.path.exists(nativeOS(r2fastq_FN)):
            self.R1fastq_FN = r1fastq_FN
            self.R2fastq_FN = r2fastq_FN
            return True  # todo checking dates and sizes

        # Now look in BAM directory for common test company names: ySeq (does not provide FASTQ's)

        # todo ask the user if they know where they are and let them set.  Implement once stats save happens.

        self.R1fastq_FN = self.R2fastq_FN = None
        return False


    def realign_BAM_filename(self, new_refgenome):
        """
        For use when realigning BAM to different build (37 to 38; or 38 to 37)
        To avoid lengthening the BAM name, see if we can simply replace the build ID already embedded in the file name
        """
        # Will never go back to hg19 in current form as do not generate Yoruba target models
        if "GRCh38" in self.file_FB:
            newBAM_FN = self.file_FBS.replace("GRCh38", "GRCh37")
        elif "GRCh37" in self.file_FB:
            newBAM_FN = self.file_FBS.replace("GRCh37", "GRCh38")
        elif "_hg38" in self.file_FB:
            newBAM_FN = self.file_FBS.replace("_hg38", "_hg37")
        elif "_hg37" in self.file_FB:
            newBAM_FN = self.file_FBS.replace("_hg37", "_hg38")
        elif "_hg19" in self.file_FB:
            newBAM_FN = self.file_FBS.replace("_hg19", "_hg38")
        elif "_hs37d5" in self.file_FB:
            newBAM_FN = self.file_FBS.replace("_hs37d5", "_hs38")
        elif "_hs38" in self.file_FB:
            newBAM_FN = self.file_FBS.replace("_hs38", "_hs37d5")
        else:  # Could not find acceptable name to patch so append the new Refgenome to the current file name.
            newBAM_FN = f'{self.file_FB}_{new_refgenome}{self.file_FS}'
        # Need to place the new file in the Output Directory
        return f'{wgse.outdir.FP}{newBAM_FN}'
