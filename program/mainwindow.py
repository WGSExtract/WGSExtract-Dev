# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
    Main window module for WGS Extract.  Everything revolves around the main window and its buttons on
     different panes. Once various support modules are initialized in the main program (wgesextract), the
     main entry point is setup_mainWindow.  Most of the button click action is handled here also.
     Notable exceptions are the BAM file handling (module bamfile) and autosomal / microarray file extraction
     (module microarray, aconv and hg38to19). The module referencegenome is quickly expanding to a stand alone
     and handle all reference model downloads (all reference human genome models, liftover file, and likely
     soon the 300 MB of microarray templates; thus dramatically reducing the footprint of the distributed program.)
"""

import os.path      # os.path.isdir, .exists
import time         # time.ctime() for results window
import webbrowser   # webbrowser.opennew()
# Original code imported tk after ttk. Thus making tk have priority.
# We know because Tk uses parameter font while Ttk uses styles like CSS; Tk has Message, Ttk only has Label
from tkinter import Tk, Toplevel, filedialog, Button, Label, LabelFrame, LEFT, CENTER, N, S, E, W
from tkinter.ttk import Notebook, Frame

from PIL import ImageTk, Image   # , ImageGrab
import pyscreenshot as ImageGrab        # PIL.ImageGrab not available in Linux

import settings as wgse
font = wgse.font
from utilities import DEBUG, nativeOS, universalOS, unquote, wgse_message

from commandprocessor import run_bash_script
from bamfiles import BAMFile, BAMContentError, BAMContentWarning
from microarray import button_select_autosomal_formats, _button_CombinedKit


# Pop-up result windows
global statsWindow, yHaploResWindow, simResWindow, selectLanguageWindow
# global errPopupWindow (utilities), selectAutosomalFormatsWindow (microarray)
# global askRefgenWindow (bamfiles), pleaseWaitWindow (commandprocessor)


def resume_main_window():
    """ General function to resume / restore the main window after some user button action(s). """
    update_action_buttons()
    wgse.window.update()
    wgse.window.deiconify()
    wgse.tempf.clean()


def button_exit():
    """ Processing user button (or window close) as main exit out of the program; do any final cleanup. """
    wgse.save_settings()
    # wgse.tempf.list.append(f'{wgse.tempf.oFP}')      # Special case; to wipe out the PID created directory
    wgse.tempf.clean(True)
    wgse.window.destroy()
    exit()


def update_action_buttons():
    """
    Called to change the status of all main window action buttons; possibly BAM data dependent.

    Newstate is "normal" IF a BAM is set and the output directory is set.  Otherwise, is "disabled"
    Cannot enable Stats nor Header button as it has a Save button which needs output set
    """
    # Special buttons that are always active but change to value set; default value initially set
    # Could maybe change out to use StringVar so updated automatically ....
    outputDirectoryButton.configure(text=wgse.outdir.FB if wgse.outdir and wgse.outdir.FP else \
                                         wgse.lang.i18n['SelectOutputDirectory'])
    reflibDirectoryButton.configure(text=wgse.reflib.FB if wgse.reflib and wgse.reflib.set else \
                                         wgse.lang.i18n['Default'])
    tempDirectoryButton.configure(text=wgse.tempf.FB if wgse.tempf and wgse.tempf.set else \
                                       wgse.lang.i18n['Default'])
    wslbwaButton.configure(text=wgse.lang.i18n['Active'] if wgse.wsl_bwa_patch else \
                                wgse.lang.i18n['Inactive'])

    # Relies on global listing "all_action_buttons" that can be set active ("normal") or deactive ("disabled")
    newstate = "normal" if wgse.outdir.FP and wgse.BAM and wgse.BAM.Stats else "disabled"
    for single_button in all_action_buttons:
        single_button.configure(state=newstate)

    # Override BAM File area buttons if BAM and Output Directory are set even if Stats is not yet set
    if wgse.BAM and wgse.outdir.FP:
        set_BAM_window_settings()

    if newstate == "disabled":
        return
    # else, if newstate == "normal" then Stats have been run and the BAM and Output Directory are both set

    # Special override: Here only if enabled buttons (state = "Normal")
    #   Disable some buttons if data is not available in the BAM (relies on Stats having been run)
    if wgse.BAM.chrom_types["A"] <= 1:
        autosomalFormatsButton.configure(state="disabled")
        autosomalVCFButton.configure(state="disabled")
        runMicroParallelButton.configure(state="disabled")

    if wgse.BAM.chrom_types["M"] <= 1:
        mitoFASTAButton.configure(state="disabled")
        mitoBAMButton.configure(state="disabled")
        mitoVCFButton.configure(state="disabled")
        haplogroupMtButton.configure(state="disabled")
        yANDmtButton.configure(state="disabled")

    # If gender is Female, disable Y buttons (note: females have >1 mapped reads on Y!)
    if wgse.BAM.chrom_types["Y"] <= 1 or wgse.BAM.gender == 'Female':
        haplogroupYButton.configure(state="disabled")
        yANDmtButton.configure(state="disabled")
        yOnlyButton.configure(state="disabled")
        yVCFButton.configure(state="disabled")

    if wgse.BAM.chrom_types["*"] <= 1:
        exportUnmappedReadsButton.configure(state="disabled")


def clear_BAM_window_settings():
    """ Called after clearing a BAM to clear all the window labels of wgse.BAM.* settings. """

    # Title bar filename label
    titleFileLabel.configure(text="")

    # Reset button to initial Select message
    bamSelectButton.configure(text=wgse.lang.i18n['SelectBamFile'])

    bamReferenceGenomeLabel.configure(text="")
    bamMapAvgReadDepthLabel.configure(text="")
    # bamAverageReadLengthLabel.configure(text="")
    bamGenderLabel.configure(text="")
    bamChromsLabel.configure(text="")
    bamFsizeLabel.configure(text="")

    bamIdxstatsButton.configure(state="disabled")
    bamHeaderButton.configure(state="disabled")
    bamSortButton.grid_remove()
    bamIndexButton.grid_remove()
    bamConvertButton.grid_remove()
    bamRealignButton.grid_remove()


def set_BAM_window_settings():
    """
    Called after processing a new BAM to (re)set window text labels dependent on wgse.BAM.* settings
     Also called when simply index / unindex a BAM as it changes the Indexed state.
     Note that sort / unsort, to CRAM / BAM, and realign buttons all cause a new BAM to be set which then call this.
     Stats button may cause delayed idxstats to run (for CRAM, unindexed BAM) so button_close also calls here
       button_close reused by other results windows so causes a few extra calls but no harm
    """

    if not wgse.BAM:
        clear_BAM_window_settings()     # Just to make sure; but likely OK to just return
        return

    # For title bar filename (so visible no matter what the tab setting; as a reminder)
    titleFileLabel.configure(text=wgse.BAM.disp_FBS)

    # Changed from FB to FBS to show CRAM file extension; use size limited one to avoid overruns in display
    bamSelectButton.configure(text=wgse.BAM.disp_FBS)
    bamReferenceGenomeLabel.configure(text=wgse.BAM.refgenome_str())

    if wgse.BAM.Stats:      # Values to report as "blank" when no Stats run to fill them in yet
        bamMapAvgReadDepthLabel.configure(text=f'{wgse.BAM.mapped_avg_read_depth_NoN} {wgse.lang.i18n["xTimes"]}')
        # bamAverageReadLengthLabel.configure(text=f'{wgse.BAM.raw_avg_read_length} {wgse.lang.i18n["bp"]}')
        bamGenderLabel.configure(
                text=wgse.lang.i18n.get(wgse.BAM.gender, wgse.BAM.gender))  # use value if key error
    else:
        bamMapAvgReadDepthLabel.configure(text="")
        # bamAverageReadLengthLabel.configure(text="")          # Removed, ran out of room
        bamGenderLabel.configure(text="")

    bamChromsLabel.configure(text=wgse.BAM.chrom_types_str())   # Will be blank if no stats; now called File Content:

    bamFsizeLabel.configure(text=wgse.BAM.filestatus_str())     # Now called File Stats as has sort and index status:

    bamIdxstatsButton.configure(state="normal");  bamIdxstatsButton.grid()     # Always show
    bamHeaderButton.configure(state="normal");    bamHeaderButton.grid()       # Always show

    # For Sort and Index buttons that serve as status indicators also; disable if already done and not DEBUG mode
    #  otherwise, if DEBUG mode, always enabled but a toggle button
    if wgse.DEBUG_MODE and wgse.BAM.Sorted:   # For developers only
        bamSortButton.configure(state="normal", text=wgse.lang.i18n["Unsort"], command=_button_unsort_BAM)
    else:
        bamSortButton.configure(state="disabled" if wgse.BAM.Sorted else "normal",
                                text=wgse.lang.i18n["Sorted"] if wgse.BAM.Sorted else wgse.lang.i18n["Sort"],
                                command=button_sort_BAM)
    bamSortButton.grid()        # Make button visible because BAM File is set now
    if wgse.DEBUG_MODE and wgse.BAM.Indexed:   # For developers only
        bamIndexButton.configure(state="normal", text=wgse.lang.i18n["Unindex"], command=_button_unindex_BAM)
    else:
        bamIndexButton.configure(state="disabled" if wgse.BAM.Indexed else "normal",
                                 text=wgse.lang.i18n["Indexed"] if wgse.BAM.Indexed else wgse.lang.i18n["Index"],
                                 command=button_index_BAM)
    bamIndexButton.grid()       # Make button visible because BAM File is set now

    # Additional buttons that are initially invisible and should now be made visible
    if wgse.BAM.file_type == "BAM":
        bamConvertButton.configure(state="normal", text=wgse.lang.i18n["ToCRAM"], command=button_BAM_to_CRAM)
    elif wgse.BAM.file_type == "CRAM":
        bamConvertButton.configure(state="normal", text=wgse.lang.i18n["ToBAM"], command=button_CRAM_to_BAM)
    bamConvertButton.grid()  # Make button visible

    bamRealignButton.configure(state="normal")
    bamRealignButton.grid()     # Make button visible


def button_select_BAM_file(newBAM=None):
    """
        Processing user button to select a BAM File. Initiates basic processing of BAM file content to
        display stats in mainWindow.  Special handling as can have errors deep in the call stack. So use
        Try...Except with custom error class here. Handle and catch all known errors at this level.
        Otherwise, no general "except:"-ion to catch all errors (simply pass on coding errors).

        This can also be called internally to replace (or set new) a BAM file.  So must diligently handle
        checking and setup before doing actual replacement.  Only internal call passes the new BAM file name.
    """
    if newBAM:
        BAM_FN = newBAM  # Called internally with new BAM file to set (after button_BAM_sort, for example)
    else:           # Ask user for BAM as one not supplied; must be mainWindow button click to set BAM
        initialdir = os.path.dirname(wgse.BAM.file_FN) if wgse.BAM and wgse.BAM.file_FN else ""
        BAM_FN = filedialog.askopenfilename(parent=wgse.window, initialdir=initialdir,
                title=wgse.lang.i18n['SelectBamFile'],
                filetypes=((wgse.lang.i18n['BamFiles'], ("*.bam", "*.cram", "*.sam")),))
        if not BAM_FN:      # User must have cancelled; simply return as nothing to do
            DEBUG("BAM/CRAM/SAM File Selection Cancelled ...")   # Do nothing; do NOT unselect current one
            resume_main_window()
            return False        # Does not really matter; had to have been called from Button here
    DEBUG(f"Selected BAM/CRAM/SAM (request): {BAM_FN}")

    # We are going to try and process the BAM but not replace an existing one, if set, until verified is OK
    save_orig_BAM = wgse.BAM;  wgse.BAM = None
    try:
        wgse.BAM = BAMFile(BAM_FN)     # Process the BAM by opening, reading header and setting up Class instance
    except BAMContentError as err:
        wgse_message("error", 'InvalidBAMTitle', True,
                     wgse.lang.i18n['errBAMFile'].replace("{{BAMF}}", BAM_FN) + err.reason)
        if wgse.BAM:                    # If set, then delete due to error and restore old value ;
            wgse.BAM = save_orig_BAM    #   garbage collect the new BAM
        # Leave old BAM in place as new one failed; fall through as deleted temporary
    except BAMContentWarning as warn:       # Todo get better title for warning box; was InvalidBAMTitle
        wgse_message("warning", 'Overview', True,
                     wgse.lang.i18n['warnBAMFile'].replace("{{BAMF}}", BAM_FN) + warn.reason)
        # Simply fall through after issuing warning
    else:
        if not (wgse.BAM.Sorted and wgse.BAM.Indexed):
            # Warn about not collecting stats on initial BAM load if not sorted and/or indexed.
            wgse_message("warning", 'Overview', True,
                         wgse.lang.i18n['warnBAMFile'].replace("{{BAMF}}", BAM_FN) +
                         wgse.lang.i18n['warnBAMNoStatsNoIndex'])
        elif wgse.BAM.file_type == "CRAM":  # Do not need to warn about both; one is enough
            # Warn without collecting stats on CRAM file.
            wgse_message("warning", 'Overview', True,
                         wgse.lang.i18n['warnBAMFile'].replace("{{BAMF}}", BAM_FN) +
                         wgse.lang.i18n['warnCRAMNoStats'])
    finally:
        if not wgse.BAM:                    # If not set, restore old setting
            wgse.BAM = save_orig_BAM
        if wgse.outdir.FP and wgse.BAM:       # Sets part of output File variables that includes BAM base name
            wgse.outdir.oFPB = wgse.outdir.oFP + wgse.BAM.file_FB
            wgse.outdir.FPB  = wgse.outdir.FP  + wgse.BAM.file_FB
        else:
            set_BAM_window_settings()         # Sets window labels of BAM file settings
        wgse.save_settings()
        resume_main_window()


def button_stats_BAM():
    """
    Processing user button to show detailed BAM Stats. We have likely already run BAM.get_samtools_stats;
    just a matter of displaying the stored table.
    """
    global statsWindow

    # There are three total, possibly long stats runs.  Each 30-60 minutes long.  We delay long stats runs until
    #  the user requests it directly with a button hit.  But, if they were run previously and we find the file still
    #  available, we will go ahread and process the file content immediately.  Each stats calculation determines itself
    #  if run before and/or will run quickly if not being forced by a button click directly.

    # Stats Button only enabled if wgse.BAM exists; but Stats may not have been calculated yet so force it now ...
    wgse.BAM.get_samtools_stats(button_directly=True)
    #stats = wgse.BAM.Stats

    # So do not force now as if a button hit directly. Instead, determine if file from previous run (of program) and
    #  load the values if so for display; otherwise display the direct buttons to generate later here.
    wgse.BAM.get_coverage_stats(wgse.window, button_directly=False)
    coverage = False if wgse.BAM.coverage is None else True

    wgse.BAM.get_WES_stats(wgse.window, button_directly=False)
    WESstats = False if wgse.BAM.coverage_WES is None else True

    #
    # Create and setup main Stats window
    #

    # wgse.window.withdraw()        # Maybe hide main window while displaying stats? Turned off for now
    statsWindow = Toplevel(wgse.window)
    statsWindow.transient(wgse.window)
    statsWindow.title(wgse.lang.i18n['BamFileStatistics'])
    statsWindow.geometry(wgse.bamstats_winsize)
    statsWindow.protocol("WM_DELETE_WINDOW", lambda win=statsWindow: button_close(win))
    statsWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    statsWindow.columnconfigure(0, weight=1)
    statsWindow.rowconfigure(0, weight=1)

    fontti = font['18b']        # Title font for each frame
    fonth1 = font['14b']        # Header font
    fontb1 = font['14']         # Body font
    fontb2 = font['13']         # Special for by-chr body only
    fontb3 = font['12']         # special for date/time/prog_version; by name table entries

    # Setup by-Chromosome Frame (main table, left)
    byChromFrame = LabelFrame(statsWindow, text=wgse.lang.i18n['ByRefSeqName'], font=fontti)
    byChromFrame.grid(row=0, column=0, padx=5, pady=2)
    byChromFrame.columnconfigure(0, weight=1)
    byChromFrame.rowconfigure(0, weight=3)

    # Setup Summary Frame (column, right; with Close button)
    summaryFrame = LabelFrame(statsWindow, text=wgse.BAM.disp_FBS, font=fontti)
    summaryFrame.grid(row=0, column=1, padx=5, pady=2)

    # Header for by-Chromosome table
    table_header1 = ['RefSeqName', 'LenInModel', 'NCntInModel', 'ReadSegsMap', 'MappedGbases', 
                     'MapAvgReadDepthNoNShort'] # Missing Coverage because may be a button instead of label
    for col in range(len(table_header1)):
        Label(byChromFrame, text=wgse.lang.i18n[table_header1[col]],
              font=fonth1).grid(column=col, row=0, rowspan=2, padx=5, pady=0)    # Single row in grid but 2x high using \n in text
    if not coverage:        # Last column Coverage not filled in yet
        Button(byChromFrame, text=wgse.lang.i18n['Coverage'], command=lambda sw=statsWindow: button_statsCov_BAM(sw),
               font=fonth1).grid(column=len(table_header1), row=0, rowspan=2, padx=5, pady=0)
    else:                   # Last column Coverage now available
        Label(byChromFrame, text=wgse.lang.i18n['Coverage'],
              font=fonth1).grid(column=len(table_header1), row=0, rowspan=2, padx=5, pady=0)

    # Body for by-Chromosome table; including final total row
    for row in range(len(wgse.BAM.stats_chroms)):
        total_row = True if (wgse.BAM.stats_chroms[row][0] == 'T') else False
        asterisk = True if "Y" in wgse.BAM.stats_chroms[row][1] or \
                          ("X" in wgse.BAM.stats_chroms[row][1] and wgse.BAM.gender == "Male") else False
        width = 0;  val = ""
        # First column is key with real start in 2nd; but want to add blank label for missing last Coverage column
        for col in range(len(wgse.BAM.stats_chroms[row])-1):
            index = col+1
            if col == 0:                            # Chromosome / Sequence label
                val = '{:>5}'.format(wgse.BAM.stats_chroms[row][index])
                width = 5
            if col == 1:                            # Use locale specific grouping (matters if K OR actual; not M)
                large = int(wgse.BAM.stats_chroms[row][index]) > 10 ** 6  # When unmapped is too large for unloc segs
                val = '{:5n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -3)/1000,
                                       wgse.lang.i18n["Thous"]) if not total_row and not large else \
                      '{:5n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -6)/(10.0 ** 6),
                                        wgse.lang.i18n["Mill"])
                width = 7                           # Room for appended " K"; possible billion number in Total
            elif col in [2, 3]:                     # Use locale specific grouping (matters if K OR actual; not M)
                large = int(wgse.BAM.stats_chroms[row][index]) > 10 ** 6    # When unmapped is too large for unloc segs
                val = '{:5n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -3)/1000,
                                       wgse.lang.i18n["Thous"]) if not total_row and not large else \
                      '{:5n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -6)/(10.0 ** 6),
                                        wgse.lang.i18n["Mill"])
                width = 7                           # Room for appended " K"; common hundred million max value
            elif col == 4:                          # Process float gbases; skip '*'
                # '*' line may have null strings? (historically)
                val = '{:6.2f}'.format(float(wgse.BAM.stats_chroms[row][index])) \
                    if wgse.BAM.stats_chroms[row][index] else '0.00'
                width = 5
            elif col == 5:                          # Avg Read Depth (add 'x')
                xstar = f'{wgse.lang.i18n["xTimes"]}' + ("*" if asterisk else "")
                large = int(wgse.BAM.stats_chroms[row][index]) > 1000       # When Mito is too large for Mapped ARDepth
                val = '{:5n}{} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), -3)/1000,
                                          wgse.lang.i18n["Thous"], xstar) if large else \
                      '{:6n} {}'.format(round(int(wgse.BAM.stats_chroms[row][index]), 0), xstar)
                width = 5
            elif col == 6:                          # Breadth of Coverage (extra run / button)
                val = '{:6.2f} %'.format(round(wgse.BAM.stats_chroms[row][index] * 100, 2)) \
                    if wgse.BAM.coverage is not None else ""
                width = 8       # mtDNA is always 100%
            Label(byChromFrame, text=val, font=fonth1 if total_row else fontb2,
                  width=width, anchor="e").grid(column=col, row=row+2, padx=0, pady=0)      # Was 10, 0

    read_type = wgse.lang.i18n['PairedEnd'] if wgse.BAM.ReadType == "Paired" else \
                wgse.lang.i18n['SingleEnd'] if wgse.BAM.ReadType == "Single" else ""

    # Summary Frame to right, top section: Mapped / Raw
    crow= 0
    Label(summaryFrame, text=wgse.lang.i18n['Mapped'].upper(), font=fonth1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text=wgse.lang.i18n['Raw'].upper(), font=fonth1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['AverageReadDepthNoN'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=f'{wgse.BAM.mapped_avg_read_depth_NoN} {wgse.lang.i18n["xTimes"]}',
          font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text=f'{wgse.BAM.raw_avg_read_depth_NoN} {wgse.lang.i18n["xTimes"]}',
          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

#    Label(summaryFrame, text=wgse.lang.i18n['AverageReadDepthFull'], justify=LEFT,
#          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
#    Label(summaryFrame, text=f'{wgse.BAM.mapped_avg_read_depth_full} {wgse.lang.i18n["xTimes"]}',
#          font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
#    Label(summaryFrame, text=f'{wgse.BAM.raw_avg_read_depth_full} {wgse.lang.i18n["xTimes"]}',
#          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['AverageReadDepthWES'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    if WESstats:            # Summary row WSE stats ready
        Label(summaryFrame, text=f'{wgse.BAM.mapped_avg_read_depth_WES} {wgse.lang.i18n["xTimes"]}',
              font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
        Label(summaryFrame, text=f'{wgse.BAM.raw_avg_read_depth_WES} {wgse.lang.i18n["xTimes"]}',
              font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1
    else:                   # Summary row WES stats not filled in
        Button(summaryFrame, text=wgse.lang.i18n['Calculate'],
               command=lambda sw=statsWindow: button_statsWSE_BAM(sw),
               font=fonth1).grid(column=1, row=crow, columnspan=3, padx=5, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['Gbases'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=str(wgse.BAM.mapped_gbases), font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text=str(round(wgse.BAM.raw_gigabases, 2)),
          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['TotalSegs'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=f'{int(round(wgse.BAM.mapped_segs_read / (10.0 ** 6), 0))} {wgse.lang.i18n["Mill"]}',
          font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text=f'{int(round(wgse.BAM.raw_segs_read / (10.0 ** 6), 0))} {wgse.lang.i18n["Mill"]}',
          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['NumReads'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=f'{int(round(wgse.BAM.mapped_reads_percent * 100, 0))} %',
          font=fontb1).grid(row=crow, column=1, padx=0, pady=0)
    Label(summaryFrame, text="100 %",
          font=fontb1).grid(row=crow, column=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text="", font=fontb1).grid(row=crow, column=0, columnspan=4, padx=0, pady=0); crow += 1

    # Summary Frame to right, second section; Other
    if coverage:
        Label(summaryFrame, text=wgse.lang.i18n['CoverageShort'], justify=LEFT,
              font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
        Label(summaryFrame, text='{:8.4f} %'.format(round(wgse.BAM.coverage * 100, 4)),
              font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    if WESstats:
        Label(summaryFrame, text=wgse.lang.i18n['CoverageWESShort'], justify=LEFT,
              font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
        Label(summaryFrame, text='{:8.4f} %'.format(round(wgse.BAM.coverage_WES * 100, 4)),
              font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['RefModel'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.BAM.refgenome_str(),
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['AverageReadLength'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=f'{wgse.BAM.raw_avg_read_length} {wgse.lang.i18n["bp"]}, {read_type}',
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['Chroms'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.BAM.chrom_types_str(),
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['Gender'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.lang.i18n[wgse.BAM.gender] if wgse.BAM.gender else "",
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['FileStats'], justify=LEFT,
          font=fonth1).grid(row=crow, column=0, padx=0, pady=0, sticky=W)
    Label(summaryFrame, text=wgse.BAM.filestatus_str(),
          font=fontb1).grid(row=crow, column=1, columnspan=3, padx=0, pady=0); crow += 1

    Button(summaryFrame, text=wgse.lang.i18n['SaveWindow'], font=fonth1,
           command=lambda:button_save(statsWindow,"stats")).grid(row=crow, column=0, padx=5, pady=2)
    Button(summaryFrame, text=wgse.lang.i18n['CloseWindow'], font=fonth1,
           command=lambda:button_close(statsWindow)).grid(row=crow, column=1, columnspan=3, padx=5, pady=2); crow += 1

    Label(summaryFrame, text=wgse.lang.i18n['StatsNote'], font=fontb3).grid(row=crow, column=0, columnspan=4, padx=0, pady=0); crow += 1

    Label(summaryFrame, text=f'{time.ctime()};  WGSE {wgse.__version__}', justify=CENTER,
          font=fontb3).grid(row=crow, column=0, columnspan=4, padx=0, pady=0); crow += 1

    # Display and wait for Close button
    statsWindow.update()
    statsWindow.grab_set()          # Toplevel equivalent of mainloop
    statsWindow.wait_window()


def button_statsCov_BAM(window):
    """
    A separate button just to get the coverage per SN.  Takes another 30 minutes so let user choose to add the column.
    Initiated by hitting button in Coverage column of main stats page.  Will replace that page with a new one once read.
    IMPORTANT: Takes the covbases column 5 only.  Uses the stored model LN and lookup N to get a better model length
    to determine the actual coverage from.  Makes a big difference on the Y chromosome!
    """
    wgse.BAM.get_coverage_stats(window, button_directly=True)      # Call wgse.BAM routine to run samtools coverage command

    try:
        window.destroy()                # Remove the current stats window (left up while calculating the new one)
    except:
        pass

    button_stats_BAM()     # Call the button to regenerate the stats window but with the Coverage data now available


def button_statsWSE_BAM(window):
    """
    A separate button just to calculate the WES average read depth and coverage of sample. Takes another 30 minutes
    so let user choose when they want to add data.  Initiated by hitting the Calculate button next to the WES
    label / title.  Will replace current stats page with new one once done and fill in set values
    """
    wgse.BAM.get_WES_stats(window, button_directly=True)       # Call wgse.BAM routine to run samtools depth on WES region

    try:
        window.destroy()                # Remove the current stats window (left up while calculating the new one)
    except:
        pass

    button_stats_BAM()     # Call the button to regenerate the stats window but with WSE data now available


def button_save(window, wtype):
    """ Processing user button to save BAM Stats windows """
    # Find Upper Left and Lower Right corners of Stats Window
    ulx = window.winfo_rootx();  uly = window.winfo_rooty()
    lrx = window.winfo_width();  lry = window.winfo_height()
    ImageGrab.grab().crop((ulx, uly, ulx+lrx, uly+lry)).convert("RGB").save(f'{wgse.outdir.oFPB}_{wtype}.jpg')


def button_close(window, top=True):
    """ Processing user button to close sub-window like BAM stats or some result """
    try:
        window.destroy()
    except:
        pass
    if top:
        set_BAM_window_settings()  # May have updated BAM stats by clicking the stats button, so ...
        resume_main_window()


def button_header_BAM():
    """ Display BAM/CRAM/SAM header in a pop-up results window and write text file """

    header_oFN = f'{wgse.outdir.oFP}{wgse.BAM.file_FB}_header.txt'
    with open(header_oFN, "wb") as f:          # 15 Mar 2020 REH Changed to binary write to keep Unix \n format
        f.write(wgse.BAM.Header.encode())

    # Todo need scrollbars on simple_result_window for header
    simple_result_window(wgse.BAM.file_FBS + " Header", wgse.BAM.Header, "header", True)


def button_sort_BAM():
    """
        Recreate BAM in coordinate sorted format. This is an exception where generated data is put where the current
        BAM file is; not in the output or temp directory.  We will use the sorted BAM going forward.
    """
    samtools = wgse.samtoolsx_qFN
    bamfile  = wgse.BAM.file_qFN
    tempdir  = f'"{wgse.tempf.FP}"'

    # Try to avoid adding conflicting names
    out_FPB =  wgse.BAM.file_FPB.replace("_unsorted", "_sorted") if "_unsorted" in wgse.BAM.file_FPB else \
            f'{wgse.BAM.file_FPB}_sorted'

    # CRAM decode requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    suffix = ".cram" if wgse.BAM.file_type == "CRAM" else ".bam"
    out_qFN = f'"{out_FPB}{suffix}"'

    if wgse.BAM.file_type == "CRAM" and not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return      # Check routine reports error if reference does not exist

    # Samtools sort cannot accept a CRAM and have a reference genome specified so ...
    # Although not clear if --reference works on sort command; does not complain but does not seem to help
    #   samtools view | sort on CRAM is 374min (67min real) on i9 - 9900 with 40GB DRAM and NVMe SSD
    #   samtools sort --reference on CRAM directly is 434m (119min read)
    #   sort on BAM is 370min (45min real) (BAM made from CRAM earlier)
    commands = (
        f'{samtools} view -uh --no-PG {cram_opt} {bamfile} | '
        f'  {samtools} sort -T {tempdir} -m 2G -@ {wgse.os_threads} -o {out_qFN} \n'
    )

    run_bash_script("GenSortedBAM", commands)

    button_select_BAM_file(unquote(out_qFN))     # Replace selected BAM and garbage collect old one

    resume_main_window()


def button_index_BAM():
    """
        Create BAM Index file. This is an exception where generated data is put where the current
        BAM file is; not in the output or temp directory.
    """
    # This is an exception where generated data is put where BAM file is; not output or temp directory
    samtools = wgse.samtoolsx_qFN
    bamfile  = wgse.BAM.file_qFN
    commands = f'{samtools} index {bamfile}\n'

    run_bash_script("GenBAMIndex", commands)

    if wgse.BAM.check_for_bam_index():  # Check if successful
        wgse.BAM.Indexed = True
        set_BAM_window_settings()  # If index became available, change its buttons state
    else:
        wgse_message("warning", 'InvalidBAMTitle', False, 'InfoBamIndexError')

    resume_main_window()


def button_CRAM_to_BAM():

    # CRAM decode requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    if wgse.BAM.file_type == "CRAM" and not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return      # Check routine reports error if reference does not exist

    CRAM_qFN = f'"{wgse.BAM.file_FN}"'
    BAM_qFN = f'"{wgse.outdir.FPB}.bam"'
    DEBUG(f"Output BAM (final): {BAM_qFN}")

    samtools = wgse.samtoolsx_qFN       # Know we are using 64 bit samtools 1.10; 64 bit required for CRAM decode
    commands  = (
        f'{samtools} view -bh {cram_opt} -@ {wgse.os_threads} -o {BAM_qFN} {CRAM_qFN} \n'
        f'{samtools} index {BAM_qFN} \n'
    )

    run_bash_script("CRAMtoBAM", commands)

    button_select_BAM_file(unquote(BAM_qFN))    # Replace selected BAM and garbage collect old one

    resume_main_window()


def button_BAM_to_CRAM():

    # CRAM encode requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}'
    if not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return      # Check routine reports error if reference does not exist

    BAM_qFN = f'"{wgse.BAM.file_FN}"'
    CRAM_qFN = f'"{wgse.outdir.FPB}.cram"'
    DEBUG(f"Output CRAM (final): {CRAM_qFN}")

    samtools = wgse.samtoolsx_qFN  # Know we are using 64 bit samtools 1.10; 64 bit required for CRAM encode
    commands = (
        f'{samtools} view -Ch {cram_opt} -@ {wgse.os_threads} -o {CRAM_qFN} {BAM_qFN} \n'
        f'{samtools} index {CRAM_qFN} \n'
    )

    run_bash_script("BAMtoCRAM", commands)

    button_select_BAM_file(unquote(CRAM_qFN))   # Replace selected CRAM and garbage collect old one

    resume_main_window()


def button_realign_BAM(paired_BAM=True):
    """
        Button to ask to realign BAM/CRAM/SAM from one reference model to another. Involves doing an unalign
        back to FASTQ files and then alignment to a new BAM/CRAM/SAM after automatically selecting a "brother"
        reference model.  This button is fixed to do to/from Build 37/38 only of the same model.
        Will generalize later to allow specific reference model to be chosen.
    """
    cpus = wgse.os_threads
    wgse_message("warning", "RealignBAMTimeWarnTitle", True,
                 wgse.lang.i18n["RealignBAMTimeWarnMesg"].replace("{{time}}",str(5+160/cpus)))

    # Try and create new BAM based on old BAM characteristics and name. If hg19, use hg38. If hs38, use hs37d5.
    new_refgenome, new_SNCnt = wgse.reflib.refgen_liftover(wgse.BAM.Refgenome)
    new_refgenome_qFN = wgse.reflib.get_reference_genome_qFN(new_refgenome, new_SNCnt)
    if not wgse.reflib.check_reference_genome_qFN(new_refgenome_qFN):       # OK if not quoted
        resume_main_window()
        return      # Check routine reports error if reference does not exist

    # Let's find a sensible new BAM filename for realignment (paired-BAM mode; default)
    newBAM_FN = wgse.BAM.realign_BAM_filename(new_refgenome)

    # Need FASTQs (unaligned) so simply make button call. Returns FASTQs set in wgse.BAM.Rxfastq variables
    # Will check if already exist and simply return (FALSE; as did not have to make them) if found.
    made_new_fastqs = button_unalign_BAM()

    # This is an 8 hour to 6 day job; split into parts
    made_new_BAM = button_align_BAM(wgse.BAM.R1fastq_FN, wgse.BAM.R2fastq_FN, new_refgenome_qFN, newBAM_FN)

    if made_new_BAM:
        button_select_BAM_file(newBAM_FN)  # Replace newly created BAM/CRAM and garbage collect old one

        # If BAM successfully replaced and made new FASTQs as part of it; then delete FASTQs (on hold for now)
        if wgse.BAM.file_FN == newBAM_FN and made_new_fastqs:
        #    wgse.tempf.list.append(nativeOS(wgse.BAM.R1fastq_FN));  wgse.BAM.R1fastq_FN = None
        #    wgse.tempf.list.append(nativeOS(wgse.BAM.R2fastq_FN));  wgse.BAM.R2fastq_FN = None
            pass

    resume_main_window()
    # Todo need paramter for when internal call to know when to refresh wgse.BAM and window


def button_align_BAM(f1_FN, f2_FN, refgen_qFN, newBAM_FN):
    """
        Button to align (create) BAM/CRAM/SAM from a reference model and FASTQs.
        If no parameters (likely for button on GUI) then generate pop-up asking for all four values

        Major Steps (each with own PleaseWait):
          (a) Index Reference Genome for alignment (if not already) (1+ hour, ~5GB)
          (b) Call ALignment tool on FASTQs to create initial raw BAM (8+ hours, ~50GB)
          (c) Cleanup (fixmate, markdups, sort) and Index new BAM (1 hour, same size)
          (d) If needed, convert to CRAM and remove BAM (1 hour)

        Inspirations:
        https://eriqande.github.io/eca-bioinf-handbook/alignment-of-sequence-data-to-a-reference-genome-and-associated-steps.html
        http://www.htslib.org/workflow/ (first section on FASTQ to BAM/CRAM)
        https://gist.github.com/tkrahn/7dfc51c2bb97a6d654378a21ea0a96d4 (although some he copied from us :)
    """
    # Todo should content be moved to bamfile.py?

    #if wgse.os_plat == "Windows" and not wgse.wsl_bwa_patch:
    #    wgse_message("warning", "AlignBAMWin10WarnTitle", False, "AlignBAMWin10WarnMesg")

    # Todo develop get_align_params() routine to query values from user if not specified (raw align from scratch)
    #if newBAM_FN or f1 or f2 or refgen:     # Something not specified, need to query for input
    #    (newBAM_FN, f1, f2, refgen) = get_align_params(newBAM_FN, f1, f2, refgen)

    bgzip = wgse.bgzipx_qFN
    samtools = wgse.samtoolsx_qFN
    cpus = wgse.os_threads
    bwa = 'wsl bwa' if wgse.wsl_bwa_patch else wgse.bwax_qFN

    # Ugly. We need SN count to an unknown ref genome. Luckily count only needed for errored / unlikely ref genome files
    # refgen_qFN = wgse.reflib.get_reference_genome_qFN(refgen, SNCnt)
    refgen_oFN = nativeOS(unquote(refgen_qFN))      # Remove quotes and change to native OS
    DEBUG(f'Reference Genome File: {refgen_oFN}')
    if not wgse.reflib.check_reference_genome_qFN(refgen_qFN):
        return False
    refgen_quFN = f'"{universalOS(unquote(refgen_qFN), wgse.wsl_bwa_patch if wgse.os_plat == "Windows" else False)}"'

    newBAM_FPB, newBAM_FS = os.path.splitext(newBAM_FN)
    newBAM_FB = os.path.basename(newBAM_FPB)
    DEBUG(f'New BAM File Basename: {newBAM_FB}, extension: {newBAM_FS}')
    newBAM_qFN = f'"{newBAM_FN}"'

    # Define a trueBAM in case newBAM is asking for a CRAM
    trueBAM_FN = newBAM_FN.replace(".cram", ".bam") if newBAM_FS == ".cram" else newBAM_FN
    trueBAM_qFN = f'"{trueBAM_FN}"'

    markdup_result_FN = f'"{newBAM_FPB}_markdup.txt"'       # Capture MarkDup report

    # Setup some intermediate files
    temp_rawalign_FN  = f'{wgse.tempf.FP}{newBAM_FB}_raw.bam'
    temp_rawalign_oFN = nativeOS(temp_rawalign_FN)
    temp_rawalign_qFN = f'"{temp_rawalign_FN}"'
    temp_fixmate_qFN  = f'"{wgse.tempf.FP}{newBAM_FB}_fixmate.bam"'
    tempdir_qFN  = f'"{wgse.tempf.FP}"'

    f1_oFN = nativeOS(f1_FN)
    f2_oFN = nativeOS(f2_FN)
    f1_qFN = f'"{f1_FN}"'
    f2_qFN = f'"{f2_FN}"'
    f1_quFN = f'"{universalOS(f1_FN, wgse.wsl_bwa_patch if wgse.os_plat == "Windows" else False)}"'
    f2_quFN = f'"{universalOS(f2_FN, wgse.wsl_bwa_patch if wgse.os_plat == "Windows" else False)}"'

    #------------------------------------------------------------------------------------------
    # (a) Check for BWA Index files -- generate if missing
    if not (os.path.isfile(refgen_oFN) and os.path.getsize(refgen_oFN) > 500000000):
        DEBUG(f'Error: missing Reference Genome File: {refgen_qFN}')
        return False
    if os.path.isfile(refgen_oFN + ".bwt") and \
       os.path.getmtime(refgen_oFN + ".bwt") > os.path.getmtime(refgen_oFN):
        DEBUG(f'Using previous BWA Indices for RefGenome: {refgen_qFN}')
    else:  # Need to index the Reference Genome
        commands = f'{bwa} index {refgen_quFN} \n'
        run_bash_script("CreateAlignIndices", commands)

    #------------------------------------------------------------------------------------------
    # (b) Lets do the alignment ... 8 hours to 6 days likely
    # todo add picard and GATK markdup function in place of samtools
    if not (os.path.isfile(f1_oFN) and os.path.isfile(f2_oFN)):
        DEBUG(f'Error: missing FASTQ files: {f1_qFN}, {f2_qFN}')
        return False
    fastq_size = os.path.getsize(f1_oFN) + os.path.getsize(f2_oFN)
    if os.path.isfile(temp_rawalign_oFN) and \
       os.path.getmtime(temp_rawalign_oFN) > os.path.getmtime(f1_oFN) and \
       os.path.getsize(temp_rawalign_oFN) / fastq_size > 0.8 and \
       os.path.getsize(temp_rawalign_oFN) / fastq_size < 1.2:
        DEBUG(f'Found previous RAW alignment to reuse: {temp_rawalign_qFN}')  # Not likely but worth a try
    else:
        # Run the alignment command; compress the output simply because it is so large otherwise
        commands  = f'{bwa} mem -t {cpus} {refgen_quFN} {f1_quFN} {f2_quFN} | {bgzip} -@ {cpus} > {temp_rawalign_qFN}\n'
        run_bash_script("ButtonAlignBAM", commands)

    #------------------------------------------------------------------------------------------
    # (c) Cleanup (fixmate, markdups, sort) and Index new BAM
    if os.path.isfile(nativeOS(trueBAM_FN)) and \
       os.path.getmtime(trueBAM_FN) > os.path.getmtime(temp_rawalign_oFN):
        DEBUG("Found previous cleaned-up and aligned BAM to reuse: {trueBAM_FN}")
    else:
        commands = (
          f'{samtools} view -uh --no-PG {temp_rawalign_qFN} | '
          f'  {samtools} fixmate -m -O bam -@ {cpus} - {temp_fixmate_qFN} \n'
          f'{samtools} sort -T {tempdir_qFN} -m 2G -@ {cpus} {temp_fixmate_qFN} | '
          f'  {samtools} markdup -s -d 2500 -@ {cpus} - {trueBAM_qFN} > {markdup_result_FN}\n'
          f'{samtools} index {trueBAM_qFN} \n'
        )
        run_bash_script("AlignCleanup", commands)
        # Note that a coordinate sorted BAM can be considerably smaller than a name or unsorted BAM. Compression
        # is more effective when the similar sequences are brought close to each other.

    #------------------------------------------------------------------------------------------
    # (d) Convert to CRAM and remove BAM (if CRAM is actually being requested)
    if newBAM_FS == ".cram":
        # Note: would not find existing CRAM that is newer as we remove old BAM when converted
        # todo modify CRAM_to_BAM() to accept parameters and simply call that here
        commands = (
            f'{samtools} view -Ch -T {refgen_qFN} -@ {cpus} -o {newBAM_qFN} {trueBAM_qFN} \n'
            f'{samtools} index {newBAM_qFN} \n'
        )
        run_bash_script("BAMtoCRAM", commands)

        wgse.tempf.list.append(nativeOS(trueBAM_FN))           # Remove intermediate BAM
        wgse.tempf.list.append(nativeOS(trueBAM_FN + ".bai"))  #  and its index

    # Internal call only currently; so do not replace BAM or resume main window
    #button_select_BAM_file(newBAM_FN)  # Replace selected BAM/CRAM and garbage collect old one
    #resume_main_window()
    # Todo need paramter for when internal call to know whether to refresh wgse.BAM and window

    return True


def button_unalign_BAM():
    """
    Function to go from BAM / CRAM to FASTQ(s). Assumes BAM already subsetted if only want a subset FASTQ.
    Looks for result in output area first before acting
    """

    # Todo handle single-end BAM / CRAM files that will generate a single fasta file

    # Let's check if FASTQ's already exist; BAMfile class can do it and store locally if found
    if wgse.BAM.find_FASTQs():
        return False    #    Meaning no need to create; already exist. Simple return as set in BAM class already

    r1fastq = f'"{wgse.outdir.FPB}_R1.fastq.gz"'    # File names to create here
    r2fastq = f'"{wgse.outdir.FPB}_R2.fastq.gz"'

    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    samtools = wgse.samtoolsx_qFN
    bamfile  = wgse.BAM.file_qFN
    tempdir = f'"{wgse.tempf.FP}"'

    if wgse.BAM.file_type == "CRAM" and not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return False      # Check routine reports error if reference does not exist

    # Sort in name order, then call fastq command to split and write FastQ's
    # Samtools sort cannot take the reference genome specification so have to view CRAM first
    commands = (
        f'{samtools} view -uh --no-PG {cram_opt} {bamfile} | '
        f'  {samtools} sort -n -T {tempdir} -m 2G -@ {wgse.os_threads} -O sam | '
        f'  {samtools} fastq -1 {r1fastq} -2 {r2fastq} -0 /dev/null -s /dev/null -n -@ {wgse.os_threads} \n'
    )

    run_bash_script('ButtonUnalignBAM', commands)

    wgse.BAM.R1fastq_FN = unquote(r1fastq)
    wgse.BAM.R2fastq_FN = unquote(r2fastq)

    resume_main_window()
    # Todo need paramter for when internal call to know whether to refresh wgse.BAM and window
    return True


def button_select_output_path(out_FP=None):
    """ Processing user button to specify Output Directory path """

    if out_FP:
        new_FP = out_FP
    else:       # If not supplied a path, query the user for one
        # Get directory to start in from the one containing the current set directory; if set. Otherwise, let OS decide
        initialdir = os.path.dirname(wgse.outdir.oFP[:-1]) if wgse.outdir and wgse.outdir.oFP else ""
        new_FP = filedialog.askdirectory(parent=wgse.window, title=wgse.lang.i18n['SelectOutputDirectory'],
                                         initialdir=initialdir)      # Returns Unix/Universal, and no trailing slash
        if not new_FP:      # If still not set ...
            DEBUG(f"Output Path not set (cancelled)")
            resume_main_window()
            return      # User likely hit exit or cancel; simply return (do not clear old value)

    new_FP += '/' if new_FP[-1] != '/' else new_FP  # Assure trailing os_slash; universalOS so always forward slash

    wgse.outdir.change(new_FP)

    if wgse.outdir.FB:
        outputDirectoryButton.configure(text=wgse.outdir.FB)
    else:
        wgse_message("error", 'InvalidPathWindowTitle', False, 'errOutputPathSpecialChars')

    wgse.save_settings()
    resume_main_window()


def button_mtdna_BAM():
    """
    Processing user button to generate mito only BAM (always BAM, not CRAM)
    Put in output directory; does not replace currently selected BAM (extract button; not BAM settings button)
    """

    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else  ""
    if wgse.BAM.file_type == "CRAM" and not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return False     # Check routine reports error if reference does not exist

    samtools = wgse.samtoolsx_qFN
    outMTbam = f'"{wgse.outdir.FPB}_chrM.bam"'

    chromos = wgse.BAM.get_chr_name('M')

    commands = f'{samtools} view -bh {cram_opt} {wgse.BAM.file_qFN} {chromos} -o {outMTbam} \n'

    run_bash_script("ButtonMitoBAM", commands)

    if not os.path.isfile(nativeOS(unquote(outMTbam))):
        wgse_message("error", 'errYTitle', False, 'errYBAMFile')

    resume_main_window()


def button_mtdna_VCF(internal=False):
    """
    Processing user button to generate the MITO VCF file (also from internal call to do the same).
    VCF is used by both mitoFASTA button and Haplogroup (haplogrep) call.  So reuse if available.
    Always place and leave in output directory for possible reuse.
    """
    if wgse.BAM.RefMito == "Yoruba" and not internal:
        wgse_message("warning", 'YorubaWindowTitle', False, 'YorubaWarning')

    chrM_qFN = f'"{wgse.outdir.FPB}_chrM.vcf.gz"'
    chrM_oFN = nativeOS(unquote(chrM_qFN))

    if not (os.path.isfile(chrM_oFN) and os.path.getmtime(chrM_oFN) > os.path.getmtime(wgse.BAM.file_oFN)):
        refgenome_qFN = wgse.BAM.Refgenome_qFN
        if not wgse.reflib.check_reference_genome_qFN(refgenome_qFN):
            if not internal:
                resume_main_window()
            return False    # Check routine reports error if reference does not exist

        bcftools = f'{wgse.bcftoolsx_qFN}'
        tabix = f'{wgse.tabixx_qFN}'
        cpus = wgse.os_threads
        chrM = wgse.BAM.get_chr_name('M')
        ploidy = "1"  # mito-only

        #   Note -B required for Nanopore long read tools; minimal effect on paired-end, mpp sequencer output
        commands = (
            f'{bcftools} mpileup -B -I -C 50 -r {chrM} -f {refgenome_qFN} -Ou {wgse.BAM.file_qFN} | '
            f'  {bcftools} call --ploidy {ploidy} -mv -P 0 --threads {cpus} -Oz -o {chrM_qFN} \n'
            f'{tabix} {chrM_qFN}\n'
        )

        run_bash_script("ButtonMitoVCF", commands)

        if not os.path.isfile(chrM_oFN):
            wgse_message("error", 'errMtTitle', False, 'errMtVCFFile')
            if not internal:
                resume_main_window()
            return False

    if not internal:
        resume_main_window()
    return True


def button_mtdna_FASTA():
    """
    Processing user button to generate the MITO FASTA file (also from internal call to do the same).
    Note: The FASTA in the industry is more like a FA and not a FASTQ. Not raw reads from the sequencer.
    But single ~16k base-pair segment of all calls. So like for Autosomal, have to mpileup and call.
    """
    if wgse.BAM.RefMito == "Yoruba":
        wgse_message("warning", 'YorubaWindowTitle', False, 'YorubaWarning')

    fasta_qFN = f'"{wgse.outdir.FPB}_mtdna.fasta"'
    fasta_oFN = nativeOS(unquote(fasta_qFN))
    if not (os.path.isfile(fasta_oFN) and os.path.getmtime(fasta_oFN) > os.path.getmtime(wgse.BAM.file_oFN)):
        chrM_qFN = f'"{wgse.outdir.FPB}_chrM.vcf.gz"'
        if not button_mtdna_VCF(internal=True):
            resume_main_window()
            return      # Generate button reports error message if does not exist

        refgenome_qFN = wgse.BAM.Refgenome_qFN
        if not wgse.reflib.check_reference_genome_qFN(refgenome_qFN):
            resume_main_window()
            return      # Check routine reports error if reference does not exist

        samtools = f'{wgse.samtoolsx_qFN}'
        bcftools = f'{wgse.bcftoolsx_qFN}'
        chrM = wgse.BAM.get_chr_name('M')

        # -B required for Nanopore long read tools; minimal effect on paired-end, mpp sequencer output
        commands = f'{samtools} faidx {refgenome_qFN} {chrM} | {bcftools} consensus {chrM_qFN} -o {fasta_qFN} \n'

        run_bash_script("ButtonMitoFASTA", commands)

        if not os.path.isfile(fasta_oFN):
            wgse_message("error", 'errMtTitle', False, 'errMtFastaFile')

    resume_main_window()


def button_mtdna_haplogroup():
    """ Process user button to call haplogrep to calculate mtDNA haplogroup. """
    if wgse.BAM.RefMito == "Yoruba":
        wgse_message("error", 'YorubaWindowTitle', False, 'YorubaWarning')
        resume_main_window()
        return

    if wgse.javax_FN is None:         # Really still needed? Is not checked in Settings at startup?
        wgse_message("error", 'MissingJavaTitle', False, 'MissingJava')
        resume_main_window()
        return

    chrM_qFN = f'"{wgse.outdir.FPB}_chrM.vcf.gz"'
    if not button_mtdna_VCF(internal=True):
        resume_main_window()
        return  # Generate button reports error message if does not exist

    haplogrep = f'{wgse.javax_FNp} "{wgse.jartools_FP}haplogrep.jar"'
    haplogroup_qFN = f'"{wgse.tempf.FP}mtdna_haplogroup.txt"'

    commands = f'{haplogrep} classify --in {chrM_qFN} --format vcf --out {haplogroup_qFN}\n'

    run_bash_script('ButtonMTHaplo', commands)

    haplogroup_oFN = nativeOS(unquote(haplogroup_qFN))
    if not os.path.isfile(haplogroup_oFN):
        wgse_message("error", 'errMtTitle', False, 'errMtHaploFile')
        resume_main_window()
        return

    mtdna_haplogroup = ""
    with open(haplogroup_oFN, "r") as source_content:
        for source_line in source_content:
            line_tabs = source_line.split("\t")
            if str(line_tabs[2]).strip('"') == "1":
                mtdna_haplogroup = line_tabs[1].strip('"')
                # Todo can we do a break here or might there be multiple with the third column equal to 1?
    report = wgse.lang.i18n['MitoHaploReport']\
        .replace("{{mthaplogroup}}", mtdna_haplogroup).replace("{{bamf}}", wgse.BAM.disp_FBS)
    simple_result_window(wgse.lang.i18n['MitochondrialDNA'], report, "mHaplo", True)

    resume_main_window()


def button_yAndMt_BAM():
    """ Processing user button to generate Y and mtDNA only BAM (e.g. for yFull) """
    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    if wgse.BAM.file_type == "CRAM" and not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return     # Check routine reports error if reference does not exist

    samtools = wgse.samtoolsx_qFN

    # Important to retain sort order with header for further processing in the future
    chromos = f'{wgse.BAM.get_chr_name("Y")} {wgse.BAM.get_chr_name("M")}'
    if wgse.BAM.Build == 38 and wgse.BAM.SNTypeC == "Chr":
        chromos += f' chrY_KI270740v1_random'       # Only Y option; added only to Build 38 based on HG model
    out_qFN = f'"{wgse.outdir.FPB}_chrYM.bam"'

    commands = f'{samtools} view -bh {cram_opt} {wgse.BAM.file_qFN} {chromos} -o {out_qFN}\n'

    run_bash_script("ButtonYandMT", commands)

    if not os.path.isfile(nativeOS(unquote(out_qFN))):
        wgse_message("error", 'errYTitle', False, 'errYandMtFile')

    resume_main_window()


def button_yOnly_BAM(internal=False):
    """
    Processing user button to generate Y only BAM (always BAM, not CRAM)
    Put in output directory; does not replace currently selected BAM (extract button; not BAM settings button)
    Can be called internally. If so, put in temp area.
    """

    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else  ""
    if wgse.BAM.file_type == "CRAM" and not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return False     # Check routine reports error if reference does not exist

    samtools = wgse.samtoolsx_qFN
    outybam = f'"{wgse.outdir.FPB}_chrY.bam"' if not internal else \
              f'"{wgse.tempf.FP}{wgse.BAM.file_FB}_chrY.bam"'

    chromos = wgse.BAM.get_chr_name('Y')
    if wgse.BAM.Build == 38 and wgse.BAM.SNTypeC == "Chr":
        chromos += f' chrY_KI270740v1_random'       # Only Y option; added only to Build 38 based on HG model

    commands = f'{samtools} view -bh {cram_opt} {wgse.BAM.file_qFN} {chromos} -o {outybam} \n'

    run_bash_script("ButtonYonly", commands)

    if not internal:
        resume_main_window()
    if not os.path.isfile(nativeOS(unquote(outybam))):
        wgse_message("error", 'errYTitle', False, 'errYBAMFile')
        return False

    return True


def button_yOnly_VCF():
    """
    Special button to create subset of VCF for Y only; used in Cladefinder.yseq.net.
    We keep the final VCF uncompressed as that is what Cladefinder can deal with.
    """

    bcftools = wgse.bcftoolsx_qFN
    tabix = wgse.tabixx_qFN
    cpus = wgse.os_threads

    refVCFtab_qFN = wgse.reflib.get_reference_vcf_qFN(wgse.BAM.Build, wgse.BAM.SNTypeC, type="Yonly")
    refgenome_qFN = wgse.BAM.Refgenome_qFN
    if not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return      # Check routine reports error if reference does not exist

    # Subset to Y-only BAM if not already available; big time saver when creating called variances
    bamfile_FN = f'{wgse.outdir.FP}{wgse.BAM.file_FB}_chrY.bam'
    bamfile_oFN = nativeOS(bamfile_FN)
    if wgse.BAM.chrom_types['A'] < 1 and wgse.BAM.chrom_types['Y'] > 1:
        bamfile_FN = wgse.BAM.file_FN   # Already a subsetted BAM; just use directly
    elif not (os.path.exists(bamfile_oFN) and os.path.getmtime(bamfile_oFN) > os.path.getmtime(wgse.BAM.file_oFN)):
        # Have to create a sub-setted y-only BAM; use Temp since buried in this routine
        bamfile_FN = f'{wgse.tempf.FP}{wgse.BAM.file_FB}_chrY.bam'   # As created by button_yOnly() in internal mode
        bamfile_oFN = nativeOS(bamfile_FN)
        if not (os.path.exists(bamfile_oFN) and os.path.getmtime(bamfile_oFN) > os.path.getmtime(wgse.BAM.file_oFN)):
            # Just check to see if already in temp directory first; otherwise create
            button_yOnly_BAM(internal=True)       # New BAM is in temp due to internal mode call

    # So use bamfile name chosen (without extension) to create various intermediate temp and final file names
    bamfile_FB = os.path.basename(os.path.splitext(bamfile_FN)[0])
    temp_called_bcf = f'"{wgse.tempf.FP}{bamfile_FB}.called.bcf.gz"'
    annotated_vcf = f'"{wgse.outdir.FP}{bamfile_FB}.annotated.vcf.gz"'      # Note: not in temp area and not compressed
    bamfile_qFN = f'"{bamfile_FN}"'

    ploidy = "1"    # Y only; not enabled for female so ...

    """  Trying to include and left-justify (normalize) InDels but did not get it working quite right
    temp_norm_bcf = f'"{wgse.tempf.FP}{bamfile_FB}.norm.bcf.gz"'
    f'{bcftools} norm -f {refgenome_qFN} {temp_called_bcf} -O b -o {temp_norm_bcf}\n'
    f'{bcftools} index {temp_norm_bcf}\n'
    """
    #   Note -B required for Nanopore long read tools; minimal effect on paired-end, mpp sequencer output
    commands = (
        f'{bcftools} mpileup -B -I -C 50 -T {refVCFtab_qFN} -f {refgenome_qFN} -Ou {bamfile_qFN} | '
        f'  {bcftools} call --ploidy {ploidy} -V indels -m -P 0 --threads {cpus} -Ob -o {temp_called_bcf}\n'
        f'{bcftools} index {temp_called_bcf}\n'
        f'{bcftools} annotate -a {refVCFtab_qFN} -c CHROM,POS,ID {temp_called_bcf} --threads {cpus} -Oz -o {annotated_vcf}\n'
        f'{tabix} -p vcf {annotated_vcf}\n'  # cladefinder cannot deal with BCF
    )
    run_bash_script("AnnotatedVCF-yOnly", commands)

    if not os.path.exists(nativeOS(unquote(annotated_vcf))):
        wgse_message("error", 'errYTitle', False, 'errYVCFFile')
    resume_main_window()


def button_ydna_haplogroup():
    """ Processing user button to run yleaf to generate Y Haplogroup call. """

    outd_FP = f'{wgse.tempf.FP}tempYleaf/'
    outf = f'"{outd_FP}haplogroups.txt"'        # note: _qFN
    outd_oFP = nativeOS(outd_FP)
    if os.path.isdir(outd_oFP):       # REH 10Mar2020 If previous run crashed ...
        wgse.tempf.list.append(outd_oFP)    # Adding it to list assures deletion even if DEBUG_MODE is True
        wgse.tempf.clean()

    position_FBS = "WGS_hg38.txt" if "38" in wgse.BAM.Refgenome else "WGS_hg19.txt"
    if wgse.BAM.file_type == "CRAM":
        filespec = f"-f {wgse.BAM.Refgenome_qFN} -cram {wgse.BAM.file_qFN}"
        if not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
            resume_main_window()
            return  # Check routine reports error if reference does not exist
    else:  # Is BAM or SAM; the same for the tool
        filespec = f"-bam {wgse.BAM.file_qFN}"

    # Each file path needs to be surrounded by quotes       # REH 10Mar2020
    ylef = f'"{wgse.yleaf_FP}yleaf.py"'
    posf = f'"{wgse.yleaf_FP}Position_files/{position_FBS}"'
    pyth = f'"{wgse.python3x_FN}"'
    samt = wgse.samtoolsx_qFN
    # Todo REH 10Mar2020 bug from wgsev2?  -r specified twice: -r 3 and -r 1?
    # Modified yleaf to accept python and samtools exec path;

    commands = f'{pyth} {ylef} {filespec} -pos {posf} -out {outf} -r 1 -q 20 -b 90 -py {pyth} -samt {samt}\n'

    run_bash_script('ButtonYHaplo', commands)

    result_oFN = nativeOS(unquote(outf))
    if not os.path.exists(result_oFN):
        print('--- FAILURE in yLeaf call (no result file)!')
        resume_main_window()
        return

    with open(result_oFN) as fp:    # MB Jan2020
        fp.readline()               # REH 13Mar2020 This would work if file named correctly; skip header
        fpdata = fp.readline()      # REH 13Mar2020 Read only first result; only supplying one file
        columns_yleaf_result = fpdata.split("\t")
        yhg = columns_yleaf_result[1]   # Second column is Y Haplogroup (in YCC Long form)

    if yhg == 'NA':
        print('--- FAILURE in yLeaf call (NA result)!')
        resume_main_window()
        return

    out_FB  = wgse.BAM.file_FB
    grep    = wgse.grepx_qFN
    awk     = wgse.awkx_qFN
    awkargs = "-F \"\\\"*\\t\\\"*\" '{print $3}'"   # Todo simplify with raw strong r""?
    outsuf  = f'"{outd_oFP}{out_FB}{wgse.os_slash}{out_FB}.out"'
    termsnp = f'"{outd_oFP}snps_terminal.txt"'
    # Grep through output listing of 64,000+ SNPs and find all with same haplogroup as terminal
    commands = f'{grep} "{yhg}" {outsuf} | {awk} {awkargs} > {termsnp}\n'

    run_bash_script('ButtonYHaplo2', commands)

    termsnp_oFN = f'{outd_oFP}snps_terminal.txt'    # Previous one is quoted, not-native
    if not os.path.exists(termsnp_oFN):
        print('--- Failure in yLeaf post-processing (no terminal SNPs file)!')
        resume_main_window()
        return
    with open(termsnp_oFN) as f:
        terminal_snps = [line.rstrip('\n') for line in f]

    yHaplo_result_window(yhg, terminal_snps)


def yHaplo_result_window(yhg, terminal_snps):
    """ Y Chromosome Haplogroup call results window pop-up. """
    global yHaploResWindow

    total_snp_count = len(terminal_snps)
    if total_snp_count > 0:
        full_snp_list = ", ".join(terminal_snps)
        snps_textline = ", ".join(terminal_snps[0:total_snp_count if total_snp_count < 3 else 3])
        if total_snp_count > 3:
            snps_textline += " ..."
    else:
        snps_textline = ""  # 1 to 3 SNPs to report directly
        full_snp_list = ""  # The full SNPs list (different than above if more than 3)

    snps_label  = wgse.lang.i18n['SNPsForThisHG'].replace('{{yhg}}', yhg).replace('{{snps}}', snps_textline)
    yFull_label = wgse.lang.i18n['findHgYFullDescription'].replace('{{BspSNP}}', terminal_snps[0])
    ftdna_label = wgse.lang.i18n['findOnFTDNADescLabel'].replace('{{BspSNP}}', terminal_snps[0]).\
        replace('{{url}}', f'{wgse.ftdna_pubytree_url}{yhg[0]}')

    # wgse.window.withdraw()        # Maybe hide main window while displaying stats? Turned off for now
    yHaploResWindow = Toplevel(wgse.window)
    yHaploResWindow.transient()
    yHaploResWindow.title(wgse.lang.i18n['YChromHaplogroup'])
    yHaploResWindow.geometry(wgse.yHgResult_winsize)
    yHaploResWindow.maxsize(width=wgse.yHgResult_maxw, height=wgse.yHgResult_maxh)
    yHaploResWindow.protocol("WM_DELETE_WINDOW", lambda win=yHaploResWindow: button_close(win))
    yHaploResWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''

    yHaploResWindow.rowconfigure(1, weight=1)
    yHaploResWindow.columnconfigure(0, weight=1)
    # On Mac and Linux, window is not sized yet
    wraplength1 = (wgse.yHgResult_maxw - 20) if yHaploResWindow.winfo_width() < 300 else \
        (yHaploResWindow.winfo_width() - 20)
    wraplength2 = int(wraplength1 * 2 / 3)

    # Upper left Haplogroup box
    yhgresFrame = LabelFrame(yHaploResWindow, text=wgse.lang.i18n['Haplogroup'], font=font['16b'])
    yhgresFrame.grid(row=0, column=0, padx=10, pady=5, sticky=(N, S, E, W))

    Label(yhgresFrame, text=f"File: {wgse.BAM.disp_FBS}", justify=LEFT,
          font=font['14']).grid(row=0, column=0, padx=5, pady=2)
    Label(yhgresFrame, text=wgse.lang.i18n['YDNAHaploReport'], justify=CENTER, wraplength=wraplength2,
          font=font['14']).grid(row=1, column=0, padx=5, pady=2)
    Label(yhgresFrame, text=f"{yhg}", justify=CENTER,
          font=font['14b']).grid(row=2, column=0, padx=5, pady=2)
    Label(yhgresFrame, text=f'{time.ctime()};  WGSE {wgse.__version__}', justify=CENTER,
          font=font['12']).grid(row=3, column=0, columnspan=4, padx=0, pady=0)

    # Upper right SNPs box
    snpsForHgFrame = LabelFrame(yHaploResWindow, text="SNPs", font=font['16b'])
    snpsForHgFrame.grid(row=0, column=1, padx=10, pady=5, sticky=(N, S, E, W))

    Label(snpsForHgFrame, text=snps_label, justify=LEFT, wraplength=snpsForHgFrame.winfo_width()-10,
          font=font['14']).grid(row=0, column=0, padx=5, pady=2)
    if total_snp_count > 3:
        Button(snpsForHgFrame, text=wgse.lang.i18n['showAllSnps'], justify=CENTER,
               command=lambda: simple_result_window(yhg + " SNPs", full_snp_list, "yHaplo_list", False),
               font=font['14']).grid(row=1, column=0, padx=5, pady=2)
    Label(snpsForHgFrame, text=wgse.lang.i18n['SNPsForThisHG2'], justify=CENTER,
          font=font['12']).grid(row=2, column=0, padx=5, pady=2)

    # Lower, spanning-columns Trees box
    findHgOtherTreesFrame = LabelFrame(yHaploResWindow, text=wgse.lang.i18n['findHgOtherTrees'], font=font['16b'])
    findHgOtherTreesFrame.grid(row=1, column=0, columnspan=2, padx=10, pady=5, sticky=(N, S, E, W))

    DEBUG(f"Wraplengths: window:{wraplength1}, 2/3 frame:{wraplength2}")
    crow = 0
    Label(findHgOtherTreesFrame, text=wgse.lang.i18n['findHgISOGGDescription'], justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Button(findHgOtherTreesFrame, text=wgse.lang.i18n['findOnISOGGButtonText'], justify=CENTER,
           command=lambda: webbrowser.open_new(wgse.isogg_tree_url),
           font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Label(findHgOtherTreesFrame, text=yFull_label, justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Button(findHgOtherTreesFrame, text=wgse.lang.i18n['findOnYFullButtonText'], justify=CENTER,
           command=lambda: webbrowser.open_new(wgse.yfull_searchsnp_url),
           font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow +=1
    Label(findHgOtherTreesFrame, text=ftdna_label, justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Button(findHgOtherTreesFrame, text=wgse.lang.i18n['findOnFTDNAButtonText'], justify=CENTER,
           command=lambda: webbrowser.open_new(f'{wgse.ftdna_pubytree_url}{yhg[0]}'),
           font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1
    Label(findHgOtherTreesFrame, text=wgse.lang.i18n['warningOtherTreesAreDeeper'], justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=crow, column=0, padx=5, pady=2); crow += 1

    Button(yHaploResWindow, text=wgse.lang.i18n['SaveWindow'],
           command=lambda: button_save(yHaploResWindow, "yHaplo"), justify=CENTER,
           font=font['14']).grid(row=2, column=0, padx=5, pady=2)
    Button(yHaploResWindow, text=wgse.lang.i18n['CloseWindow'],
           command=lambda: button_close(yHaploResWindow), justify=CENTER,
           font=font['14']).grid(row=2, column=1, padx=5, pady=2)

    yHaploResWindow.update()
    yHaploResWindow.grab_set()      # Toplevel equivalent of mainloop
    yHaploResWindow.wait_window()


def simple_result_window(title, simple_text, file_ext, top=False):
    """
    Generic simple result window handler.  Provides results given in simple_text and a save and exit button.
    """
    # wgse.window.withdraw()        # Maybe hide main window while displaying stats? Turned off for now
    simResWindow = Toplevel(wgse.window)
    simResWindow.transient()
    simResWindow.title(title)
    simResWindow.geometry(wgse.simResult_winsize)
    simResWindow.protocol("WM_DELETE_WINDOW", lambda: button_close(simResWindow))
    simResWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    simResWindow.maxsize(width=wgse.simResult_maxw, height=wgse.simResult_maxh)

    simResWindow.rowconfigure(1, weight=1)
    simResWindow.columnconfigure(0, weight=1)
    # On Mac and Linux, window is not sized yet
    wraplength1 = (wgse.simResult_maxw - 20) if simResWindow.winfo_width() < 300 else \
        (simResWindow.winfo_width() - 20)

    Label(simResWindow, text=f'{title}:', justify=LEFT,
          font=font['14b']).grid(row=0, column=0, columnspan=2, padx=5, pady=2)
    Label(simResWindow, text=simple_text, justify=LEFT, wraplength=wraplength1,
          font=font['14']).grid(row=1, column=0, columnspan=2, padx=5, pady=2)

    Label(simResWindow, text=f'{time.ctime()};  WGSE {wgse.__version__}', justify=CENTER,
          font=font['12']).grid(row=2, column=0, columnspan=2, padx=0, pady=0)

    Button(simResWindow, text=wgse.lang.i18n['SaveWindow'],
           command=lambda: button_save(simResWindow, file_ext), justify=CENTER,
           font=font['14']).grid(row=3, column=0, padx=5, pady=2)
    Button(simResWindow, text=wgse.lang.i18n['CloseWindow'],
           command=lambda: button_close(simResWindow, top), justify=CENTER,
           font=font['14']).grid(row=3, column=1, padx=5, pady=2)

    simResWindow.update()
    simResWindow.grab_set()     # Toplevel equivalent of mainloop
    simResWindow.wait_window()


def button_export_unmapped_reads(internal=False):
    """ Process user button to export unmapped reads from BAM. """
    # Todo give option for BAM and/or paired-end FASTQ generation (now does both)
    # Todo use BAM Read Type to determine if paired-end (two files) or single-end (one file)
    # Todo include key decoys added that should be considered unmapped (EBV, hs37d5?, etc)

    # CRAM file requires reference genome be specified with VIEW command
    cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
    if wgse.BAM.file_type == "CRAM" and not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
        resume_main_window()
        return      # Check routine reports error if reference does not exist

    samtools = wgse.samtoolsx_qFN
    # bgzip    = wgse.bgzipx_qFN
    bamfile  = wgse.BAM.file_qFN
    unmapbam     = f'"{wgse.outdir.FPB}_unmap.bam"'
    unmapR1fastq = f'"{wgse.outdir.FPB}_unmapR1.fastq.gz"'
    unmapR2fastq = f'"{wgse.outdir.FPB}_unmapR2.fastq.gz"'
    tempdir = f'"{wgse.tempf.FP}"'

    # bam2fq deprecated in 1.9; now use fastq
    # Subset BAM to unmapped, sort in name order, then call fastq command to split and write FastQ's
    commands = (
        f"{samtools} view -bh {cram_opt} -@ {wgse.os_threads} -o {unmapbam} {bamfile} '*' \n"
        f"{samtools} sort -n -T {tempdir} -m 2G -@ {wgse.os_threads} {unmapbam} | "
        f"  {samtools} fastq -1 {unmapR1fastq} -2 {unmapR2fastq} -0 /dev/null -s /dev/null -n -@ {wgse.os_threads}\n"
    )

    #wgse.window.withdraw()
    run_bash_script('ButtonUnmappedReads', commands)

    if internal:        # Meaning, internal call not from a user button hit
        wgse.tempf.list.append(f'{wgse.outdir.oFPB}_unmap.bam')
    else:               # On user button, give a result pop-up
        # Want just file names; not full path
        unmapbam = f'{wgse.BAM.file_FB}_unmap.bam'
        unmapR1fastq = f'{wgse.BAM.file_FB}_unmapR1.fastq.gz"'
        unmapR2fastq = f'{wgse.BAM.file_FB}_unmapR2.fastq.gz"'
        temp_descript = wgse.lang.i18n['DescriptionUploadCosmosId']\
            .replace("{{fq1}}", unmapR1fastq).replace("{{fq2}}", unmapR2fastq).replace("{{obam}}", unmapbam)
        simple_result_window(wgse.lang.i18n['HowToContinue'], temp_descript, "unmap", True)

    resume_main_window()


def button_extract_WES():
    """
    samtools view -bh mybam*.bam -l TruSeq_Exome*bed -b -o mybam_WES.bam
    samtools index mybam_WES.bam
    BED file is chromosome name dependent; one here from Illumina is HG
    """
    wgse_message("info", 'ComingSoonTitle', False, 'ComingSoonBody')


def button_all_VCF():
    wgse_message("info", 'ComingSoonTitle', False, 'ComingSoonBody')


def _button_unindex_BAM():
    if wgse.BAM and wgse.BAM.Indexed and wgse.BAM.file_type in ["BAM", "CRAM"]:
        bam = wgse.BAM.file_oFN + (".bai" if wgse.BAM.file_type == "BAM" else ".crai") # if wgse.BAM.file_type == "CRAM"
        f'"{wgse.outdir.oFPB}_unmapped.bam"'
        wgse.tempf.list.append(bam)
        wgse.BAM.Indexed = False
        set_BAM_window_settings()       # Do not need full new BAM file setup / read-in; just process for Index file gone

    resume_main_window()


def _button_unsort_BAM():
    """ Special internal call enabled only in DEBUG_MODE; simply updates the existing BAM Header to unsorted. """
    if wgse.BAM and wgse.BAM.Sorted:
        samtools = wgse.samtoolsx_qFN
        sed = wgse.sedx_qFN
        newhead = f'"{wgse.tempf.FP}newhead.sam"'
        bam = wgse.BAM.file_qFN
        if "_sorted" in bam:    # Try to avoid adding conflicting names
            unsortbam = bam.replace("_sorted", "_unsorted")
        else:                   # Just append unsorted status
            unsortbam = f'"{wgse.BAM.file_FPB}_unsorted.bam"'

        # CRAM decode requires reference genome be specified with VIEW command
        cram_opt = f'-T {wgse.BAM.Refgenome_qFN}' if wgse.BAM.file_type == "CRAM" else ""
        if wgse.BAM.file_type == "CRAM" and not wgse.reflib.check_reference_genome_qFN(wgse.BAM.Refgenome_qFN):
            resume_main_window()
            return  # Check routine reports error if reference does not exist

        # Will only work where sed is in the PATH; cannot use quoted, path-based sed call in Samtools it seems (Win10)
        # cramopts = "-i" if wgse.BAM.file_type == "CRAM" else ""    # will do in place but we want a copy
        # commands = f"""{samtools} reheader {cramopts} -c 'sed "s/coordinate/unsorted/"' {bam} > {unsortbam}\n"""
        commands = (
            f'{samtools} view -H {cram_opt} {bam} | {sed} "s/SO:coordinate/SO:unsorted/" > {newhead}'
            f'{samtools} reheader {newhead} {bam} > {unsortbam}\n'
        )

        run_bash_script("UnsortBAM", commands)

        button_select_BAM_file(unquote(unsortbam))     # Replace current BAM and garbage collect old

    elif wgse.DEBUG_MODE:
        DEBUG("Hey, if running in DEBUG_MODE, you should know to select a Sorted BAM file first!")
    else:
        pass    # No BAM defined nor DEBUG_MODE on; how did we get here?

    resume_main_window()


def button_language_reload():
    """
    Setup in DEBUG_MODE section for developers to reload the languages.csv file (presumably after changes)
    and recreate the mainwindow like occurs after a language change. This along with language setting button allows
    for quick loops of language translation development (test by quick change and reload)
    :return:
    """
    from utilities import LanguageStrings       # Localized import to prevent module import loop

    wgse.save_settings()        # Make sure we have saved the latest settings; like on exit

    # Restart the language subsystem by doing a new class init and assigning to global setting; garbage collect old
    wgse.lang = LanguageStrings(wgse.language_oFN)  # (re)Start language subsystem (utilities.py); reread language.csv
    wgse.lang.change_language(wgse.lang.language)   # re-setup selected language
    resume_main_window()


def button_change_language():
    """
    User button that displays currently set language. But when clicked, requests the new language setting from the user
    It then proceeds to do a switch_language setting which resets the language setting class and causes a new
    mainWindow to be created and drawn for the new language
    :return:
    """
    # Todo move GUI / get_language from utilities/Language to this module. Then simply call get_language.
    #   Call get_language in this module if settings does not include language at startup like before.
    button_set_language()     # Do everything locally in the language class

    resume_main_window()


def button_set_language():
    """
       Language subsystem user pop-up query for language setting to use. Note, if no stored language in settings,
       then this is called before anything else so we have a language for the initial program window. Hence always
       creating our own root window.
    """

    DEBUG("Displaying Language Selection Pop-up")
    if wgse.window:             # Withdraw main window as will destroy and create new one with new language
        wgse.window.withdraw()
        selectLanguageWindow = Toplevel(wgse.window)
        selectLanguageWindow.transient()
    else:
        selectLanguageWindow = Tk()
    selectLanguageWindow.title(wgse.lang.i18n["RequestLanguage"] if wgse.lang.language else "Please Select Language")
    selectLanguageWindow.geometry("")
    selectLanguageWindow.protocol("WM_DELETE_WINDOW", lambda win=selectLanguageWindow: button_close(win, False))
    selectLanguageWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    selectLanguageWindow.columnconfigure(0, weight=1)
    selectLanguageWindow.rowconfigure(0, weight=1)

    crow = 0
    for lang in wgse.lang.avail_langs:       # Need all languages so go back to original class dictionary
        Label(selectLanguageWindow, text=wgse.lang.request[lang],
              font=font['14']).grid(column=0, row=crow, padx=5, pady=2)
        Button(selectLanguageWindow, text=lang, font=font['14'],
              command=lambda l=lang: wgse.lang.switch_language(l)).grid(column=1, row=crow, padx=5, pady=2)
        crow += 1

    wgse.lang.selectLanguageWindow = selectLanguageWindow   # So we can close window in switch_language()
    selectLanguageWindow.update()
    selectLanguageWindow.grab_set()     # Toplevel equivalent of mainloop
    if wgse.window:           # If wgse.dnaImage set, have run setup_mainWindow
        selectLanguageWindow.wait_window(selectLanguageWindow)
    else:
        selectLanguageWindow.mainloop()  # Run because we ask for language before setting up main window and its loop
    DEBUG("Finally returned from button_set_language")


def button_set_referencelibrary():
    """
    Button to allow override of default Reference_Library location.  Note: user selects but must already exist and
    have content. We only check for the process_reference_genomes.sh script file.
    """
    initialdir = os.path.dirname(wgse.reflib.FP[:-1]) if wgse.reflib and wgse.reflib.FP else \
                 os.path.dirname(wgse.reflib.default_FP[:-1]) if wgse.reflib and wgse.reflib.default_FP else \
                 wgse.install_FP
    reflib_FP = filedialog.askdirectory(parent=wgse.window, title=wgse.lang.i18n['SelectReferenceLibrary'],
                                        initialdir=initialdir)  # Returns Unix/Universal, and no trailing slash
    if not reflib_FP:      # If still not set ...
        DEBUG(f"New Reference Library Directory not set (cancelled)")
    else:
        reflib_FP += '/' if reflib_FP[-1] != '/' else ''    # Note: not NativeOS so forward slash always

        wgse.reflib.change(reflib_FP)
        reflibDirectoryButton.configure(text=wgse.reflib.FB if wgse.reflib.set else wgse.lang.i18n['Default'])

    resume_main_window()


def button_set_tempdir():
    """
    Button to allow override of default Temporary Files location.  Note: user selects and must already exist.
    No files to check for as should be empty to start.  Will clean out old one before replacing.
    """
    initialdir = os.path.dirname(wgse.tempf.FP[:-1]) if wgse.tempf and wgse.tempf.FP else \
                 os.path.dirname(wgse.tempf.default_FP[:-1]) if wgse.tempf and wgse.tempf.default_FP else \
                 wgse.install_FP
    temp_FP = filedialog.askdirectory(parent=wgse.window, title=wgse.lang.i18n['SelectTemporaryFiles'],
                                        initialdir=initialdir)  # Returns Unix/Universal, and no trailing slash
    if not temp_FP:  # If still not set ...
        # User likely hit exit or cancel; simply return (do not clear old value)
        DEBUG(f"New Temporary File Directory not set (cancelled)")
    else:
        # Assure trailing slash
        temp_FP += '/' if temp_FP[-1] != '/' else ''    # Note: universalOS so forward slash always

        wgse.tempf.clean(True)                      # Clean old area before leaving
        wgse.tempf.change(temp_FP)                  # Setup new Temporary Files area (checking validity first)
        tempDirectoryButton.configure(text=wgse.tempf.FB if wgse.tempf.set else wgse.lang.i18n['Default'])

    resume_main_window()


def button_wsl_bwa_patch():
    """ Toggle wsl_bwa_patch setting on or off.  True is On / Active / use WSL BWA."""

    wgse.wsl_bwa_patch = not wgse.wsl_bwa_patch
    wslbwaButton.configure(text=wgse.lang.i18n['Active'] if wgse.wsl_bwa_patch else wgse.lang.i18n['Inactive'])
    resume_main_window()


def init_mainWindow():
    """
    mainWindow (all GUI processing) subsystem.  Not yet a class but preparing for that.  So this is the early __init__.
    withdraw until cal setup_mainWindow to fill it in later
    """
    # Create root / main window so toplevel dialogs can be made
    root = Tk()         # Like for mainloop, should only be one in whole application. This is it.
    root.title("WGS Extract")
    root.geometry(wgse.mainwindow_size)
    root.protocol("WM_DELETE_WINDOW", button_exit)
    root.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    root.columnconfigure(0, weight=1)
    root.rowconfigure(1, weight=1)

    # Let's hide window until setup_mainWindow fills it in. Created now so can have dialogs like language selection.
    root.withdraw()
    return root


def setup_mainWindow():
    """
    Setup of TK main window with all the pane's, buttons and labels. Will be looped over for rest of program life.
    Could be renamed to change_mainWindow() to be similar to other classes.  Only called again when recreating
    after a language change.

    wgse.window not null says init_mainWindow() call made.  wgse.dnaImage not null says this setup_mainWindow called.
    """
    global all_action_buttons

    # Need to globally reset labels and buttons for Settings (Output Dir, BAM File)
    global documentationButton, exitButton, titleFileLabel
    global outputDirectoryLabel, outputDirectoryButton, langSelectLabel, langSelectButton, lreloadButton
    global reflibDirectoryLabel, reflibDirectoryButton, tempDirectoryLabel, tempDirectoryButton
    global bamSelectedLabel, bamSelectButton, bamReferenceGenomeLabel
    # Need to globally adjust summary stats of BAM (and enable/disable its detailed stats button)
    global bamMapAvgReadDepthLabel, bamGenderLabel, bamChromsLabel, bamFsizeLabel      # bamAverageReadLengthLabel,
    global bamIdxstatsButton, bamHeaderButton, bamSortButton, bamIndexButton, bamConvertButton, bamRealignButton
    # Need to globally adjust the action buttons on Extract and Analysis tabs
    global autosomalFormatsButton, autosomalVCFButton, mitoFASTAButton, mitoBAMButton, mitoVCFButton
    global yANDmtButton, yOnlyButton, yVCFButton
    global haplogroupYButton, haplogroupMtButton, exportUnmappedReadsButton
    # Debug_MODE only buttons
    global wslbwaButton, bamUnalignButton, runMicroParallelButton

    # dna frame (Program Title area above Tabs)
    dnaFrame = Frame(wgse.window, relief="solid")
    dnaFrame.grid(row=0, padx=5, pady=2)
    dnaFrame.columnconfigure(0, weight=1)
    dnaFrame.rowconfigure(0, weight=1)

    # Logo
    wgse.dnaImage = ImageTk.PhotoImage(Image.open(wgse.image_oFP))      # Needs stability outside this call
    dnaLabel = Label(dnaFrame, compound="center", text="", image=wgse.dnaImage)
    dnaLabel.grid(row=0, column=0, padx=1, pady=1)  # command=lambda: webbrowser.open_new(wgse.manual_url)

    # Title (program name, version, manual / exit buttons, current set file name)
    headlineFrame = Frame(dnaFrame)
    headlineFrame.grid(row=0, column=1, padx=1)

    Label(headlineFrame, text="WGS Extract", font=font['32']).grid(column=0, row=0, columnspan=2)
    Label(headlineFrame, text=wgse.__version__, font=font['12']).grid(column=0, row=1, columnspan=2)

    documentationButton = Button(headlineFrame, text=wgse.lang.i18n['WGSExtractManual'],
                                 font=font['13'], command=lambda: webbrowser.open_new(wgse.manual_url))
    documentationButton.grid(column=0, row=2)
    exitButton = Button(headlineFrame, text=wgse.lang.i18n['ExitProgram'], font=font['13'], command=button_exit)
    exitButton.grid(column=1, row=2)

    titleFileLabel = Label(headlineFrame, text="", font=font['12'])
    titleFileLabel.grid(column=0, row=3, columnspan=2, pady=3)

    # Setup Tabs area (notebook sub-frames) for Main Window (TTK feature)
    tabParent = Notebook(wgse.window)
    tabParent.grid(row=1, padx=10, pady=5)
    tabParent.columnconfigure(0, weight=1)
    tabParent.rowconfigure(0, weight=1)
    tabParent.enable_traversal()

    # Setup each Tab as a Frame in the Parent Notebook (TTK Feature)
    tabSettings = Frame(tabParent)
    tabExtract = Frame(tabParent)
    tabAnalyze = Frame(tabParent)
    tabParent.add(tabSettings, text=wgse.lang.i18n['Settings'])
    tabParent.add(tabExtract, text=wgse.lang.i18n['ExtractData'])
    tabParent.add(tabAnalyze, text=wgse.lang.i18n['Analyze'])

    ##################################################################################################################
    # Settings Tab

    # Settings frame  # Todo generalize to add temp dir, language (re)set
    fileSelectFrame = LabelFrame(tabSettings, text=wgse.lang.i18n['Settings'], font=font['14'])
    fileSelectFrame.grid(row=1, padx=10, pady=5, sticky=(N, S, E, W))
    fileSelectFrame.columnconfigure(1, weight=1)
    fileSelectFrame.rowconfigure(3, weight=1)
    crow = 0

    # Output Directory setting
    outputDirectoryLabel = Label(fileSelectFrame, text=wgse.lang.i18n['OutputDirectory'],
                                 font=font['14'])
    outputDirectoryLabel.grid(column=0, row=crow, padx=5, pady=2)
    butlabel = wgse.outdir.FB if wgse.outdir.oFP else wgse.lang.i18n['SelectOutputDirectory']
    outputDirectoryButton = Button(fileSelectFrame, text=butlabel, font=font['14'], command=button_select_output_path)
    outputDirectoryButton.grid(column=1, columnspan=2, row=crow, padx=5, pady=2, sticky=W); crow += 1

    # Reference Library Directory setting (default: reference_library in installation directory)
    reflibDirectoryLabel = Label(fileSelectFrame, text=wgse.lang.i18n['ReferenceLibrary'],
                                 font=font['14'])
    reflibDirectoryLabel.grid(column=0, row=crow, padx=5, pady=2)
    butlabel = wgse.reflib.FB if wgse.reflib.set else wgse.lang.i18n['Default']
    reflibDirectoryButton = Button(fileSelectFrame, text=butlabel, font=font['14'], command=button_set_referencelibrary)
    reflibDirectoryButton.grid(column=1, columnspan=2, row=crow, padx=5, pady=2, sticky=W); crow += 1

    # Temporary Files Directory setting (default: temp directory in installation directory)
    tempDirectoryLabel = Label(fileSelectFrame, text=wgse.lang.i18n['TempDirectory'],
                               font=font['14'])
    tempDirectoryLabel.grid(column=0, row=crow, padx=5, pady=2)
    butlabel = wgse.tempf.FB if wgse.tempf.set else wgse.lang.i18n['Default']
    tempDirectoryButton = Button(fileSelectFrame, text=butlabel, font=font['14'], command=button_set_tempdir)
    tempDirectoryButton.grid(column=1, columnspan=2, row=crow, padx=5, pady=2, sticky=W); crow += 1

    # Language setting
    langSelectLabel = Label(fileSelectFrame, text=wgse.lang.i18n["LanguageSetting"], font=font['14'])
    langSelectLabel.grid(column=0, row=crow, padx=5, pady=2)
    langSelectButton = Button(fileSelectFrame, text=wgse.lang.language, font=font['14'],
                              command=button_change_language)
    if not wgse.DEBUG_MODE:
        langSelectButton.grid(column=1, columnspan=2, row=crow, padx=5, pady=2, sticky=W); crow += 1
    else:
        langSelectButton.grid(column=1, row=crow, padx=5, pady=2, sticky=W)
        lreloadButton = Button(fileSelectFrame, text=wgse.lang.i18n["LanguageReload"], font=font['14'],
                               command=button_language_reload)
        lreloadButton.grid(column=2, row=crow, padx=5, pady=2); crow += 1

    # BAM File Selection and Stats
    fileInfoFrame = LabelFrame(tabSettings, text=wgse.lang.i18n['Overview'], font=font['14'])
    fileInfoFrame.grid(row=2, padx=10, pady=5, sticky=(N, S, E, W))
    fileInfoFrame.columnconfigure(1, weight=1)
    fileInfoFrame.rowconfigure(0, weight=1)
    crow = 0

    # We columnspan all the below as we now want four columns for the buttons at the bottom.
    fileInfoFrame.grid_columnconfigure(1, weight=1, uniform="mybuttons")
    fileInfoFrame.grid_columnconfigure(2, weight=1, uniform="mybuttons")
    fileInfoFrame.grid_columnconfigure(3, weight=1, uniform="mybuttons")
    fileInfoColSpan = 3

    # May have restored BAM from settings before calling here; or recreating window due to language change
    bamSelectedLabel = Label(fileInfoFrame, text=wgse.lang.i18n['FilenameOfSourceBam'], font=font['14'])
    bamSelectedLabel.grid(column=0, row=crow, padx=5, pady=2)
    bamSelectButton = Button(fileInfoFrame, text=wgse.lang.i18n['SelectBamFile'], font=font['14'],
                             command=button_select_BAM_file)
    bamSelectButton.grid(column=1, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan, sticky=W); crow += 1

    infoLabelReferenceGenomeOfBam = Label(fileInfoFrame,
                                          text=wgse.lang.i18n['BAMAlignedToReferenceGenome'], font=font['14'])
    infoLabelReferenceGenomeOfBam.grid(column=0, row=crow, padx=5, pady=2)
    bamReferenceGenomeLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamReferenceGenomeLabel.grid(column=1, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    infoLabelMapAvgReadDepth = Label(fileInfoFrame, text=wgse.lang.i18n['MapAvgReadDepthNoN'] + ":", font=font['14'])
    infoLabelMapAvgReadDepth.grid(column=0, row=crow, padx=5, pady=2)
    bamMapAvgReadDepthLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamMapAvgReadDepthLabel.grid(column=1, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    infoLabelGender = Label(fileInfoFrame, text=wgse.lang.i18n['Gender'] + ":", font=font['14'])
    infoLabelGender.grid(column=0, row=crow, padx=5, pady=2)
    bamGenderLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamGenderLabel.grid(column=1, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    infoLabelChroms = Label(fileInfoFrame, text=wgse.lang.i18n['Chroms'] + ":", font=font['14'])
    infoLabelChroms.grid(column=0, row=crow, padx=5, pady=2)
    bamChromsLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamChromsLabel.grid(column=1, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    infoLabelFsize = Label(fileInfoFrame,
                           text=f'{wgse.lang.i18n["FileStats"]}:', font=font['14'])
    infoLabelFsize.grid(column=0, row=crow, padx=5, pady=2)
    bamFsizeLabel = Label(fileInfoFrame, text="", font=font['14b'])
    bamFsizeLabel.grid(column=1, row=crow, padx=5, pady=2, columnspan=fileInfoColSpan); crow += 1

    bamIdxstatsLabel = Label(fileInfoFrame, text=wgse.lang.i18n['BamCheckSamtoolsDescription'], font=font['14'])
    bamIdxstatsLabel.grid(column=0, row=crow, padx=5, pady=2)

    # Stats button is special. Normally we call stats when BAM selected and this just shows results. But if not
    #  sorted or indexed, or a CRAM file, then stats is not run by default. This button overrides that and runs
    #  stats anyway. But we encourage the user through the warning to hit the sort and/or index button or convert the
    #  CRAM to a BAM first instead of forcing the stats to run otherwise here.
    bamIdxstatsButton = Button(fileInfoFrame, text=wgse.lang.i18n['BamCheckSamtoolsButton'], font=font['14'],
                               command=button_stats_BAM, state="disabled")
    bamIdxstatsButton.grid(column=1, row=crow, padx=5, pady=2)

    # Header view button
    bamHeaderButton = Button(fileInfoFrame, text=wgse.lang.i18n['BamHeader'], font=font['14'],
                             command=button_header_BAM, state="disabled")
    bamHeaderButton.grid(column=1, row=crow+1, padx=5, pady=2)

    # Sort and Index buttons are special.  Disabled and set as done if BAM already that way.
    # Will be enabled and set as to do when BAM is not that way.
    bamIndexButton = Button(fileInfoFrame, text=wgse.lang.i18n['Indexed'], font=font['14'],
                            command=button_index_BAM, state="disabled")
    bamIndexButton.grid(column=2, row=crow, padx=5, pady=2)
    bamIndexButton.grid_remove()

    bamSortButton = Button(fileInfoFrame, text=wgse.lang.i18n['Sorted'], font=font['14'],
                           command=button_sort_BAM, state="disabled")
    bamSortButton.grid(column=2, row=crow+1, padx=5, pady=2)
    bamSortButton.grid_remove()

    # Just setup now as if a BAM file specified; will set appropriately before grid() to display
    bamConvertButton = Button(fileInfoFrame, text=wgse.lang.i18n['ToCRAM'], font=font['14'],
                              command=button_BAM_to_CRAM, state="disabled")
    bamConvertButton.grid(column=3, row=crow, padx=5, pady=2)
    bamConvertButton.grid_remove()

    bamRealignButton = Button(fileInfoFrame, text=wgse.lang.i18n['Realign'], font=font['14'],
                              command=button_realign_BAM, state="disabled")
    bamRealignButton.grid(column=3, row=crow+1, padx=5, pady=2)
    bamRealignButton.grid_remove()


    ##################################################################################################################
    # Extract Tab

    autosomesLabel = Label(tabExtract, text=wgse.lang.i18n['SuggestedFor'], font=font['14u'])
    autosomesLabel.grid(column=0, row=0, padx=5, pady=2)

    autosomesLabel = Label(tabExtract, text=wgse.lang.i18n['GenerateFile'], font=font['14u'])
    autosomesLabel.grid(column=1, row=0, padx=5, pady=2)

    #------------ Microarray / Autosomal Frame -------------------------
    autosomesFrame = LabelFrame(tabExtract, text=wgse.lang.i18n['Autosomes'], font=font['16'])
    autosomesFrame.grid(row=2, columnspan=2, padx=10, pady=2, sticky=(N, S, E, W))
    autosomesFrame.columnconfigure(0, weight=1)
    autosomesFrame.rowconfigure(1, weight=1)

    autosomesLabel = Label(autosomesFrame, text=wgse.lang.i18n['AutosomesDescr'], font=font['14'])
    autosomesLabel.grid(column=0, row=0, rowspan=2, padx=5, pady=2)

    autosomalFormatsButton = Button(autosomesFrame, text=wgse.lang.i18n['GoToSelectAutosomalFormats'], font=font['14'],
                                    command=button_select_autosomal_formats, state="disabled")
    autosomalFormatsButton.grid(column=1, row=0, padx=5, pady=2)

    autosomalVCFButton = Button(autosomesFrame, text=wgse.lang.i18n['GenerateVCFSNP'], font=font['14'],
                                    command=button_all_VCF, state="disabled")
    autosomalVCFButton.grid(column=1, row=1, padx=5, pady=2)

    autosomesTrailerLabel = Label(autosomesFrame, text=wgse.lang.i18n['AutosomesTrailer'], font=font['14'])
    autosomesTrailerLabel.grid(column=0, row=3, columnspan=2, padx=5, pady=2)

    #------------ Mitochondial Frame -------------------------
    mitoFrame = LabelFrame(tabExtract, text=wgse.lang.i18n['MitochondrialDNA'], font=font['16'])
    mitoFrame.grid(row=3, columnspan=2, padx=10, pady=2, sticky=(N, S, E, W))
    mitoFrame.columnconfigure(0, weight=1)
    mitoFrame.rowconfigure(0, weight=1)

    mtLabel = Label(mitoFrame, text=wgse.lang.i18n['MitochondrialDescriptionNeededFor'], font=font['14'])
    mtLabel.grid(column=0, row=0, padx=5, pady=2)
    mitoFASTAButton = Button(mitoFrame, text=wgse.lang.i18n['GenerateFastaMtdna'], font=font['14'],
                             command=button_mtdna_FASTA, state="disabled")
    mitoFASTAButton.grid(column=1, row=0, padx=5, pady=2)

    mtLabel2 = Label(mitoFrame, text=wgse.lang.i18n['RequiredForMitoverse'], font=font['14'])
    mtLabel2.grid(column=0, row=1, rowspan=2, padx=5, pady=2)

    mitoBAMButton = Button(mitoFrame, text=wgse.lang.i18n['GenerateBAMOnlyWithMT'], font=font['14'],
                             command=button_mtdna_BAM, state="disabled")
    mitoBAMButton.grid(column=1, row=1, padx=5, pady=2)

    mitoVCFButton = Button(mitoFrame, text=wgse.lang.i18n['GenerateVCFOnlyWithMT'], font=font['14'],
                             command=button_mtdna_VCF, state="disabled")
    mitoVCFButton.grid(column=1, row=2, padx=5, pady=2)

    #mtLabel2 = Label(mitoFrame, text=wgse.lang.i18n['MitochondrialDescriptionNeededFor2'], font=font['14'])
    #mtLabel2.grid(column=0, row=3, padx=5, pady=2)

    #------------ Y Chromosome Frame -------------------------
    yFrame = LabelFrame(tabExtract, text=wgse.lang.i18n['YDNA'], font=font['16'])
    yFrame.grid(row=4, columnspan=2, padx=10, pady=2, sticky=(N, S, E, W))
    yFrame.columnconfigure(0, weight=1)
    yFrame.rowconfigure(1, weight=1)

    Label(yFrame, text=wgse.lang.i18n['RequiredForUploadingToYfull'], font=font['14']).grid(column=0, row=0, padx=5, pady=2)
    yANDmtButton = Button(yFrame, text=wgse.lang.i18n['GenerateBAMwithYandMtdna'], font=font['14'], command=button_yAndMt_BAM,
                          state="disabled")
    yANDmtButton.grid(column=1, row=0, padx=5, pady=2)

    Label(yFrame, text=wgse.lang.i18n['RequiredForYDNAWarehouse'], font=font['14']).grid(column=0, row=1, padx=5, pady=2)
    yOnlyButton = Button(yFrame, text=wgse.lang.i18n['GenerateBAMOnlyWithY'], font=font['14'], command=button_yOnly_BAM,
                         state="disabled")
    yOnlyButton.grid(column=1, row=1, padx=5, pady=2)

    Label(yFrame, text=wgse.lang.i18n['RequiredForCladefinder'], font=font['14']).grid(column=0, row=2, padx=5, pady=2)
    yVCFButton = Button(yFrame, text=wgse.lang.i18n['GenerateVCFOnlyWithY'], font=font['14'], command=button_yOnly_VCF,
                        state="disabled")
    yVCFButton.grid(column=1, row=2, padx=5, pady=2)

    Label(yFrame, text=wgse.lang.i18n['RequiredForUploadingToYfull2'],
          font=font['14']).grid(columnspan=2, row=3, padx=5, pady=2)

    ##################################################################################################################
    # Analyze Tab
    # Moved BAM Detailed stats to main settings tab where BAM is selected

    haplogroupFrame = LabelFrame(tabAnalyze, text=wgse.lang.i18n['Haplogroups'], font=font['16'])
    haplogroupFrame.grid(row=2, padx=10, pady=5, sticky=(N, S, E, W))
    haplogroupFrame.columnconfigure(0, weight=1)
    haplogroupFrame.rowconfigure(1, weight=1)

    Label(haplogroupFrame, text=wgse.lang.i18n['DetermineHaplogroups'],
          font=font['14']).grid(column=0, row=0, padx=5, pady=2)
    haplogroupYButton = Button(haplogroupFrame, text=wgse.lang.i18n['YDNA'], font=font['14'],
                               command=button_ydna_haplogroup, state="disabled")
    haplogroupYButton.grid(column=1, row=0, padx=5, pady=2)
    haplogroupMtButton = Button(haplogroupFrame, text=wgse.lang.i18n['MitochondrialDNA'], font=font['14'],
                                command=button_mtdna_haplogroup, state="disabled")
    haplogroupMtButton.grid(column=2, row=0, padx=5, pady=2)
    Label(haplogroupFrame, text=wgse.lang.i18n['DetermineHaplogroups2'],
          font=font['14']).grid(columnspan=3, row=1, padx=5, pady=2)

    oralMicrobiomeFrame = LabelFrame(tabAnalyze, text=wgse.lang.i18n['OralMicrobiome'], font=font['14'])
    oralMicrobiomeFrame.grid(row=3, padx=10, pady=5, sticky=(N, S, E, W))
    oralMicrobiomeFrame.columnconfigure(0, weight=1)
    oralMicrobiomeFrame.rowconfigure(1, weight=1)

    Label(oralMicrobiomeFrame, text=wgse.lang.i18n['ForUploadsCosmosId'],
          font=font['14']).grid(column=0, row=0, padx=5, pady=2)
    exportUnmappedReadsButton = Button(oralMicrobiomeFrame, text=wgse.lang.i18n['ExportUnmappedReads'], font=font['14'],
                                       command=button_export_unmapped_reads, state="disabled")
    exportUnmappedReadsButton.grid(column=1, row=0, padx=5, pady=2)

    # Debug Frame will only appear (get gridded) if DEBUG_MODE is on; but still need to create first here
    debugFrame = LabelFrame(tabAnalyze, text="DEBUG_MODE", font=font['14'])
    debugFrame.columnconfigure(0, weight=1)
    debugFrame.rowconfigure(1, weight=1)

    # Ugly; for translators only. Static language selection buttons in case they really screw up the current language
    lselLabel = Label(debugFrame, text=wgse.lang.i18n["RequestLanguage"], font=font['14'])
    lselButton = [x for x in range(len(wgse.lang.avail_langs))]
    column = 0
    for lang in wgse.lang.avail_langs:  # Need all languages so go back to original _lang dictionary
        # lselLabel[column] = Label(debugFrame, text=wgse.lang._lang[lang]["RequestLanguage"], font=font['14'])
        # noinspection PyTypeChecker
        lselButton[column] = Button(debugFrame, text=wgse.lang.avail_langs[column], font=font['14'],
                                    command=lambda l=lang: wgse.lang.switch_language(l))
        column += 1

    # Windows 10 WSL2 BWA Aligner Override
    wslbwaLabel = Label(debugFrame, text=wgse.lang.i18n['Win10WslBwaOverride'],
                                 font=font['14'])
    wsllabel = wgse.lang.i18n['Active'] if wgse.wsl_bwa_patch else wgse.lang.i18n['Inactive']
    wslbwaButton = Button(debugFrame, text=wsllabel, font=font['14'], command=button_wsl_bwa_patch)

    bamUnalignLabel = Label(debugFrame, text=wgse.lang.i18n['UnalignTitle'], font=font['14'])
    bamUnalignButton = Button(debugFrame, text=wgse.lang.i18n['Unalign'], font=font['14'],
                              command=button_unalign_BAM, state="disabled")     # Action button; enabled with others

    runMicroParallelLabel = Label(debugFrame, text="Generate CombinedKit (Parallel)", font=font['14'])
    runMicroParallelButton = Button(debugFrame, text=wgse.lang.i18n['buttonGenerateSelectedFiles'],
                font=font['14'], command=lambda win=wgse.window: _button_CombinedKit(win, parallel=True))

    # debugFrame will only become visible if we grid so ...
    if wgse.DEBUG_MODE:
        debugFrame.grid(row=4, padx=10, pady=5, sticky=(N, S, E, W))
        lselLabel.grid(row=0, column=0, padx=5, pady=2)

        # Create rows of 4 buttons each for various languages; large lanugage names span columns
        crow = 0; ccol = 0
        for numlangs in range(len(wgse.lang.avail_langs)):
            # Figure out if button needs more than one column (each column is 10 characters max); max columns is 4
            # noinspection PyUnresolvedReferences
            span = min((len(lselButton[numlangs].cget("text")) // 10) + 1, 4)
            if ccol + span > 4:     # Not enough room in row left for this button; start new row
                crow += 1
                ccol = 0
            # noinspection PyUnresolvedReferences
            lselButton[numlangs].grid(row=crow, column=ccol+1, columnspan=span, padx=5, pady=2)
            ccol = (ccol + span) % 4
            crow += 1 if ccol == 0 else 0

        crow += 1
        wslbwaLabel.grid(row=crow, column=0, columnspan=3, padx=5, pady=2)
        wslbwaButton.grid(row=crow, column=3, padx=5, pady=2, sticky=W)

        crow += 1
        bamUnalignLabel.grid(row=crow, column=0, columnspan=3, padx=5, pady=2)
        bamUnalignButton.grid(row=crow, column=3, padx=5, pady=2)

        crow += 1
        runMicroParallelLabel.grid(row=crow, column=0, columnspan=3, padx=5, pady=2)
        runMicroParallelButton.grid(row=crow, column=3, padx=5, pady=2)

    # Todo add theme setting and button to change it  (in development; may change GUI from tkinter before completing)
    ''' Developing code to allow a theme change setting
    style = ttk.Style()
    cur_theme = style.theme_use()
    themes = StringVar(value=style.named_themes())

    Label(debugFrame, text=f"Change Theme (current highlighted):", font=font['14']).grid(column=0, row=0, padx=5, pady=2)
    select_theme = Listbox(debugFrame, listvariable=themes, selectmode="browse", height=1)
    select_theme.grid(column=1, row=0, padx=5, pady=2)
    select_theme.selection_set(themes.index[cur_theme])
    select_theme.bind('<<ListboxSelect>>', change_theme)
    '''

    # Listing of main window action buttons that may need the state changed based on content of the loaded BAM file
    # Does not include Settings buttons nor BAM file button itself which are handled separately.
    all_action_buttons = [
        # Settings Tab (inside BAM File setting)
        bamIdxstatsButton,
        bamHeaderButton,
        bamIndexButton,
        bamSortButton,
        bamConvertButton,
        bamRealignButton,
        # Extract Data Tab
        autosomalFormatsButton,
        autosomalVCFButton,
        mitoFASTAButton,
        mitoBAMButton,
        mitoVCFButton,
        yANDmtButton,
        yOnlyButton,
        yVCFButton,
        # Analyze Tab
        haplogroupYButton,
        haplogroupMtButton,
        exportUnmappedReadsButton,
        bamUnalignButton,
        runMicroParallelButton
    ]

    # We may have called setup_mainWindow() after already loading saved settings (including BAMFile) or
    #   due to a language change.  In these cases, need to modify default setup for already set values
    set_BAM_window_settings()       # If BAM not set from stored settings; simply makes sure settings are cleared
    resume_main_window()            # Includes call to update_action_buttons()
    return wgse.window
