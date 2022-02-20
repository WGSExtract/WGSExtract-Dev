# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import sys      # for argv[] if called as separate process / program
import os       # for os.path
from tkinter import Toplevel, BooleanVar, W, Checkbutton, Button, Label, LabelFrame, Frame

import settings as wgse
font = wgse.font
from utilities import nativeOS, DEBUG, wgse_message
from commandprocessor import run_bash_script
# from mainwindow import resume_main_window      # Localized inside cancel_autosomal_formats_window() due to loop
# Todo Microarray Format Window update -- as tab to main window? New standalone program and class?

"""###################################################################################################################
  Microarray File Export subsystem / window (aka Autosomal File Formats)
    This microarray files generator developed with WGSE Extract by City Farmer.  The largest, most
    complete and unique contribution of the WGS Extract tool.  Grew out of the basic concept in Extract23 that took a
    30x WGS and generated a 23andMe v3 file.  Now based on a generated CombinedKit template file (in 37/38 reference
    format) which is then subsetted for each of 12+ formats (using a new module "aconv" developed here). Also does 
    a liftover (based on wrapping pyliftover rewrite of UCSC's liftover) from Build 38 to 37.
"""
checkbuttons_results = []

combinedButton      = None
g23andmev3Button    = None
g23andmev4Button    = None
g23andmev5Button    = None
g23andmeAPIButton   = None
ancestryV1Button    = None
ancestryV2Button    = None
ftdnav2Button       = None
ftdnav3Button       = None
ldnav1Button        = None
ldnav2Button        = None
myheritagev1Button  = None
myheritagev2Button  = None

combinedButtonResult      = None
g23andmev3ButtonResult    = None
g23andmev4ButtonResult    = None
g23andmev5ButtonResult    = None
g23andmeAPIButtonResult   = None
ancestryV1ButtonResult    = None
ancestryV2ButtonResult    = None
ftdnav2ButtonResult       = None
ftdnav3ButtonResult       = None
ldnav1ButtonResult        = None
ldnav2ButtonResult        = None
myheritagev1ButtonResult  = None
myheritagev2ButtonResult  = None


def button_select_autosomal_formats():
    """
        Main entry into the microarray file generator from WGSE (or as a standalone with locally recreated settings).
        Generates main selection window and then waits on User direction. Utilizes selected (in settings) BAM and
        Output directory area.  Final output are microarray files selected by user.
        Note: if finds an existing CombinedKit file, will start with that instead of recreating it
    """
    global combinedButton
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton
    global ancestryV1Button, ancestryV2Button, ftdnav2Button, ftdnav3Button
    global ldnav1Button, ldnav2Button, myheritagev1Button, myheritagev2Button
    global combinedButtonResult
    global g23andmev3ButtonResult, g23andmev4ButtonResult, g23andmev5ButtonResult, g23andmeAPIButtonResult
    global ancestryV1ButtonResult, ancestryV2ButtonResult, ftdnav2ButtonResult, ftdnav3ButtonResult
    global ldnav1ButtonResult, ldnav2ButtonResult, myheritagev1ButtonResult, myheritagev2ButtonResult
    global generateSelectedFilesButton, checkbuttons_results, selectAutosomalFormatsWindow

    # Create main selection window for the Autosomal microarray file generator
    #wgse.window.withdraw()
    selectAutosomalFormatsWindow = Toplevel(wgse.window)
    selectAutosomalFormatsWindow.transient()
    selectAutosomalFormatsWindow.protocol("WM_DELETE_WINDOW", cancel_autosomal_formats_window)
    selectAutosomalFormatsWindow.title(wgse.lang.i18n['TitleSelectAutosomalFormats'])
    selectAutosomalFormatsWindow.geometry(wgse.microarr_winsize)
    selectAutosomalFormatsWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    selectAutosomalFormatsWindow.columnconfigure(0, weight=1)
    selectAutosomalFormatsWindow.rowconfigure(0, weight=1)

    questionAutosomalFormats = Label(selectAutosomalFormatsWindow, text=wgse.lang.i18n['QuestionWhichAutosomalFormats'],
                                     font=font['16b'], anchor="nw", justify="left")
    questionAutosomalFormats.grid(column=0, row=0, padx=5, pady=5)

    # Create Frames for each company format and buttons across the bottom
    combinedFrame = LabelFrame(selectAutosomalFormatsWindow, text=wgse.lang.i18n['everythingCombined'], font=font['14'])
    combinedFrame.grid(row=1, padx=5, pady=3, sticky='W')
    combinedFrame.columnconfigure(0, weight=1)
    combinedFrame.rowconfigure(1, weight=1)
    g23andmeFrame = LabelFrame(selectAutosomalFormatsWindow, text="23andMe", font=font['14'])
    g23andmeFrame.grid(row=2, padx=5, pady=3, sticky='W')
    g23andmeFrame.columnconfigure(0, weight=1)
    g23andmeFrame.rowconfigure(1, weight=1)
    ancestryFrame = LabelFrame(selectAutosomalFormatsWindow, text="Ancestry", font=font['14'])
    ancestryFrame.grid(row=3, padx=5, pady=3, sticky='W')
    ancestryFrame.columnconfigure(0, weight=1)
    ancestryFrame.rowconfigure(1, weight=1)
    ftdnaFrame = LabelFrame(selectAutosomalFormatsWindow, text="Family Tree DNA", font=font['14'])
    ftdnaFrame.grid(row=4, padx=5, pady=3, sticky='W')
    ftdnaFrame.columnconfigure(0, weight=1)
    ftdnaFrame.rowconfigure(1, weight=1)
    ldnaFrame = LabelFrame(selectAutosomalFormatsWindow, text="Living DNA", font=font['14'])
    ldnaFrame.grid(row=5, padx=5, pady=3, sticky='W')
    ldnaFrame.columnconfigure(0, weight=1)
    ldnaFrame.rowconfigure(1, weight=1)
    myheritageFrame = LabelFrame(selectAutosomalFormatsWindow, text="MyHeritage", font=font['14'])
    myheritageFrame.grid(row=6, padx=5, pady=3, sticky='W')
    myheritageFrame.columnconfigure(0, weight=1)
    myheritageFrame.rowconfigure(1, weight=1)
    buttonsFrame = Frame(selectAutosomalFormatsWindow)
    buttonsFrame.grid(row=7, padx=5, pady=3, sticky='W')
    buttonsFrame.columnconfigure(0, weight=1)
    buttonsFrame.rowconfigure(1, weight=1)

    # Generate select box buttons for each available format
    checkbuttons_results = []

    combinedSpacer = Label(combinedFrame, text="     ", font=font['14'], anchor="w")
    combinedSpacer.grid(row=0, column=0, padx=5, pady=2)

    combinedButtonResult = BooleanVar()
    combinedButtonResult.set(False)
    checkbuttons_results.append(False)
    combinedButton = Checkbutton(combinedFrame, variable=combinedButtonResult, font=font['14'], anchor=W,
                text=wgse.lang.i18n['CombinedFileAllSNPs'],
                command=lambda:adna_checkbutton_click(0))
    combinedButton.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    combinedTrailingSpacer = Label(combinedFrame, text="     ", font=font['14'], anchor="w")
    combinedTrailingSpacer.grid(row=0, column=2, padx=5, pady=2)

    g23andmeSpacer = Label(g23andmeFrame, text="     ", font=font['14'], anchor="w")
    g23andmeSpacer.grid(row=0, column=0, padx=5, pady=2)

    g23andmev3ButtonResult = BooleanVar()
    g23andmev3ButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmev3Button = Checkbutton(g23andmeFrame, variable=g23andmev3ButtonResult, font=font['14'], anchor=W,
                text='23andMe v3 (11/2010 - 11/2013) (' + wgse.lang.i18n['RECOMMENDED'] + ')',
                command=lambda:adna_checkbutton_click(1))
    g23andmev3Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    g23andmev4ButtonResult = BooleanVar()
    g23andmev4ButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmev4Button = Checkbutton(g23andmeFrame, variable=g23andmev4ButtonResult, font=font['14'], anchor=W,
                text='23andMe v4 (11/2013 - 08/2017)',
                command=lambda:adna_checkbutton_click(2))
    g23andmev4Button.grid(row=1, column=1, padx=1, pady=1, sticky='W')

    g23andmev5ButtonResult = BooleanVar()
    g23andmev5ButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmev5Button = Checkbutton(g23andmeFrame, variable=g23andmev5ButtonResult, font=font['14'], anchor=W,
                text='23andMe v5 (' + wgse.lang.i18n['Since'] + ' 08/2017) (' + wgse.lang.i18n['RECOMMENDED'] + ')',
                command=lambda:adna_checkbutton_click(3))
    g23andmev5Button.grid(row=2, column=1, padx=1, pady=1, sticky='W')

    g23andmeAPIButtonResult = BooleanVar()
    g23andmeAPIButtonResult.set(False)
    checkbuttons_results.append(False)
    g23andmeAPIButton = Checkbutton(g23andmeFrame, variable=g23andmeAPIButtonResult, font=font['14'], anchor=W,
                text=wgse.lang.i18n['23andMeFutureSNPs'],
                command=lambda:adna_checkbutton_click(4))
    g23andmeAPIButton.grid(row=3, column=1, padx=1, pady=1, sticky='W')

    ancestrySpacer = Label(ancestryFrame, text="     ", font=font['14'], anchor="w")
    ancestrySpacer.grid(row=0, column=0, padx=5, pady=2)

    ancestryV1ButtonResult = BooleanVar()
    ancestryV1ButtonResult.set(False)
    checkbuttons_results.append(False)
    ancestryV1Button = Checkbutton(ancestryFrame, variable=ancestryV1ButtonResult, font=font['14'], anchor=W,
                text='Ancestry v1 (01/2012 - 05/2016)',
                command=lambda:adna_checkbutton_click(5))
    ancestryV1Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    ancestryV2ButtonResult = BooleanVar()
    ancestryV2ButtonResult.set(False)
    checkbuttons_results.append(False)
    ancestryV2Button = Checkbutton(ancestryFrame, variable=ancestryV2ButtonResult, font=font['14'], anchor=W,
                text='Ancestry v2 (' + wgse.lang.i18n['Since'] + ' 05/2016)',
                command=lambda:adna_checkbutton_click(6))
    ancestryV2Button.grid(row=1, column=1, padx=1, pady=1, sticky='W')

    ftdnaSpacer = Label(ftdnaFrame, text="     ", font=font['14'], anchor="w")
    ftdnaSpacer.grid(row=0, column=0, padx=5, pady=2)

    ftdnav2ButtonResult = BooleanVar()
    ftdnav2ButtonResult.set(False)
    checkbuttons_results.append(False)
    ftdnav2Button = Checkbutton(ftdnaFrame, variable=ftdnav2ButtonResult, font=font['14'], anchor=W,
                text='FTDNA v2 (02/2011 - 04/2019)',
                command=lambda:adna_checkbutton_click(7))
    ftdnav2Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    ftdnav3ButtonResult = BooleanVar()
    ftdnav3ButtonResult.set(False)
    checkbuttons_results.append(False)
    ftdnav3Button = Checkbutton(ftdnaFrame, variable=ftdnav3ButtonResult, font=font['14'], anchor=W,
                text='FTDNA v3 (' + wgse.lang.i18n['Since'] + ' 04/2019)',
                command=lambda:adna_checkbutton_click(8))
    ftdnav3Button.grid(row=1, column=1, padx=1, pady=1, sticky='W')

    ldnaSpacer = Label(ldnaFrame, text="     ", font=font['14'], anchor="w")
    ldnaSpacer.grid(row=0, column=0, padx=5, pady=2)

    ldnav1ButtonResult = BooleanVar()
    ldnav1ButtonResult.set(False)
    checkbuttons_results.append(False)
    ldnav1Button = Checkbutton(ldnaFrame, variable=ldnav1ButtonResult, font=font['14'], anchor=W,
                text='Living DNA v1 (09/2016 - 10/2018)',
                command=lambda:adna_checkbutton_click(9))
    ldnav1Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    ldnav2ButtonResult = BooleanVar()
    ldnav2ButtonResult.set(False)
    checkbuttons_results.append(False)
    ldnav2Button = Checkbutton(ldnaFrame, variable=ldnav2ButtonResult, font=font['14'], anchor=W,
                text='Living DNA v2 (' + wgse.lang.i18n['Since'] + ' 10/2018)',
                command=lambda:adna_checkbutton_click(10))
    ldnav2Button.grid(row=1, column=1, padx=1, pady=1, sticky='W')

    myheritageSpacer = Label(myheritageFrame, text="     ", font=font['14'], anchor="w")
    myheritageSpacer.grid(row=0, column=0, padx=5, pady=2)

    myheritagev1ButtonResult = BooleanVar()
    myheritagev1ButtonResult.set(False)
    checkbuttons_results.append(False)
    myheritagev1Button = Checkbutton(myheritageFrame, variable=myheritagev1ButtonResult, font=font['14'], anchor=W,
                text='MyHeritage v1 (11/2016 - 03/2019)',
                command=lambda:adna_checkbutton_click(11))
    myheritagev1Button.grid(row=0, column=1, padx=1, pady=1, sticky='W')

    myheritagev2ButtonResult = BooleanVar()
    myheritagev2ButtonResult.set(False)
    checkbuttons_results.append(False)
    myheritagev2Button = Checkbutton(myheritageFrame, variable=myheritagev2ButtonResult, font=font['14'], anchor=W,
                text='MyHeritage v2 (' + wgse.lang.i18n['Since'] + ' 03/2019)',
                command=lambda:adna_checkbutton_click(12))
    myheritagev2Button.grid(row=1, column=1, padx=1, pady=1, sticky='W')

    # Generate action buttons across the bottom row
    selectAllFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonSelectEverything'],
                font=font['14'], command=button_select_every_autosomal_test)
    selectAllFileFormatsButton.grid(column=0, row=0, padx=5, pady=3)

    selectRecFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonSelectRecommended'],
                font=font['14'], command=button_select_recommended_autosomal_test)
    selectRecFileFormatsButton.grid(column=1, row=0, padx=5, pady=3)

    deSelectAllFileFormatsButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonDeselectEverything'],
                font=font['14'], command=button_deselect_every_autosomal_test)
    deSelectAllFileFormatsButton.grid(column=2, row=0, padx=5, pady=3)

    generateSelectedFilesButton = Button(buttonsFrame, text=wgse.lang.i18n['buttonGenerateSelectedFiles'],
                font=font['14'], command=button_generate_selected_autosomal, state="disabled")
    generateSelectedFilesButton.grid(column=4, row=0, padx=5, pady=3)

    # Had previously assumed user knew they could click window exit in upper bar; now provide explicit exit button
    exitButton = Button(buttonsFrame, text=wgse.lang.i18n['Exit'], font=font['14'],
                        command=cancel_autosomal_formats_window)
    exitButton.grid(column=5, row=0, padx=5, pady=3)

    # Select Recommended by default and wait for an action button (generate or exit)
    button_select_recommended_autosomal_test()
    selectAutosomalFormatsWindow.update()
    selectAutosomalFormatsWindow.deiconify()    # Update not showing recommended at start after grab_set() added
    selectAutosomalFormatsWindow.grab_set()
    selectAutosomalFormatsWindow.wait_window()


# noinspection PyUnresolvedReferences
def button_select_recommended_autosomal_test():
    """ Routine to setup Recommended selections (default on entry; can be restored by user button).  """
    global checkbuttons_results, generateSelectedFilesButton
    global combinedButton, g23andmev3Button, g23andmev5Button

    # Unfortunately, no key index. So position dependent. Happens to be 0, 1 and 3.
    button_deselect_every_autosomal_test()
    combinedButton.select();   checkbuttons_results[0] = True
    g23andmev3Button.select(); checkbuttons_results[1] = True
    g23andmev5Button.select(); checkbuttons_results[3] = True
    generateSelectedFilesButton.configure(state="normal")


# noinspection PyUnresolvedReferences
def button_select_every_autosomal_test():
    global checkbuttons_results, generateSelectedFilesButton
    global combinedButton
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton
    global ancestryV1Button, ancestryV2Button, ftdnav2Button, ftdnav3Button
    global ldnav1Button, ldnav2Button, myheritagev1Button, myheritagev2Button

    combinedButton.select()
    g23andmev3Button.select()
    g23andmev4Button.select()
    g23andmev5Button.select()
    g23andmeAPIButton.select()
    ancestryV1Button.select()
    ancestryV2Button.select()
    ftdnav2Button.select()
    ftdnav3Button.select()
    ldnav1Button.select()
    ldnav2Button.select()
    myheritagev1Button.select()
    myheritagev2Button.select()
    for i in range(len(checkbuttons_results)):
        checkbuttons_results[i] = True
    generateSelectedFilesButton.configure(state="normal")


# noinspection PyUnresolvedReferences
def button_deselect_every_autosomal_test():
    global checkbuttons_results, generateSelectedFilesButton
    global combinedButton
    global g23andmev3Button, g23andmev4Button, g23andmev5Button, g23andmeAPIButton
    global ancestryV1Button, ancestryV2Button, ftdnav2Button, ftdnav3Button
    global ldnav1Button, ldnav2Button, myheritagev1Button, myheritagev2Button

    combinedButton.deselect()
    g23andmev3Button.deselect()
    g23andmev4Button.deselect()
    g23andmev5Button.deselect()
    g23andmeAPIButton.deselect()
    ancestryV1Button.deselect()
    ancestryV2Button.deselect()
    ftdnav2Button.deselect()
    ftdnav3Button.deselect()
    ldnav1Button.deselect()
    ldnav2Button.deselect()
    myheritagev1Button.deselect()
    myheritagev2Button.deselect()
    
    for i in range(len(checkbuttons_results)):
        checkbuttons_results[i] = False
    generateSelectedFilesButton.configure(state="disabled")


def adna_checkbutton_click(i):
    global checkbuttons_results, generateSelectedFilesButton

    checkbuttons_results[i] = not checkbuttons_results[i]

    number_of_selected_buttons = 0
    for i in range(len(checkbuttons_results)):
        if (checkbuttons_results[i] == True):
           number_of_selected_buttons += 1
    generateSelectedFilesButton.configure(state="normal" if number_of_selected_buttons > 0 else "disabled")


def get_target_type_suffix(target_type_name_all):
    if "FTDNA" in target_type_name_all or "MyHeritage" in target_type_name_all:
        return ".csv"
    elif "LDNA" in target_type_name_all or "23andMe" in target_type_name_all or "Ancestry" in target_type_name_all:
        return ".txt"
    else:
        return ".err"                                 # REH 20Mar2020 Refactor


def button_generate_selected_autosomal():
    """ Run Microarray format generator (originally based on Extract23) to generate selected formats. """
    global selectAutosomalFormatsWindow

    if wgse.BAM.chrom_types['A'] == 0:
        # Need at least Autosomes present to create microarray file. Button should not have been available but ...
        # Todo should provide a pop-up warning message that doing nothing here
        cancel_autosomal_formats_window()
        return

    if wgse.BAM.Refgenome != "hs37d5":        # Why the warning? Need to document in the manual. HG19 not good?
        wgse_message("warning", 'NoHs37d5WarningTitle', True,
                     wgse.lang.i18n['NoHs37d5WarningText'].replace("{{refgenome}}", wgse.BAM.Refgenome))

    # selectAutosomalFormatsWindow.withdraw()

    # Call internal button to reuse or create a new CombinedKit file in the Output Diretory from the current BAM file
    reuse = _button_CombinedKit(selectAutosomalFormatsWindow)

    # By this point, we should have a .zip and .txt CombinedKit file

    # Make list of requested output formats from checked buttons
    possible_output_formats = ['CombinedKit', '23andMe_V3', '23andMe_V4', '23andMe_V5', '23andMe_SNPs_API',
                               'Ancestry_V1', 'Ancestry_V2', 'FTDNA_V2', 'FTDNA_V3', 'LDNA_V1', 'LDNA_V2',
                               'MyHeritage_V1', 'MyHeritage_V2']  # NOTE: IMPORTANT, order dependent
    output_formats_to_use = []
    for i in range(len(checkbuttons_results)):
        if checkbuttons_results[i] == True:
            output_formats_to_use.append(possible_output_formats[i])

    out_FPB = wgse.outdir.FPB
    out_FPB_cmbkit = f'{wgse.outdir.FPB}_CombinedKit'
    CombinedKitTXT_oFN = nativeOS(f'{out_FPB_cmbkit}.txt')
    CombinedKitZIP_oFN = nativeOS(f'{out_FPB_cmbkit}.zip')

    if os.path.exists(CombinedKitTXT_oFN) and \
        os.path.getmtime(CombinedKitTXT_oFN) > os.path.getmtime(wgse.BAM.file_oFN) and \
        os.path.getsize(CombinedKitTXT_oFN) > 5000000 :
        # Now for EACH format selected (other than CombinedKit), add a run of aconv.py (our script)
        #     to take the CombinedKit content and subset it to create the respective format
        # Note that we are simply running aconv as a stand-alone python program; not importing and calling it.

        python = f'"{wgse.python3x_FN}"'
        aconv = f'"{wgse.prog_FP}aconv.py"'
        cmdzip = wgse.zipx_qFN

        commands = ""
        # Todo modify below to simply call aconv routine, within same python thread, and use python zip command
        for single_output_format in output_formats_to_use:
            if single_output_format != "CombinedKit":
                suffix = get_target_type_suffix(single_output_format)
                out_uncompressed_FN = f'"{out_FPB}_{single_output_format}{suffix}"'
                out_compressed_FN = f'"{out_FPB}_{single_output_format}.zip"'
                commands += (
                    f'echo "Generating microarray file for format {single_output_format}"\n'
                    f'{python} {aconv} {single_output_format} "{out_FPB_cmbkit}.txt" "{out_FPB}"\n'
                    f'{cmdzip} -mj {out_compressed_FN} {out_uncompressed_FN}\n'  # Do a move; delete uncompressed
                )

        run_bash_script("ButtonMicroarrayDNA", commands, parent=selectAutosomalFormatsWindow)

    # Handle CombinedKit file cleanup based on whether requested and/or reused an existing copy
    wgse.tempf.list.append(CombinedKitTXT_oFN)           # Always delete the uncompressed version;
    if not("CombinedKit" in output_formats_to_use or reuse) or os.path.getsize(CombinedKitZIP_oFN) < 5000000:
        wgse.tempf.list.append(CombinedKitZIP_oFN)       # Not reused, not asked to save, nor error as too small

    cancel_autosomal_formats_window()


def _button_CombinedKit(parent_window, parallel=False):
        """
          CombinedKit merged Microarray File Creator (internal button)
          Based on the (original) Extract23 Windows script: https://github.com/tkrahn/extract23/blob/master/extract23.sh
          With parallelization extension shown in: https://gist.github.com/tkrahn/ef62cfaab678f447ea53ddee09ce0eb2
          Original Extract23 did this only for a single 23andMe v3 template.  Here we created a CombinedKit template
          of the merger of all known formats.  Then wrote our own aconv.py module to subset the CombinedKit for each
          desired file format.
        """

        out_FPB_cmbkit = f'{wgse.outdir.FPB}_CombinedKit'
        out_FB_cmbkit = f'{wgse.BAM.file_FB}_CombinedKit'
        out_FP = f'{wgse.outdir.FP}'

        python = f'"{wgse.python3x_FN}"'
        lifthg38 = f'"{wgse.prog_FP}hg38tohg19.py"'
        bamfile = f'"{wgse.BAM.file_FN}"'
        refgen_qFN = wgse.BAM.Refgenome_qFN

        # Select reference genome. Add liftover command to modify CombinedKit if GRCh38/HG38
        refVCFtab_qFN = wgse.reflib.get_reference_vcf_qFN(wgse.BAM.Build, wgse.BAM.SNTypeC, type="Microarray")

        if wgse.BAM.Build == 38:
            if wgse.BAM.SNTypeC == "Chr":
                liftover_hg38tohg19 = f'{python} {lifthg38} "{out_FPB_cmbkit}" hg38\n'
            elif wgse.BAM.SNTypeC == "Num":
                liftover_hg38tohg19 = f'{python} {lifthg38} "{out_FPB_cmbkit}" GRCh38\n'
            else:
                DEBUG("***Internal ERROR: Unknown reference model type")
                liftover_hg38tohg19 = ""
        else:
            liftover_hg38tohg19 = ""

        # Setup additional filenames and paths
        temp_called_vcf = f'"{wgse.tempf.FP}CombKit_called.vcf.gz"'
        temp_annotated_vcf = f'"{wgse.tempf.FP}CombKit_annotated.vcf.gz"'
        temp_result_tab = f'"{wgse.tempf.FP}CombKit_result.tab"'
        temp_sorted_result_tab = f'"{wgse.tempf.FP}CombKit_result_sorted.tab"'
        # In case DEBUG mode on and Temp directory not being cleared out; lets clear these files to make sure not reused
        wgse.tempf.list.append(temp_called_vcf) if os.path.exists(temp_called_vcf) else ""
        wgse.tempf.list.append(temp_annotated_vcf) if os.path.exists(temp_annotated_vcf) else ""
        wgse.tempf.list.append(temp_result_tab) if os.path.exists(temp_result_tab) else ""
        wgse.tempf.list.append(temp_sorted_result_tab) if os.path.exists(temp_sorted_result_tab) else ""

        temp_head_qFN = f'"{wgse.microarray_FP}raw_file_templates/head/23andMe_V3.txt"'
        ploidy_qFN = f'"{wgse.microarray_FP}ploidy.txt"'

        #samtools = wgse.samtoolsx_qFN
        bcftools = wgse.bcftoolsx_qFN
        tabix = wgse.tabixx_qFN
        sed = wgse.sedx_qFN
        sort = wgse.sortx_qFN
        cat = wgse.catx_qFN
        cmdzip = wgse.zipx_qFN
        cmdunzip = wgse.unzipx_qFN
        cpus = wgse.os_threads

        CombinedKitZIP_oFN = nativeOS(f'{out_FPB_cmbkit}.zip')

        if os.path.exists(CombinedKitZIP_oFN) and \
                os.path.getmtime(CombinedKitZIP_oFN) > os.path.getmtime(wgse.BAM.file_oFN) and \
                os.path.getsize(CombinedKitZIP_oFN) > 5000000:
            # Uncompress and reuse existing, newer CombinedKit file in the Output directory
            commands = f'{cmdunzip} -oj "{out_FPB_cmbkit}.zip" -d "{out_FP}" "{out_FB_cmbkit}.txt"\n'
            reuse = True
        else:  # Have to create a CombinedKit from scratch; use Outdir for result but delete if not requested later
            if not wgse.reflib.check_reference_genome_qFN(refgen_qFN):
                return False       # check() reports error if reference genome does not exist

            '''
              Based on the (original) Extract23 Windows script: https://github.com/tkrahn/extract23/blob/master/extract23.sh
              With parallelization extension shown in: https://gist.github.com/tkrahn/ef62cfaab678f447ea53ddee09ce0eb2 
              Original Extract23 did this only for a single 23andMe v3 template.  Here we created a CombinedKit template
              of the merger of all known formats.  Then wrote our own aconv.py module to subset the CombinedKit for each
              desired file format.
            '''
            if parallel:
                # NOTE: Parallelization method does not work on any platform (tried Win10, Ubuntu, MacOS; Intel, AMD, M1)
                #   Generates about 1% of values in error; especially chr 13-16 and 22 at start. Identical error on all
                #   platforms.
                commands = (
                    f'set +x\n'  # Turn off command echo; otherwise get echo of echo. Comments are not printed unless -v
                    f'echo "Generating CombinedKit file from BAM (takes upwards of an hour)" \n'
                    f'echo "Parallelizing each chromosome mpileup / variant-call; then will merge the results" \n'
                    f'set -x\n'
                )
                # Proper ploidy setting messed up CombinedKit generation for males; so override to avoid warning message
                #  and later error by using special ploidy.txt to give diploid for now.
                #ploidy = "GRCh38" if wgse.BAM.Build == 38 else "GRCh37" if wgse.BAM.Build == 37 or wgse.BAM.Build == 19 else "X"
                for chrom in wgse.valid_chromosomes:    # Using internal list; always same chromosome name convention
                    chrtmp = f'"{wgse.tempf.FP}{chrom}.vcf.gz"'  # Put each chromosome VCF in temp for later cleanup
                    # For GRCh numbering, need to remove "chr" from valid_chromosome name and add T to Mito name
                    fchrom = chrom.replace('chr', '').replace("M", "MT") if wgse.BAM.SNTypeC == "Num" else chrom
                    # Old versus new pileup command; old generates warning that should switch to new
                    #   {samtools} mpileup -B    -C 50 -r {fchrom} -l {reftab_qFN} -f {refgenome} -hu {bamfile}
                    #   {bcftools} mpileup -B -I -C 50 -r {fchrom} -T {reftab_qFN} -f {refgenome} -Ou {bamfile}
                    # Note -B required for Nanopore long read tools; minimal effect on paired-end sequencer output
                    commands += (
                        f'{bcftools} mpileup -B -I -C 50 -r {fchrom} -T {refVCFtab_qFN} -f {refgen_qFN} -Ou {bamfile} | '
                        f'  {bcftools} call --ploidy-file {ploidy_qFN} -V indels -m -P 0 --threads {cpus} -Oz -o {chrtmp} &\n'
                    )
                # Wait for all 25 forked processes to finish; then do a subshell cd otherwise shell line exceeds limit
                #  when using 25 fully qualified file names; all files are in temp directory and will get deleted at end
                commands += (
                    f'wait\n'
                    f'(cd {wgse.tempf.FP}; {bcftools} concat -Oz -o {temp_called_vcf} '
                    f' chr[1-9].vcf.gz chr[1-2][0-9].vcf.gz chr[M,X,Y].vcf.gz )\n'
                )
            else:
                # Old versus new pileup command; old generates warning that should switch to new
                #   {samtools} mpileup -B    -C 50 -r {fchrom} -l {reftab_qFN} -f {refgenome} -hu {bamfile}
                #   {bcftools} mpileup -B -I -C 50 -r {fchrom} -T {reftab_qFN} -f {refgenome} -Ou {bamfile}
                # Note -B required for Nanopore long read tools; minimal effect on paired-end sequencer output
                commands = (
                    f'{bcftools} mpileup -B -I -C 50 -T {refVCFtab_qFN} -f {refgen_qFN} -Ou {bamfile} | '
                    f'  {bcftools} call --ploidy-file {ploidy_qFN} -V indels -m -P 0 --threads {cpus} -Oz -o {temp_called_vcf}\n'
                )

            # Note: liftover_hg38to19 will be Null if not a Build38 BAM so simply puts an empty line when not used
            # Todo more simplication and reduction of intermediate files possible? sed/sort/cat on one line? Avoid more
            #  compression and tabix calls by piping directly?
            commands += (
                f'{tabix} -p vcf {temp_called_vcf}\n'
                f'{bcftools} annotate -Oz -a {refVCFtab_qFN} -c CHROM,POS,ID {temp_called_vcf} > {temp_annotated_vcf}\n'
                f'{tabix} -p vcf {temp_annotated_vcf}\n'
                f'{bcftools} query -f \'%ID\t%CHROM\t%POS[\t%TGT]\n\' {temp_annotated_vcf} -o {temp_result_tab}\n'
                f'{sed} \'s/chr//; s/\tM\t/\tMT\t/g; s/\///; s/\.\.$/--/; s/TA$/AT/; s/TC$/CT/; s/TG$/GT/; s/GA$/AG/; '
                f'    s/GC$/CG/; s/CA$/AC/\' {temp_result_tab} | {sort} -t $\'\t\' -k2,3 -V > {temp_sorted_result_tab}\n'
                f'{cat} {temp_head_qFN} {temp_sorted_result_tab} > "{out_FPB_cmbkit}.txt"\n'
                f'{liftover_hg38tohg19}'
                f'{cmdzip} -j "{out_FPB_cmbkit}.zip" "{out_FPB_cmbkit}.txt"\n'
            )
            reuse = False
        run_bash_script("ButtonCombinedKit", commands, parent=parent_window)
        # Todo really should check if completed OK; but general issue with all run_bash_script calls

        return reuse


def cancel_autosomal_formats_window():
    ''' Remove Autosomal formats window in preparation for restoring WGSE Main Window '''
    from mainwindow import resume_main_window   # Have to localize due to import loop

    global selectAutosomalFormatsWindow

    try:        # May have been destroyed before this call
        selectAutosomalFormatsWindow.destroy()
    except:
        pass
    if __name__ != '__main__':
        resume_main_window()    # sets up for return to main window


# ***************MAIN PROGRAM (stand-alone)*******************
# Todo standalone program in development; not really been tested of late
if __name__ == '__main__':
    from mainwindow import button_select_BAM_file, button_select_output_path       # Need local import to avoid loop

    """ Microarray start as an independent task (main program) """
    wgse.init(False)

    # Could potentially just get parameters from stored settings JASON file? Or make args optional as can get from there
    if not button_select_BAM_file(sys.argv[1]):
        exit()      # Already reported issue in pop-up
    button_select_output_path(sys.argv[2])

    button_select_autosomal_formats()   # What gets called by WGS Extract mainwindow.py to enter here

    # Simply exit out which closes the iconified mainWindow that was never filled by setup_mainWindow call


