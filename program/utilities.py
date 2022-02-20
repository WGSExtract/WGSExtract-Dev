# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
  General Utilities module support for the WGS Extract system

   This file is included by all others; import
  So imports of library functions are buried in function definitions where needed.
   (None are called often so there is little penalty for doing so.)

  Classes: TemporaryFiles, LanguageStrings
"""

import settings as wgse
font = wgse.font


def DEBUG(msg):
    """ Old C Macro hack for programmer debug system; but remains as executed code in Python """
    if wgse.DEBUG_MODE:                  # Could use __debug__ but does not save much as calls have if statement then
        print(f"DEBUG: {msg}")


def wgse_message(type, title, bodyx, body):
    from tkinter import messagebox

    ttitle = wgse.lang.i18n[title]
    tbody  = body if bodyx else wgse.lang.i18n[body]    # Is body translated yet? If not, then translate

    print(tbody.replace("\n\n", "\n"))

    if type == "error":
        messagebox.showerror(ttitle, tbody)
    elif type == "warning":
        messagebox.showwarning(ttitle, tbody)
    elif type == "info":
        messagebox.showinfo(ttitle, tbody)
    # Have not gotten to sophistication of question / response yet; but likely will


# Technically, Cygwin BASH will accept both forms.  And in fact, WinPython returns Win10 disk drive letters to
# start a path. Which most programs (CygWin bash included) seem to still accept even if they require a forward slash.
def nativeOS(path_FP):
    """
    Replace os_slash's on path; from universal forward to back slash if Win10 """
    if wgse.os_plat == "Windows" and path_FP:
        if path_FP[0:10] == "/cygdrive/":           # Just in case
            path_FP = f'{path_FP[11].upper()}:{path_FP[12:]}'    # Convert to Win10 Drive letter format
        elif path_FP[0:4] == "/mnt/":               # In case used WSL mechanism
            path_FP = f'{path_FP[5].upper()}:{path_FP[6:]}'      # Convert to Win10 Drive letter format
        elif path_FP[1] == ":":                    # Often drive letters still come through on Win10 systems
            pass
        elif path_FP[0:1] == "//":                 # Network drive
            pass
        path_oFP = path_FP.replace("/", "\\")       # Change to backward slash
    else:
        path_oFP = path_FP      # Nothing to do for Linux, Unix, MacOS
    #DEBUG(f'nativeOS: before {path_FP}, after {path_oFP}')
    return path_oFP


def universalOS(path_oFP, wsl_fix=False):
    """
    From universal to OS Specific (forward slash to back slash in Win10).
    Is a no-op except for forward / back slash change due to inconsistencies between WinPython and CygWin
    wsl bwa is the only place that will NOT accept win10 drive letters; locally patch it there
    """
    if wgse.os_plat == "Windows":
        path_FP = path_oFP.replace("\\", "/")       # Change to forward slash
        if path_FP[1] == ":":                      # Win10 Drive letter
            if wsl_fix:
                path_FP = f'/mnt/{path_FP[0].lower()}{path_FP[2:]}'  # Change driver letter to path
            else:
                #path_FP = f'/cygdrive/{path_FP[0].lower()}{path_FP[2:]}'   # Change driver letter to path
                # Leave as is; most programs do not follow or understand Cygwin introduced mount points
                # and Cygwin Bash understands windows disk letter specifiers
                pass
        elif path_FP[0:1] == "//":                  # Network drive; leave alone
            pass
    else:
        path_FP = path_oFP                          # Nothing to do for Linux, Unix, MacOS
    #DEBUG(f'universalOS: before {path_oFP}, after {path_FP}')
    return path_FP


def unquote(path):
    """ Simple routine to remove quotes surrounding a path; instead of simply calling .replace('"', '') """
    return path.replace('"', '')


# Todo Still needed?  Historical, before quoted path names everywhere.  Added space as legal character now.
def is_legal_path(path_to_check):
    """
    Check file path string for illegal characters. Space no longer illegal as have quoted paths everywhere.
    (but with quoted paths is this still needed?)
    """
    import re
    return re.match(r'^[a-zA-Z0-9_/\-:~.\s\\]+$', path_to_check)


class Error(Exception):
    pass

class Warning(Exception):
    pass


######################################################################################################################
# Todo Temp Location should become part of a settings class; that also saves / restores settings
class TemporaryFiles:
    """
    A class to manage temporary files.
    Includes self.list for list of files and folders to clean; likely in the user specified output directory
    Optionally clears out user specified temp directory of all files.

    As of 1 Apr 2021, we have modified to use the os.getpid() of the process instance to create a subdirectory
    within the specified Temporary Files directory.  That way, multiple instances of the program can be safely
    started and run.
    """

    def __init__(self, default_FP, top_level=False):

        self.set  = None
        self.FP   = None
        self.oFP  = None
        self.FB   = None            # Now will always be the PID

        self.default_FP = None      # For the root "default"

        self.change(default_FP)     # Initial call uses this value to define default_FP

        self.list = []  # files to be cleaned that are not inside TempFiles area (helps limit os.remove calls)
        if self.FP:
            self.clean(True) if top_level else ""  # Start with a clean temporary files directory
        else:
            # Setting temporary files directory failed; really fatal and should probably exit (raise exception)
            DEBUG(f'***ERROR: Cannot set default Temporary Files directory: {default_FP}')
            del self


    def change(self, dir_FP):       # For when changing temporary files directory; including initialization
        """
        The real guts of the Temporary Files Directory setup.

        Note: no longer require the directory be empty.  We will setup a subdirectory based on the PID which will,
        by definition, be empty when we start.  OK if other files in the directory to start. This way we now support
        multiple instances of the program at the same time.
        """
        import os

        if not dir_FP or dir_FP == "/":
            # Directory must exist and not be the root file system; universal naming here
            return      # Simply return and ignore as null may be due to restore settings or user cancel

        dir_FP += '/' if dir_FP[-1] != '/' else ''  # Assure trailing slash; universalOS so always forward slash
        dir_oFP = nativeOS(dir_FP)

        # If just initializing default, do that. Otherwise, check for common conflicts and validity of temp files dir
        if not self.default_FP:     # If just initializsing, default is not set yet
            self.default_FP = dir_FP
        elif dir_FP == self.FP:     # No change; nothing to do. Silently return. If initilization, will not be set
            return

        if not (is_legal_path(dir_FP) and os.path.isdir(dir_oFP)) or \
           dir_FP in [wgse.install_FP, wgse.prog_FP] or \
           (wgse.reflib and (dir_FP == wgse.reflib.default_FP or dir_FP == wgse.reflib.FP)) or \
           (wgse.BAM and dir_FP == wgse.BAM.file_FP) or \
           (wgse.outdir and dir_FP == wgse.outdir.FP):
            # Ouch, do not allow the reference library for temporary files else we wipe it out
            # do not allow BAM source or output directory areas (if set) for similar reasons
            # Or maybe simply a bad path (not a directory, does not exist, etc)
            # Not setting default so must be from mainwindow button to set TempFile area; so tkinter initialized
            # Trying to be extra careful as we delete any files or directories in the Temporary Files directory
            DEBUG(f'Bad Temporary File Path: {dir_FP}')
            wgse_message("error", 'InvalidTempDirTitle', False, 'errTempDirPath') if wgse.lang.i18n else ""
            return
        else:
            # Want to start with an empty temporary files directory to make sure not the wrong area we wipe out
            # So warn user if not empty. Todo create new subdirectory based on PID so unique and empty.
            dir_contents = [x for x in os.listdir(dir_oFP) if not x.startswith('.')]
            if len(dir_contents) > 0:
                DEBUG(f'New Temporary File area already has content (# files: {len(dir_contents)})\n{dir_contents}')
                wgse_message("error", 'InvalidTempDirTitle', False, 'errTempDirNotEmpty') if wgse.lang.i18n else ""

        # OK, finally, everything is good and we can set the values
        self.set = True if dir_FP != self.default_FP else False
        self.FP  = dir_FP
        self.oFP = dir_oFP
        self.FB  = os.path.basename(dir_FP[:-1])    # Trick; remove trailing slash so basename returns directory name
        # self.list[] is independent of the temporary files directory setting


    def clean(self, incl_temp_dir=True):        # To clean the list[] and temp dir
        import os       # for os.path.*, os.remove, os.listdir
        import shutil   # shutil.rmtree()

        # First, cleanup files on the clean-up list (ever specify a directory? maybe remove)
        for oFN in self.list:
            DEBUG(f"Clean Temp List has {len(self.list)} files or directories to clean:\n{self.list}")
            try:
                if os.path.isdir(oFN):
                    oFN += wgse.os_slash
                    shutil.rmtree(oFN) if not wgse.DEBUG_MODE else ''  # Careful: remove hieararchical tree
                elif os.path.isfile(oFN):
                    os.remove(oFN) if not wgse.DEBUG_MODE else ''
                else:
                    DEBUG(f"Cleanup object is not a file or directory? {oFN} (ignoring)")
            except:
                pass
        self.list.clear()

        if not incl_temp_dir:
            return

        if not (self.oFP and os.path.isdir(self.oFP)):  # Just being extra careful; directory path not set or invalid
            DEBUG(f"***FATAL ERROR: No valid Temporary Directory set; inside wgse.tempf.clean(): {self.oFP}")
            wgse_message("error", 'InvalidTempDirTitle', False, 'errTempDirPath')
            exit()

        # Now clean-up the temp directory no matter what was listed in the cleanup global
        #  but do not clean temp directory if DEBUG mode is on -- want to save files for analysis
        for oFBS in os.listdir(self.oFP):
            oFN = self.oFP + oFBS  # make absolute
            try:
                if os.path.isdir(oFN):  # yleaf tempYleaf dirctory is only commonly occuring directory
                    oFN += wgse.os_slash
                    DEBUG(f"Clean_Temp dir: {oFN} (DEBUG_MODE = {wgse.DEBUG_MODE})")
                    shutil.rmtree(oFN) if not wgse.DEBUG_MODE else ''  # Must be careful; hierarchical tree removal
                elif os.path.isfile(oFN):
                    DEBUG(f"Clean_Temp file: {oFN} (DEBUG_MODE = {wgse.DEBUG_MODE})")
                    os.remove(oFN) if not wgse.DEBUG_MODE else ''
                else:
                    DEBUG(f"Cleanup type? {oFN} (ignoring)")
            except:
                pass


######################################################################################################################
class LanguageStrings:
    """
      Internationalization subsystem providing language specific labels, buttons and messages via table lookup.

      Entry Points:
        __init__ (Class creation entry; takes in language.csv file spec. Read at class init)
        change_language  (waterfall call from read_language_setting; also from mainWindow language button)
        switch_language  (to set or switch a language mid-stream; creates new mainWindow if language change)

        In DEBUG_MODE, for language translators, there is a reload languages.csv button to recreate this class instance
         and start with change_language_setting again (like at program start). Purposely do not save set language
         so query pop-up happens again to help verify new languages.csv file loaded correctly.

      Note:
        Moved get_lanugage GUI to to mainwindow.py button_set_language (keep all GUI stuff there)
        Moved store / read language setting to generalized functions in settings module (added more settings)
        Greatly simplified this class then.
    """

    def __init__(self, language_oFN):
        """ Read in file language_oFN as init action. Is a CSV dictionary file of language translations
            wait till actually read_language_setting to then set language and active wgse.lang.i18n dictionary. """
        import openpyxl

        # Start by reading in the language.csv file to determine languages supported and their translations
        self._lang = {};   self.request = {}
        try:
            wb = openpyxl.load_workbook(language_oFN)
            ws = wb.active
            dim = ws.calculate_dimension()
            DEBUG(f'Language.xlsx Dimension: {dim}')

            first = True  ;  avail_langs = []
            for row in ws.iter_rows(values_only=True):
                if first:       # To get header row; does not seem to be a next in this package
                    num_of_lang = int(row[0])
                    avail_langs = [row[i+1] for i in range(num_of_lang)]    # iterate over # langs to read language names / key
                    self._lang = {key: {} for key in avail_langs}   # iterate langs to create empty translation dictionaries
                    first=False
                else:
                    # Rest of file is language translations with the Key in column 1 and each native string following
                    # A special entry is RequestLang that has language specific string to ask for a language setting
                    # _lang is a dictionary of dictionaries. Language as first key. First column in each row as next key.
                    # Final dictionary has string representation (i18n) of the second key identifier.

                    k = row[0]  # initial element (col 0) is key index into language dictionary
                    for i in range(0, len(avail_langs)):
                        self._lang[avail_langs[i]][k] = row[i + 1].replace("^", "\n") if row[i + 1] else ""

            # Directory with language request strings for get_language pop-up; no reliance on internal _lang{}{}
            for i in range(0, len(avail_langs)):
                self.request[avail_langs[i]] = self._lang[avail_langs[i]]["RequestLanguage"]
        except:
            print(f"***Error reading from system Language Strings file: {language_oFN}")
            raise ValueError
        self.language_oFN = language_oFN
        self.avail_langs = avail_langs
        # self.request = request
        # self._lang = _lang
        self.language = None    # Keep null; language not set yet; just the language translations
        self.i18n = {}          # only set to particular translation dictionary once language is selected
        self.selectLanguageWindow = None        # Used to save language window popup to cancel from here


    def change_language(self, language=None):
        """
        Called after loading settings from user file. If not set or not valid, query user.
        Local import to avoid loop. All GUI pushed to mainwindow module.
        """
        from mainwindow import button_set_language

        if language not in self.avail_langs:      # If requested language not set or not understand
            button_set_language()           # Ask user for language; calls switch_language() directly
        else:
            self.switch_language(language)  # Switch language


    def switch_language(self, language):
        """ Language Subsytem function to set/switch language; externally called to change language. """
        from mainwindow import init_mainWindow, setup_mainWindow, resume_main_window

        # If got here from the Language Selection button; destroy its pop-up window (ignore errors)
        # We saved the window ID in the class' local variables
        if self.selectLanguageWindow:
            try:
                self.selectLanguageWindow.destroy()
            except:
                pass
            self.selectLanguageWindow = None

        if language == self.language:             # If no change, then simply return
            pass
        elif language not in self.avail_langs:    # If not valid language to switch too ...
            DEBUG(f'*** INTERNAL ERROR: Unknown language: "{language}"; ignoring')
        else:                                       # Everything OK; setup new language translation for selected lang
            self.i18n = self._lang[language]  # Set wgse.lang.i18n to selected language dictionary
            self.language = language          # Save current language in class
            DEBUG(f"New Language: {language}")

            if wgse.window and wgse.dnaImage:       # If setup main window already, then recreate with new language
                wgse.dnaImage = None                # Only way we know setup_mainWindow was run
                wgse.window.destroy()               # Kill current mainwindow
                wgse.window = init_mainWindow()     # Create new mainWindow in new language (what about mainloop?)
                setup_mainWindow()                  # Fill in new mainWindow (sets wgse.dnaImage)

        if wgse.window and wgse.dnaImage:
            resume_main_window()


######################################################################################################################
class OutputDirectory:
    """
    A class to handle the Output Directory setup.  Mostly for changing it in a consistent way from the users button
    as well as the stored setting.  And to move the functionality out of mainwindow so it can focus on UI.
    """

    def __init__(self):
        self.oFP   = None  # Output directory path (user specified)
        self.FP    = None  # Output directory path (user specified) (Unix/Universal format)
        self.oFPB  = None  # Output directory path with BAM basename (both user specified)
        self.FPB   = None  # Output directory path with BAM Basename (both user specified) (Unix/Universal format)
        self.FB    = None  # Output directory path (user specified) basename (not BAM but last directory name in path)

    def change(self, new_FP):
        """
          Called to change (or initially set) the Output Directory.
        """
        import os

        if new_FP == self.FP:
            return      # No change; simply return. Especially if stll None

        new_oFP = nativeOS(new_FP)

        # If a valid and legal output path; then set the global variables for this new setting
        if is_legal_path(new_FP) and os.path.exists(new_oFP) and os.path.isdir(new_oFP):
            self.FP = new_FP
            self.oFP = new_oFP
            DEBUG(f"Output Path (final): {self.FP}")
            if wgse.BAM:
                self.FPB  = new_FP  + wgse.BAM.file_FB
                self.oFPB = new_oFP + wgse.BAM.file_FB
            # Like the Unix basename() for directories; save last name in path
            self.FB = new_FP.split('/')[-2] if len(new_FP.split('/')) > 1 else new_FP