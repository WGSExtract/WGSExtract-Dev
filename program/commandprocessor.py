# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
    Module commandprocessor.  Main handling of subprocess calls -- currently all are bash shell runs.
    Will eventually expand to async multiprocessor handling with a list of processes.  To allow queue'ed calls and
    the main program window to be non-blocking.  Handles new ProcessWait pop-up with timer to terminate errant
    subprocess calls. Still not terminate for running call though.

    v1 /v2 was 50% BATCH and 50% BASH calls.  Moving to all BASH, code is more OS
    independent (as long as Win10 environment has BASH/Unix calls available) and internal file / path handling easier
    as most are left in universal/Unix format and not nativeOS specific. Only the BASH script filename itself needs
    to be in nativeOS form. Historically, BATCH files were needed as dos2unix and similar calls were made.  But all
    dos2unix needs were replaced by simply writing binary ('wb') and using .encode() on strings written in binary
    form; alleviating the need to convert back to common POSIX \n format.
"""

import os           # os.chmod()
import stat         # stat.S_IXxxx flags
import time         # time.time()
import subprocess   # Popen, run, etc
# from tkinter.ttk import Button, Label
from tkinter import Toplevel, Label

import settings as wgse
font = wgse.font
from utilities import DEBUG, wgse_message
# from mainwindow import resume_main_window     # embedded in button_continueAfterBatchJob to break loop


####################################################################################################################
# Command Execution Subsection

pleaseWaitWindow = None


def is_command_available(command, opts, internal=True, min_version=0):
    """ Verify a shell command is available to run as specified.  Possibly get version to check as well. """
    try:
        result = subprocess.run([command, opts], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except:
        if not internal:
            wgse_message("warning", 'NoCmdWindowTitle', True,
                         wgse.lang.i18n['NoCmdInstalledErrorMessage'].replace("{{cmd}}", command))
        return False
    else:
        if " 11 " in result.stdout.decode('utf-8'):
            # Todo put in version check; not just check if it is available
            pass
        return True


def time_label(etime):
    time_exp = (etime / 3600., 'hours')      if etime > 3600 else \
                (int(etime / 60), 'minutes') if etime > 60 else \
                (etime, 'seconds')
    return f"{time_exp[0]:3.1f} {time_exp[1]}" if time_exp[1] == 'hours' else f"{time_exp[0]:3d} {time_exp[1]}"


def run_external_program(script_FBS, command):
    """ Run an external batch program with a time-limit"""

    # Rough std deviation (max) from expected time (x2); more robust for table out-of-date
    maxtime = wgse.expected_time.get(script_FBS,3600) * 2
    command_str = os.path.basename(command[-1])     # " ".join(command).strip()
    start_time = time.time()
    start_ctime = time.ctime()
    print(f'--- Exec: {command_str}, started @ {start_ctime}')     # REH 14Mar2020 Debug
    try:
        # subprocess.run(command_to_run_and_args, timeout=maxtime*2)
        with subprocess.Popen(command) as p: # , stderr=subprocess.STDOUT, stdout=subprocess.PIPE, text=True) as p:
            # stdout, stderr = p.communicate(timeout=maxtime)
            stdout, stderr = p.communicate()
    except subprocess.TimeoutExpired:   # todo This specific exception not being caught here?
        # todo need to see if due to a timeout and handle appropriately; also user initiated abort
        p.kill()
        stdout, stderr = p.communicate()
        stop_ctime = time.ctime()
        tlabel = time_label(maxtime)
        DEBUG(f"--- FAILURE to finish before timeout {tlabel}: {command_str} (@ {stop_ctime})----")
        raise
    except subprocess.CalledProcessError:
        # todo need to handle different process termination conditions
        DEBUG(f"--- FAILURE to execute: {command_str} ----")
        raise
    else:
        stop_ctime = time.ctime()
        tlabel = time_label(round(time.time() - start_time))
        print(f'--- SUCCESS: {tlabel} to run: {command_str} (finished @ {stop_ctime}')


def show_please_wait_window(reason, parent):
    global pleaseWaitWindow

    # Make code more robust to errors; external language file and expected time table may not be up-to-date
    etime = wgse.expected_time.get(reason, 0)
    # etime = round(min(1.5, etime * wgse.BAM.relfsize)) # See note below on adjusting time relative to 60 GB nominal
    tlabel = time_label(etime) if etime else f"unknown ({reason}?)"
    verbose_reason = wgse.lang.i18n.get(reason, f"(Internal Translation Error: {reason} unknown)")
    DEBUG(f'In Please Wait: {verbose_reason}, Expected Wait: {tlabel}, Start')
    """ 
        We will adjust expected wait time by amount relative to nominal 45 GB 30x WGS file size. Floor of 0.01 adj.
        Using linear but likely not linear.  Idea is smaller, subsetted BAMs may take less time? True for operations
        that have to scan the whole file but not others?  Need to study this further and drop if nonsense.
        Does not seem to work as written; commenting out for now
    """

    pleaseWaitWindow = Toplevel(wgse.window)
    pleaseWaitWindow.transient(parent)
    pleaseWaitWindow.title(wgse.lang.i18n['PleaseWait'])
    pleaseWaitWindow.geometry("")
    pleaseWaitWindow.protocol("WM_DELETE_WINDOW", cancelWait)
    # pleaseWaitWindow.attributes('-topmost', 'true')
    # pleaseWaitWindow.parent.withdraw()
    pleaseWaitWindow.iconbitmap(wgse.icon_oFP) if wgse.os_plat == "Windows" else ''
    pleaseWaitWindow.columnconfigure(0, weight=1)
    pleaseWaitWindow.rowconfigure(0, weight=1)
    Label(pleaseWaitWindow, text=str(wgse.lang.i18n['PleaseWait']), font=font['28b']).grid(column=0, row=0, padx=56, pady=28)
    Label(pleaseWaitWindow, text=verbose_reason, font=font['14']).grid(column=0, row=1, padx=1, pady=1)
    Label(pleaseWaitWindow, text=f'{wgse.lang.i18n["ExpectedWait"]} {tlabel}.').grid(column=0, row=2, padx=1, pady=1)
    time_now = time.ctime()
    Label(pleaseWaitWindow, text=f'{wgse.lang.i18n["StartedAt"]} {time_now}').grid(column=0, row=3, padx=1, pady=1)
    # TODO add user abort button that is caught, kills job waiting on, and takes one back where?
    pleaseWaitWindow.update()
    pleaseWaitWindow.grab_set()


def cancelWait():
    global pleaseWaitWindow

    #DEBUG("cancelWait; unfortunately not yet implemented!")
    # TODO Setup exception that can be caught by threads executing pleaseWait
    # noinspection PyUnresolvedReferences
    pleaseWaitWindow.destroy()


def abortWait():
    # Todo implement abortWait()
    #DEBUG("abortWait called; not yet implemented.")
    cancelWait()


def finishWait():
    # Todo implement finishWait()
    #DEBUG("finishWait called; not yet implemented.")
    cancelWait()


def run_bash_script(script_title, script_contents, parent=wgse.window):
    """
        Main entry point for module.  Create the BASH script file in Temp from the supplied content, then
        run the bash job directly.
    """
    script_oFN = f'{wgse.tempf.oFP}{script_title}.sh'
    if os.path.isfile(script_oFN):             # Todo really needed? Won't we just overwrite if still exists?
        os.remove(script_oFN)
    with open(script_oFN, "wb") as f:          # 15 Mar 2020 REH Changed to binary write to keep Unix \n format
        f.write("#!/bin/bash -x\n".encode())
        f.write(script_contents.encode())
    os.chmod(script_oFN, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # Historically added chmod and dos2unix to command here

    DEBUG(f"Starting Bash Script: {script_oFN}")
    show_please_wait_window(script_title, parent)
    # full path specified so do not need to pre-pend './'; Windows needs BASH command to start .sh files
    # All commands are BASH script files; never a user command directly. So no need to parse script content using shlex
    command = [wgse.bashx_oFN, "-x", script_oFN] if wgse.os_plat == "Windows" else [script_oFN]
    run_external_program(script_title, command)
    finishWait()