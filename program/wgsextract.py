#!/usr/bin/env python3
# coding: utf8
# Copyright (C) 2018-2020 City Farmer
# Copyright (C) 2020-2021 Randy Harr
#
# License: GNU General Public License v3 or later
# A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

"""###################################################################################################################
  Main WGS Extract program module. As can see, is pretty devoid of effort.  Simply gets global system setup via
  settings module. Then calls the mainWindow processing subsystem setup and loop to wait on a user button input.
  Included by all other modules that need settings. Hence kept simple as well with local imports inside wgse_main()
"""
from argparse import ArgumentParser


def get_arguments():
    """
    For stand-alone, command line invocation (or even double click in GUI after associating BAM/CRAM file
    types with program)
    Todo implement this feature; template / fragment now
    Todo processing if standalone file name without other args
    Todo help args need internationalization
    :return:
    """
    parser = ArgumentParser()

    parser.add_argument("-b", "--bam",
            dest="Bamfile", required=False,
            help="input BAM file", metavar="PATH")

    parser.add_argument("-c", "--cram",
            dest="Cramfile", required=False,
            help="input CRAM file", metavar="PATH")

    parser.add_argument("-d", "--debug",
            dest="DebugMode", required=False,
            help="Turns on Debug Mode")

    parser.add_argument("-f", "--fasta-ref",
            dest="reference", required=False,
            help="fasta reference genome sequence ", metavar="PATH")

    parser.add_argument("-o", "--outputdir",
            dest="Outputdir", required=False,
            help="Folder name containing outputs", metavar="STRING")

    parser.add_argument("-t", "--Threads",
            dest="threads", required=False, type=int, default=2,
            help="Set number of additional processor threads")

    args = parser.parse_args()
    return args


def wgse_main(module):
    """ WGS Extract main program start as independent task. Gives a direct call. """
    import settings as wgse
    from mainwindow import setup_mainWindow

    # Start WGSE subsystem with all the main program settings (settings.py)
    #   Used to be a class; now just module to remove an additional layer of naming
    wgse.init(True)     # Includes subsystem init calls including init_mainWindow

    if module == '__main__':
        # Print explanation to command script window, setup main window, let's go
        print(wgse.lang.i18n['ExplainWhyTerminalWinIsOpen'])
        wgse.window = setup_mainWindow()
        wgse.window.mainloop()  # Go wait for a button click on the main window ...


# ***************MAIN PROGRAM*******************
if __name__ == '__main__':

    args = get_arguments()      # Place holder; not implemented yet. Intent is to allow BAM file in command line.
    wgse_main(__name__)


"""
  To use the global settings in other modules, we need the following code:
    import settings as wgse
  
  Instead of a class Settings in settings.py, we instead use a simple module settings.py with all the global 
   variables defined and initialized at the top level there. Then simply define an init() function that is called
   from the wgsextract main area. Using a class gives us an extra level of naming indirection.

  Both techniques still require the additional line:
    font = wgse.font
  which exists to get a global variable without any context / module path.  Simply shorten the common name in the GUI.
  Works because is a static value.
"""