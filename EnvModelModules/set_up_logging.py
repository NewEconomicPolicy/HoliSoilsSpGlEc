"""
#-------------------------------------------------------------------------------
# Name:        set_up_logging.py
# Purpose:     script to read read and write the setup and configuration files
# Author:      Mike Martin
# Created:     31/07/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
"""

__prog__ = 'set_up_logging.py'
__version__ = '0.0.0'

# Version history
# ---------------
#
from os import makedirs, remove
from os.path import isdir, join, normpath, getmtime

import logging
from time import strftime, sleep
from glob import glob
try:
    from PyQt5.QtGui import QTextCursor
except ModuleNotFoundError as err:
    print(str(err))

ERROR_STR = '*** Error *** '
MAX_LOG_FILES = 10
sleepTime = 5

def set_up_logging(form, appl_name, lggr_flag = False):
    """
    set up logging
    """
    if not hasattr(form, 'settings'):
        print(ERROR_STR + 'could not setup logging, form object must have settings attribute')
        return

    if 'log_dir' not in form.settings:
        print(ERROR_STR + 'could not setup logging, settings dictionary must have log_dir attribute')
        return

    log_dir = form.settings['log_dir']
    if not isdir(log_dir):
        makedirs(log_dir)

    # string to uniquely identify log files
    # =====================================
    date_stamp = strftime('_%Y_%m_%d_%I_%M_%S')
    log_file_name = appl_name + date_stamp + '.log'


    log_fname = join(log_dir, log_file_name)

    # Setup up initial logger to handle logging prior to setting up the full logger using config options
    # ==================================================================================================
    lggr = logging.getLogger(appl_name)
    lggr.setLevel(logging.INFO)

    fhand = logging.FileHandler(log_fname)    # send log recs (created by loggers) to appropriate destination

    fhand.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')    # alternatively: Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fhand.setFormatter(formatter)                   # specify layout of records
    lggr.addHandler(fhand)

    # if more than 15 or so files then delete oldest
    # ==============================================
    log_flist = glob(normpath(log_dir + '/' + appl_name + '*.log'))
    list_len = len(log_flist)
    num_to_delete = list_len - MAX_LOG_FILES
    if num_to_delete > 0:
        log_flist.sort(key = getmtime)    # date order list
        for ifile in range(num_to_delete):
            try:
                remove(log_flist[ifile])
                lggr.info('removed log file: ' + log_flist[ifile])
            except (OSError, IOError) as err:
                print('Failed to delete log file: {}'.format(err))

    if lggr_flag:
        form.lggr = lggr
    else:
        form.lgr = lggr

    return

class OutLog:
    def __init__(self, edit, out=None, color=None):
        """
        can write stdout, stderr to a QTextEdit widget
        edit = QTextEdit
        out = alternate stream (can be the original sys.stdout)
        color = alternate color (i.e. color stderr a different color)
        """
        self.edit = edit
        self.out = None
        self.color = color

    def flush(self):

        pass

    def write(self, mstr):
        if self.color:
            tc = self.edit.textColor()
            self.edit.setTextColor(self.color)

        self.edit.moveCursor(QTextCursor.End)
        self.edit.insertPlainText( mstr )

        if self.color:
            self.edit.setTextColor(tc)

        if self.out:
            self.out.write(mstr)
