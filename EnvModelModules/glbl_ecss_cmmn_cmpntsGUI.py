"""
#-------------------------------------------------------------------------------
# Name:        glbl_ecss_cmmn_funcs.py
# Purpose:     module comprising miscellaneous functions common to some Global Ecosse variations
# Author:      Mike Martin
# Created:     31/07/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
"""

__prog__ = 'glbl_ecss_cmmn_funcs.py'
__version__ = '0.0.0'

# Version history
# ---------------
#
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QLabel, QComboBox

from shape_funcs import calculate_area, format_bbox

sleepTime = 5

WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '
DUMMY_STR = 'dummy'

GRANULARITY = 120

WDGT_SIZE_60 = 60

RESOLUTIONS = {120:'30"', 30:'2\'', 20:'3\'', 10:'6\'', 8:'7\' 30"', 6:'10\'', 4:'15\'', 3:'20\'', 2:'30\''}
REVERSE_RESOLS = {}
for key in RESOLUTIONS:
    REVERSE_RESOLS[RESOLUTIONS[key]] = key

def grid_resolutions(form, grid, irow):
    """

    """
    irow += 1
    lbl16 = QLabel('Grid resolution:')
    lbl16.setAlignment(Qt.AlignRight)
    helpText = 'The size of each grid cell is described in arc minutes and arc seconds. The smallest cell resolution \n' \
        + 'corresponds to that of the HWSD database (30 arc seconds) and the largest to that used by the climate data ' \
        + '(30 arc minutes)'
    lbl16.setToolTip(helpText)
    grid.addWidget(lbl16, irow, 0)

    combo16 = QComboBox()
    for resol in sorted(RESOLUTIONS,reverse = True):
        combo16.addItem(str(RESOLUTIONS[resol]))
    combo16.setToolTip(helpText)
    combo16.setFixedWidth(WDGT_SIZE_60)
    form.combo16 = combo16
    combo16.currentIndexChanged[str].connect(form.resolutionChanged)
    form.combo16 = combo16
    grid.addWidget(combo16, irow, 1)

    form.lbl16a = QLabel('')
    form.lbl16a.setToolTip(helpText)
    grid.addWidget(form.lbl16a, irow, 2, 1, 3)

    return irow

def calculate_grid_cell(form, granularity = GRANULARITY):
    """

    """

    # use current lower left latitude for reference
    # =============================================
    latitude = 52.0     # default
    if hasattr(form, 'w_ll_lat'):
        latText = form.w_ll_lat.text()
        try:
            latitude = float(latText)
        except ValueError:
            print(latText)

    resol = form.combo16.currentText()
    try:
        granul = REVERSE_RESOLS[resol]
    except KeyError as err:
        print(str(err))
        return

    resol_deg = 1.0/float(granul)  # units of decimal degrees
    bbox = list([0.0, latitude, resol_deg, latitude + resol_deg])
    area = calculate_area(bbox)
    form.lbl16a.setText(format_bbox(bbox,area,2))

    form.req_resol_upscale = int(granularity/granul)    # number of base granuals making up one side of a cell
    form.req_resol_granul = granul                      # number of cells per degree
    form.req_resol_deg = resol_deg                      # size of each trapizoid cell in decimal degrees

    return resol_deg

def print_resource_locations(setup_file, config_dir, hwsd_dir, wthr_dir, lta_nc_fname, sims_dir, log_dir):
    """
    report settings
    """
    print('\nResource locations:')
    print('\tsetup file:          ' + setup_file)
    print('\tconfiguration files: ' + config_dir)
    print('\tHWSD database:       ' + hwsd_dir)
    print('\tweather data sets:   ' + wthr_dir)

    if lta_nc_fname is not None:
        print('\tLTA weather file:    ' + lta_nc_fname)

    print('\tsimulations:         ' + sims_dir)
    print('\tlog_dir:             ' + log_dir)
    print('')

    return