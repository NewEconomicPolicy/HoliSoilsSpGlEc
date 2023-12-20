#-------------------------------------------------------------------------------
# Name:        run_mode_params.py
# Purpose:     Class to organise run parameters
# Author:      Mike Martin
# Created:     28/09/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'run_mode_params.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import lexists
from os import makedirs
from numpy import arange, dtype, zeros, int32, int64

class dummy(object,):

    def __init__(self, form, clim):

        self.lgr = form.lgr
        self.sims_dir = form.sims_dir
        self.fut_clim_scen = clim.fut_clim_scen
        self.study = form.study
        self.fstudy = form.fstudy
        self.kml_flag = form.kmlFlag
        self.default_model_switches = form.default_model_switches
        self.req_resol_deg = form.req_resol_deg
        self.hist_ave_file = clim.hist_ave_file
        self.clim_dir = form.sims_dir + '/' + clim.fut_clim_scen + '/'
        if not lexists(self.clim_dir):
            makedirs(self.clim_dir)