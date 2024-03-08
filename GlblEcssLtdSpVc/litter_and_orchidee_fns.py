#-------------------------------------------------------------------------------
# Name:        litter_and_orchidee_fns.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------

__prog__ = 'litter_and_orchidee_fns.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import exists, isfile
from math import isnan
from netCDF4 import Dataset, num2date
from pandas import read_excel, DataFrame

from mngmnt_fns_and_class import OchideeSet, check_ochidee_dset

ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '
CNVRSN_FACT = 10000 * 365 / 1000  # converts gC/m**2/day to kgC/ha/yr

def resize_yrs_pi(sim_strt_yr, sim_end_yr, yrs_pi):
    """
    patch to enable adjust yrs_pi to correspond to user specified simulation period
    """
    if yrs_pi is None:
        return None

    yr_frst = yrs_pi['yrs'][0]
    yr_last = yrs_pi['yrs'][-1]

    pi_frst = yrs_pi['pis'][0]
    pi_last = yrs_pi['pis'][-1]

    sim_yrs = list(range(sim_strt_yr, sim_end_yr + 1))

    sim_pis = []
    for iyr, yr in enumerate(sim_yrs):
        if yr < yr_frst:
            sim_pis.append(pi_frst)
        elif yr > yr_last:
            sim_pis.append(pi_last)
        else:
            sim_pis.append(yrs_pi['pis'][iyr])

    new_yrs_pi = {'yrs': sim_yrs, 'pis': sim_pis}

    return new_yrs_pi

def fetch_nc_litter(form, fname):
    """
    currently permit only a single cell
    """
    if not exists(fname):
        mess = WARN_STR + 'ORCHIDEE NetCDF litter file name ' + fname
        if fname.isspace() or fname == '':
            mess += ' must not be blank'
        else:
            mess += '  does not exist'

        print(mess)
        return None

    if not check_ochidee_dset(fname):
        return None

    pfts = form.pfts

    lttr_defn = OchideeSet(fname, pfts)
    form.litter_defn = lttr_defn
    form.w_nc_extnt.setText(lttr_defn.nc_extnt)

    # reset PFT combo widget
    # ======================
    form.w_combo_pfts.clear()
    for pft in lttr_defn.aves:
        form.w_combo_pfts.addItem(form.pfts[pft])

    # report average value
    # ====================
    pft_name = form.w_combo_pfts.currentText()
    pft_key = list({elem for elem in pfts if pfts[elem] == pft_name})[0]
    if pft_key in lttr_defn.aves:
        ave_val = lttr_defn.aves[pft_key]
    else:
        ave_val = 0.0

    mess = 'veget type: ' + pft_key + '  average value: ' + str(round(float(ave_val), 2))
    form.w_ave_val.setText(mess)

    return None

def orchidee_pfts():
    """
    Vegetation and corresponding plant functional type as defined in OCHIDEE model
    """
    pfts = {'01': 'SoilBareGlobal',
            '02': 'BroadLeavedEvergreenTropical',
            '03': 'BroadLeavedRaingreenTropical',
            '04': 'NeedleleafEvergreenTemperate',
            '05': 'BroadLeavedEvergreenTemperate',
            '06': 'BroadLeavedSummergreenTemperate',
            '07': 'NeedleleafEvergreenBoreal',
            '08': 'BroadLeavedSummergreenBoreal',
            '09': 'LarixSpBoreal',
            '10': 'C3GrassTemperate',
            '11': 'C4GrassTemperate',
            '12': 'C3AgricultureTemperate',
            '13': 'C4AgricultureTemperate',
            '14': 'C3GrassTropical',
            '15': 'C3GrassBoreal'}

    return pfts
