"""
#-------------------------------------------------------------------------------
# Name:
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""

__prog__ = 'plant_input_fns.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

from os.path import exists
from math import sqrt
from pandas import read_csv, DataFrame
from netCDF4 import Dataset

sleepTime = 5
max_lines = 10

def check_plant_input_nc(form, pi_nc_fname, var_name='PlantInput'):
    """
    check if plant inputs netCDF file is present and valid
    adjust w_combo15 items appropiately
    """
    warn_mess = '*** Warning *** '
    form.w_combo15.clear()
    if not exists(pi_nc_fname):
        form.w_use_pi_nc.setCheckState(0)
        print(warn_mess + ' plant inputs NC file ' + pi_nc_fname + ' does not exist')
        return False

    nc_dset = Dataset(pi_nc_fname)
    vars_found = []
    for var in nc_dset.variables:
        if var.find(var_name) >= 0:
            vars_found.append(var)
    nc_dset.close()

    if len(vars_found) > 0:
        print('Found variables:' + str(vars_found) + ' in ' + pi_nc_fname + '\n')
        # populate drop down
        # =================
        if form.check_plant_input_nc is not None:
            for var in vars_found:
                form.w_combo15.addItem(var)
        return True
    else:
        if hasattr(form, 'w_use_pi_nc'):
            form.w_use_pi_nc.setCheckState(0)
        print('No required variables found in: ' + pi_nc_fname + '\n')
        return False

def associate_yield_nc(logger_info, latitude, longitude, ltd_data, yield_defn, yield_dset, pi_var,
                                                                                            write_to_logger = False):
    """
    modify ltd_data object with plant inputs - assumes NC dataset is already open
    """

    # get plant inputs nearest to this lat/lon
    # ========================================

    lat_indx, lon_indx, ret_code = yield_defn.get_nc_coords(latitude, longitude)
    pi_vals = yield_dset[pi_var][:, lat_indx, lon_indx]

    if write_to_logger:
        mess = 'Plant input: {} at latitude: {}\tlongitude: {}'.format(str(pi_vals),
                                                                       round(latitude, 3), round(longitude, 3))
        logger_info(mess)

    # only overwrite if arable
    # ========================
    for iyr, land_luse in enumerate(ltd_data.landUses):
        if land_luse == 1 or land_luse == 2:
            ltd_data.plantInput[iyr] = round(float(pi_vals[iyr]),3)

    return

def fetch_yields(form):
    """
    read the CSV of yields
    ======================
    """

    yield_map_fname = form.yield_map_fname
    if yield_map_fname is None:
        return

    data_frame = read_csv(yield_map_fname, sep = ',')

    nlines = len(data_frame)
    print('Read {} lines from yield file {}'.format(nlines, yield_map_fname))
    if 'lon' not in data_frame.columns or 'lat' not in data_frame.columns:
        print('Yields file ' + yield_map_fname + ' must have fields lon and lat')
        return DataFrame()

    # create set of points to assist with location
    # ============================================
    data_frame['point'] = [(y, x) for y, x in zip(data_frame['lat'], data_frame['lon'])]

    lats = sorted(data_frame['lat'].unique())
    lons = sorted(data_frame['lon'].unique())
    resols = []
    lat1 = lats[0]
    for lat2 in lats[1:]:
        step = lat2 - lat1
        resols.append(step)
        lat1 = lat2
    resols.sort()
    print('Resolution: {}'.format(resols[0]))

    return data_frame

def associate_yield(logger_info, latitude, longitude, ltd_data, yield_frame):
    """
    modify ltd_data object with plant inputs for arable land-use
    """
    if yield_frame is None:
        return

    # get fertiliser nearest to this lat/lon
    # ======================================
    yield_frame['cdist'] = [( sqrt((latitude - pnt[0])**2 + (longitude - pnt[1])**2) ) for pnt in yield_frame['point']]
    recid_clos = yield_frame['cdist'].idxmin()  # closest
    recid_frth = yield_frame['cdist'].idxmax()  # furthest
    yield_val = yield_frame['yield'][recid_clos]
    mess = 'Value of yield: {} at latitude: {} {}\tlongitude: {} {}\trecord nearest: {}\tfurthest: {}'\
                .format(yield_val, round(latitude,3), yield_frame['lat'][recid_clos], round(longitude,3),
                                                            yield_frame['lon'][recid_clos], recid_clos, recid_frth)
    logger_info(mess)

    # step through land-uses looking for arable
    # =========================================
    for iyr, land_luse in enumerate(ltd_data.landUses):
        if land_luse == 1:
            ltd_data.plantInput[iyr] = yield_val

    return
