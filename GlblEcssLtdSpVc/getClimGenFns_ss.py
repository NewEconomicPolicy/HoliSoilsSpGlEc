"""
#-------------------------------------------------------------------------------
# Name:        getClimGenFns.py
# Purpose:     additional functions for getClimGenNC.py
# Author:      s03mm5
# Created:     08/02/2018
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------
"""
__prog__ = 'getClimGenFns.py'
__author__ = 's03mm5'

from numpy.ma.core import MaskedConstant, MaskError
from numpy import arange, seterr, ma
from netCDF4 import Dataset
from warnings import filterwarnings
import math

GRANULARITY = 120

ERROR_STR = '*** Error *** '
WARNING = '*** Warning *** '

def genLocalGrid(wthr_set, bbox):
    """
    return the weather indices for the area which encloses the supplied bounding box
    this function does not alter the ClimGenNC (self) object
    """
    # junk = seterr(all='ignore') # switch off warning messages

    lats = wthr_set['latitudes']
    num_lats = len(lats)
    resol_lat = wthr_set['resol_lat']
    lat_min = round(abs(num_lats * resol_lat) + wthr_set['lat_ll'], 8)

    lons= wthr_set['longitudes']
    num_lons = len(lons)
    resol_lon = wthr_set['resol_lon']
    lon_min = round(num_lons*resol_lon + wthr_set['lon_ll'], 8)

    bbLonMin, bbLatMin, bbLonMax, bbLatMax = bbox

    # determine bounds for climate grid which will enclose the supplied bounding box
    # ==============================================================================
    lat_indices = []
    clim_lat_min = lat_min
    num_lats = math.ceil(abs((bbLatMax - clim_lat_min) / resol_lat))
    latMax = round(abs(num_lats * resol_lat) + clim_lat_min, 8)  # rounding introduced for NCAR_CCSM4
    lat_indices.append(lats.index(latMax))

    num_lats = int(abs((bbLatMin - clim_lat_min) / resol_lat))
    latMin = round(abs(num_lats * resol_lat) + clim_lat_min, 8)  # rounding introduced for NCAR_CCSM4
    lat_indices.append(lats.index(latMin))

    # lons
    # ==========
    lon_indices = []
    resol_lon = resol_lon
    clim_lon_min = lon_min
    num_lons = math.ceil((bbLonMax - clim_lon_min) / resol_lon)
    lonMax = round(num_lons * resol_lon + clim_lon_min, 8)
    lon_indices.append(lons.index(lonMax))

    num_lons = int((bbLonMin - clim_lon_min) / resol_lon)
    lonMin = round(num_lons * resol_lon + clim_lon_min, 8)
    lon_indices.append(lons.index(lonMin))

    # generate ClimGen grid    NB need to add one when slicing!!!
    # =====================    ==================================
    alons = arange(lonMin, lonMax, resol_lon)
    alats = arange(latMin, latMax, resol_lat)

    granlons = [(180.0 + lon)*GRANULARITY for lon in alons]
    granlons.sort()

    granlats = [(90.0 - lat)*GRANULARITY for lat in alats]
    granlats.sort()

    # must be in correct order
    # ========================
    lat_indices.sort()
    lon_indices.sort()

    aoi_indices = lat_indices + lon_indices

    return aoi_indices

def _apply_start_year_correction(sim_strt_yr, hist_dset_defn, pettmp):
    """
    assume monthly datasets
    check and, if necessary, correct situation where sim_strt_yr is before historic dataset start year
    """
    repeat_period = (hist_dset_defn['year_start'] - sim_strt_yr)*12
    if repeat_period <= 0:
        return pettmp

    new_pettmp = {}
    for metric in pettmp:
        new_pettmp[metric] = pettmp[metric][0:repeat_period] + pettmp[metric]

    return new_pettmp

def _fetch_wthrset_indices(wthr_set_defn, sim_strt_yr, sim_end_yr):
    """
    get indices for simulation years for monthly weather set
    """
    wthr_yr_strt = wthr_set_defn['year_start']
    wthr_yr_end = wthr_set_defn['year_end']

    # simulation end year is before start of this dataset - nothing to do
    # ===================================================================
    if sim_end_yr < wthr_yr_strt:
        return 3*[None]

    # simulation start year is after end of this dataset - nothing to do
    # ===================================================================
    if wthr_yr_strt > sim_end_yr:
        return 3 * [None]

    indx_strt = max(0, (sim_strt_yr - wthr_yr_strt)*12)

    if sim_end_yr >= wthr_yr_end:

        # simulation end year is in future and beyond this dataset end year
        # =================================================================
        indx_end = -1
        next_strt_yr = wthr_yr_end + 1
    else:
        # simulation end year is before this dataset end year
        # ===================================================
        indx_end = (sim_end_yr - wthr_yr_end)*12
        next_strt_yr = -1

    return indx_strt, indx_end, next_strt_yr

def join_hist_fut_to_sim_wthr(climgen, pettmp_hist, pettmp_fut):
    """
    join historic and future weather
    TODO: can be made more efficient by doing this once
    """
    sim_strt_yr = climgen.sim_start_year
    sim_end_yr = climgen.sim_end_year
    indx_hist_strt, indx_hist_end, next_strt_yr = _fetch_wthrset_indices(climgen.hist_wthr_set_defn,
                                                                            sim_strt_yr, sim_end_yr)
    indx_fut_strt, indx_fut_end, dummy = _fetch_wthrset_indices(climgen.fut_wthr_set_defn,
                                                                            next_strt_yr, sim_end_yr)

    pettmp_sim = {}
    for metric in pettmp_hist:
        if indx_hist_end is not None:
            if indx_hist_end == -1:
                hist_seg = pettmp_hist[metric][indx_hist_strt:]
            else:
                hist_seg = pettmp_hist[metric][indx_hist_strt:indx_hist_end]

            pettmp_sim[metric] = hist_seg
            del hist_seg

        if indx_fut_end is not None:
            if indx_fut_end == -1:
                fut_seg = pettmp_fut[metric][indx_fut_strt:]
            else:
                fut_seg = pettmp_fut[metric][indx_fut_strt:indx_fut_end]

            pettmp_sim[metric] += fut_seg
            del fut_seg

    pettmp_sim = _apply_start_year_correction(sim_strt_yr, climgen.hist_wthr_set_defn, pettmp_sim)

    return pettmp_sim

def open_wthr_NC_sets(climgen):
    """
    C
    """
    hist_wthr_dsets = {}
    fut_wthr_dsets = {}

    for metric, ds_fname in zip(list(['precip', 'tas']), list(['ds_precip', 'ds_tas'])):
        hist_wthr_dsets[metric] = Dataset(climgen.hist_wthr_set_defn[ds_fname])
        fut_wthr_dsets[metric] = Dataset(climgen.fut_wthr_set_defn[ds_fname])

    return hist_wthr_dsets, fut_wthr_dsets

def fetch_WrldClim_data(lgr, lat, lon, climgen, nc_dsets, lat_indx, lon_indx, hist_flag=False):
    """
    C
    """
    filterwarnings("error")

    pettmp = {}
    for metric in list(['precip', 'tas']):
        if hist_flag:
            varname = climgen.hist_wthr_set_defn[metric]
            try:
                slice = nc_dsets[metric].variables[varname][:, lat_indx, lon_indx]
            except BaseException as err:
                print(ERROR_STR + str(err))
                return None
        else:
            varname = climgen.fut_wthr_set_defn[metric]
            slice = nc_dsets[metric].variables[varname][:, lat_indx, lon_indx]

        # test to see if cell data is valid, if not then this location is probably sea
        # =============================================================================
        if type(slice[0]) is MaskedConstant:
            pettmp = None
            mess = 'No data at lat: {} {}\tlon: {} {}\thist_flag: {}'.format(lat, lat_indx, lon, lon_indx, hist_flag)
            lgr.info(mess)
            print(mess)
        else:
            pettmp[metric] = [float(val) for val in slice]

    return pettmp

def get_wthr_nc_coords(dset_defn, latitude, longitude):
    """
    C
    """
    lon_frst = dset_defn['lon_frst']
    lat_frst = dset_defn['lat_frst']
    resol_lat = dset_defn['resol_lat']
    resol_lon = dset_defn['resol_lon']
    max_lat_indx = len(dset_defn['latitudes']) - 1
    max_lon_indx = len(dset_defn['longitudes']) - 1

    lat_indx = int(round((latitude - lat_frst)/resol_lat))
    lon_indx = int(round((longitude - lon_frst)/resol_lon))

    if lat_indx < 0 or lat_indx > max_lat_indx:
        print(WARNING + 'latitude index {} out of bounds for latitude {}\tmax indx: {}'
                                                            .format(lat_indx, round(latitude, 4), max_lat_indx))
        return -1, -1

    if lon_indx < 0 or lon_indx > max_lon_indx:
        print(WARNING + 'longitude index {} out of bounds for longitude {}\tmax indx: {}'
                                                            .format(lon_indx, round(longitude, 4), max_lon_indx))
        return -1, -1

    return lat_indx, lon_indx

def check_clim_nc_limits(form, bbox_aoi = None, wthr_rsrce = 'CRU') :
    """
    this function checks that the specified bounding box lies within extent of the requested weather dataset
    """
    func_name =  __prog__ + ' check_clim_nc_limits'

    limits_ok_flag = True
    '''
    if hasattr(form, 'combo10w'):
        wthr_rsrce = form.combo10w.currentText()
    
    if wthr_rsrce == 'NASA' or wthr_rsrce == 'CRU':
        return limits_ok_flag
    '''
    lon_ll_aoi = float(form.w_ll_lon.text())
    lat_ll_aoi = float(form.w_ll_lat.text())
    lon_ur_aoi = float(form.w_ur_lon.text())
    lat_ur_aoi = float(form.w_ur_lat.text())

    wthr_rsrce = wthr_rsrce + '_hist'      # was + '_Day'
    lat_ur_dset = form.wthr_sets[wthr_rsrce]['lat_ur']
    lon_ur_dset = form.wthr_sets[wthr_rsrce]['lon_ur']
    lat_ll_dset = form.wthr_sets[wthr_rsrce]['lat_ll']
    lon_ll_dset = form.wthr_sets[wthr_rsrce]['lon_ll']

    # similar functionality in lu_extract_fns.py in LU_extract project
    # ================================================================
    if (lon_ll_dset < lon_ll_aoi and lon_ur_dset > lon_ur_aoi) and \
                    (lat_ll_dset < lat_ll_aoi and lat_ur_dset > lat_ur_aoi):
        print('AOI lies within ' + wthr_rsrce + ' weather datasets')
    else:
        print('AOI lies outwith ' + wthr_rsrce + ' weather datasets - LL long/lat: {} {}\tUR long/lat: {} {}'
              .format(lon_ll_dset, lat_ll_dset, lon_ur_dset, lat_ur_dset))
        limits_ok_flag = False

    return limits_ok_flag
