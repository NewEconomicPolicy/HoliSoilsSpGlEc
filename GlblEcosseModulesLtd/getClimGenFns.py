#-------------------------------------------------------------------------------
# Name:        getClimGenFns.py
# Purpose:     additional functions for getClimGenNC.py
# Author:      s03mm5
# Created:     08/02/2018
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'getClimGenFns.py'
__author__ = 's03mm5'

import time
import sys
from glob import glob
import netCDF4 as cdf
import numpy as np
import math

set_spacer_len = 12
ngranularity = 120
delta_deg = 0.001 # ensures creation of enclosing climate grid

def _write_coords_for_key(mess, climgen, proximate_keys, lookup_key, func_name):

    # should go in log file - TODO
    gran_lat, gran_lon = proximate_keys[lookup_key]
    latitude  = 90 - gran_lat/ngranularity
    longitude = gran_lon/ngranularity - 180
    mess += ' Lat: {}\tGran lat: {}\tLon: {}\tGran lon: {}\tin {}'.format(latitude, gran_lat, longitude, gran_lon,
                                                                                                        func_name)
    climgen.lgr.info(mess)

    return

def associate_climate(site_rec, climgen, pettmp_hist, pettmp_fut):
    """
    this function associates each soil grid point with the most proximate climate data grid cell
    at the time of writing (Dec 2015) HWSD soil data is on a 30 arc second grid whereas climate data is on 30 or 15 or
     7.5 arc minute grid i.e. 0.5 or 0.25 or 0.125 of a degree
    """
    func_name =  __prog__ + ' associate_climate'

    proximate_keys = {}
    gran_lat_cell, gran_lon_cell, latitude, longitude, dummy, dummy = site_rec
    metric_list = pettmp_fut.keys()

    # TODO: find a more elegant methodology
    # =====================================
    for lookup_key in pettmp_hist['precipitation']:
        if pettmp_fut['precipitation'][lookup_key] == None or pettmp_fut['temperature'][lookup_key] == None or \
                    pettmp_hist['precipitation'][lookup_key] == None or pettmp_fut['temperature'][lookup_key] == None:
            continue
        else:
            slat, slon = lookup_key.split('_')
            gran_lat = int(slat)
            gran_lon = int(slon)

            # situation where grid cell is coincidental with weather cell
            # ===========================================================
            if gran_lat == gran_lat_cell and gran_lon == gran_lon_cell:
                climgen.lgr.info('Cell with lookup key ' + lookup_key + ' is coincidental with weather cell')
                pettmp_out = {}
                for metric in pettmp_fut.keys():
                    pettmp_out[metric] = list([pettmp_hist[metric][lookup_key], pettmp_fut[metric][lookup_key]])
                return pettmp_out
            else:
                proximate_keys[lookup_key] = list([gran_lat, gran_lon])

    # return empty dict if no proximate keys (unlikely)
    # =================================================
    if len(proximate_keys) == 0:
        print('\nNo weather keys assigned for site record with granular coordinates: {} {}\tand lat/lon: {} {}'
                                        .format(gran_lat_cell, gran_lon_cell, round(latitude,4), round(longitude,4)))
        return {}

    '''
    use the minimum distance to assign weather for specified grid cell
    '''

    # first stanza: calculate the squares of the distances in granular units between the grid cell and weathers cells
    # =============
    dist = {}
    total_dist = 0
    for lookup_key in proximate_keys:
        gran_lat, gran_lon = proximate_keys[lookup_key]
        # _write_coords_for_key('\t', proximate_keys, lookup_key, func_name)

        dist[lookup_key] = (gran_lat - gran_lat_cell)**2 + (gran_lon - gran_lon_cell)**2
        total_dist += dist[lookup_key]

    # find key corresponding to the minimum value using conversions to lists
    # ======================================================================
    minval = sorted(dist.values())[0]
    lookup_key = list(dist.keys())[list(dist.values()).index(minval)]
    _write_coords_for_key('Selected weather key', climgen, proximate_keys, lookup_key, func_name)

    #
    pettmp_final = {}
    for metric in metric_list:
        pettmp_final[metric] = list([pettmp_hist[metric][lookup_key], pettmp_fut[metric][lookup_key]])

    return pettmp_final

def check_clim_nc_limits(form, weather_resource, bbox_aoi = None):

    """
    this function makes sure that the specified bounding box lies within extent of the requested weather dataset
    NB lats run from North to South
        lons run from West to East
    """
    func_name =  __prog__ + ' check_clim_nc_limits'

    limits_ok_flag = True

    # CRU covers the globe
    # ====================
    if weather_resource == 'CRU':
        return limits_ok_flag
    #
    wthr_set_name = form.weather_set_linkages[weather_resource][0]

    if bbox_aoi is None:
        lon_ur_aoi = float(form.w_ur_lon.text())
        lat_ur_aoi = float(form.w_ur_lat.text())
        if form.version == 'HWSD_grid':
            lon_ll_aoi = float(form.w_ll_lon.text())
            lat_ll_aoi = float(form.w_ll_lat.text())
        else:
            lon_ll_aoi = lon_ur_aoi
            lat_ll_aoi = lat_ur_aoi
    else:
        lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi = bbox_aoi

    lat_ur_dset = form.weather_sets[wthr_set_name]['lat_ur']
    lon_ur_dset = form.weather_sets[wthr_set_name]['lon_ur']
    lat_ll_dset = form.weather_sets[wthr_set_name]['lat_ll']
    lon_ll_dset = form.weather_sets[wthr_set_name]['lon_ll']

    # similar functionality in lu_extract_fns.py in LU_extract project
    # ================================================================
    if (lon_ll_dset < lon_ll_aoi and lon_ur_dset > lon_ur_aoi) and \
                    (lat_ll_dset < lat_ll_aoi and lat_ur_dset > lat_ur_aoi):
        print('AOI lies within ' + weather_resource + ' weather dataset')
    else:
        print('AOI lies outwith ' + weather_resource + ' weather dataset - LL long/lat: {} {}\tUR long/lat: {} {}'
              .format(lon_ll_dset, lat_ll_dset, lon_ur_dset, lat_ur_dset))
        limits_ok_flag = False

    return limits_ok_flag

def update_progress_clim_soil(last_time, nsoilres, pt_key, ncsv_lines, skipped = 0, failed = 0):

    """Update progress bar."""
    new_time = time.time()
    if new_time - last_time > 5:

        mess = '\rSize of soil list: {}\tpt_key: {}\tNumber of sites remaining: {}'\
                                            .format(nsoilres, pt_key, ncsv_lines - nsoilres)
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time

def update_progress_clim(last_time, ngrid_cells, total_num, no_data):
    """Update progress bar."""
    new_time = time.time()
    if new_time - last_time > 2:
        mess = '\rCompleted checking of: {} climate cells\tNo data: {}\tRemaining: {}'\
                                            .format(ngrid_cells, no_data, 1 + total_num - ngrid_cells)
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time
