"""
#-------------------------------------------------------------------------------
# Name:        hwsd_glblecsse_fns.py
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
# Description:
#   comprises two functions:
#       def _generate_soil_files(form, climgen,  mask_defn, num_band)
#       def generate_banded_sims(form)
#-------------------------------------------------------------------------------
#
"""

__prog__ = 'glbl_ecsse_high_level_sp.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

import time
from operator import itemgetter
from copy import copy
from netCDF4 import Dataset
from PyQt5.QtWidgets import QApplication

import make_ltd_data_files
import hwsd_bil

from hwsd_mu_globals_fns import gen_grid_cells_for_band
from prepare_ecosse_files import update_progress, make_ecosse_file
from glbl_ecss_cmmn_cmpntsGUI import calculate_grid_cell

def _generate_soil_files(form, num_band):
    """
    Main loop for generating soil data outputs
    """
    func_name =  __prog__ + '\t_generate_soil_files'

    study = form.study
    print('Gathering soil and climate data for study {}...\t\tin {}'.format(study,func_name))
    snglPntFlag = False

    # instantiate a soil grid and climate objects
    hwsd = hwsd_bil.HWSD_bil(form.lgr, form.hwsd_dir)

    # add requested grid resolution attributes to the form object
    calculate_grid_cell(form, hwsd.granularity)
    bbox = form.bbox

    # create grid of mu_globals based on bounding box
    # ===============================================
    nvals_read = hwsd.read_bbox_hwsd_mu_globals(bbox, form.hwsd_mu_globals, form.req_resol_upscale)

    # retrieve dictionary consisting of mu_globals (keys) and number of occurences (values)
    # =====================================================================================
    print('\nRetrieving soil data for band ' + str(num_band))
    QApplication.processEvents()

    mu_globals = hwsd.get_mu_globals_dict()
    if mu_globals is None:
        print('No soil records for AOI: {}\n'.format(bbox))
        return

    mess = 'Retrieved {} values  of HWSD grid consisting of {} rows and {} columns: ' \
          '\n\tnumber of unique mu_globals: {}'.format(nvals_read, hwsd.nlats, hwsd.nlons, len(mu_globals))
    form.lgr.info(mess)

    # for each grid point in the band) return a quintuples of integer Lat/Lon 30 seconds coords, Lat/Lons and mu_global
    hwsd.bad_muglobals = form.hwsd_mu_globals.bad_mu_globals
    aoi_res, bbox = gen_grid_cells_for_band(hwsd, form.req_resol_upscale)
    if form.w_use_high_cover.isChecked():
        aoi_res =  _simplify_aoi(aoi_res)

    lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi = bbox
    num_meta_cells = len(aoi_res)
    print('Band aoi LL lon/lat: {} {}\tUR lon/lat: {} {}\t# meta cells: {}'
                            .format(lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi, num_meta_cells))
    QApplication.processEvents()

    if num_meta_cells == 0:
        mess = 'No aoi_res recs therefore unable to create simulation files... \n'
        print(mess); form.lgr.info(mess)
        return

    mess = 'Generated {} Area of Interest grid cells for band {} '.format(num_meta_cells, num_band)
    form.lgr.info(mess); print(mess)

    print('Writing soil data for band {}...'.format(num_band))
    QApplication.processEvents()

    # open land use NC dataset
    # ========================
    last_time = time.time()
    start_time = time.time()
    completed = 0
    skipped = 0
    landuse_yes = 0
    landuse_no = 0
    warn_count = 0
    no_pis = 0

    # each soil can have one or more dominant soils
    # =======================================================================
    for site_indx, site_rec in enumerate(aoi_res):

        gran_lat, gran_lon, lat, long, area, mu_globals_props = site_rec

        # create limited data object
        # ==========================
        # TODO: ltd_data = make_ltd_data_files.MakeLtdDataFiles(form, comments=True)
        # TODO: make_ecosse_file(form, ltd_data, site_rec, study)
        completed += 1

        last_time = update_progress(last_time, start_time, completed, num_meta_cells, skipped, warn_count)
        QApplication.processEvents()

    mess = '\nBand: {}\tLU yes: {}  LU no: {}\t'.format(num_band, landuse_yes, landuse_no)
    mess += 'skipped: {}\tcompleted: {}\tno plant inputs: {}'.format(skipped, completed, no_pis)
    print(mess); QApplication.processEvents()

    print('')   # spacer
    return

def generate_soil_outputs(form):
    '''
    called from GUI
    '''
    if form.hwsd_mu_globals == None:
        print('Undetermined HWSD aoi - please select a valid HSWD csv file')
        return

    if form.w_use_dom_soil.isChecked():
        use_dom_soil_flag = True
    else:
        use_dom_soil_flag = False

    # make sure bounding box is correctly set
    # =======================================
    lon_ll = -10.25
    lat_ll = 35.25
    lon_ur = 34.75
    lat_ur = 69.75
    form.bbox =  list([lon_ll, lat_ll, lon_ur, lat_ur])

    # lat_ll_aoi is the floor i.e. least latitude, of the HWSD aoi which marks the end of the banding loop
    # ====================================================================================================
    lat_ll_aoi = form.hwsd_mu_globals.lat_ll_aoi
    lon_ll_aoi = form.hwsd_mu_globals.lon_ll_aoi
    lat_ur_aoi = form.hwsd_mu_globals.lat_ur_aoi
    lon_ur_aoi = form.hwsd_mu_globals.lon_ur_aoi
    bbox_aoi = list([lon_ll_aoi,lat_ll_aoi,lon_ur_aoi,lat_ur_aoi])

    # check overlap - study too far to west or east or too far south or north of AOI file
    # ===================================================================================
    if (lon_ur < lon_ll_aoi) or (lon_ll > lon_ur_aoi) or  (lat_ur < lat_ll_aoi) or (lat_ll > lat_ur_aoi):
        print('Error: Study bounding box and HWSD CSV file do not overlap - no simulations are possible')
        return

    # ============================ for each PFT end =====================================
    # print('Study bounding box and HWSD CSV file overlap')
    #        ============================================
    start_at_band = form.start_at_band
    print('Starting at band {}'.format(start_at_band))

    # extract required values from the HWSD database and simplify if requested
    # ========================================================================
    hwsd = hwsd_bil.HWSD_bil(form.lgr, form.hwsd_dir)

    # TODO: patch to be sorted
    # ========================
    mu_global_pairs = {}
    for mu_global in form.hwsd_mu_globals.mu_global_list:
        mu_global_pairs[mu_global] = None

    soil_recs = hwsd.get_soil_recs(mu_global_pairs)  # list is already sorted

    # TODO: patch to be sorted
    # ========================
    for mu_global in hwsd.bad_muglobals:
        del(soil_recs[mu_global])

    form.hwsd_mu_globals.soil_recs = simplify_soil_recs(soil_recs, use_dom_soil_flag)
    form.hwsd_mu_globals.bad_mu_globals = [0] +  hwsd.bad_muglobals
    del(hwsd); del(soil_recs)

    # main banding loop
    # =================
    lat_step = 0.5
    nsteps = int((lat_ur-lat_ll)/lat_step) + 1
    for isec in range(nsteps):
        lat_ll_new = lat_ur - lat_step
        num_band = isec + 1
        '''
        if num_band > 2:       # TODO remove when no longer needed
            print('Exiting from processing after {} bands'.format(num_band - 1))
            break
        '''
        # if the latitude floor of the band has not reached the ceiling of the HWSD aoi then skip this band
        if lat_ll_new > form.hwsd_mu_globals.lat_ur_aoi or num_band < start_at_band:
            print('Skipping out of area band {} of {} with latitude extent of min: {}\tmax: {}\n'
              .format(num_band, nsteps, round(lat_ll_new,6), round(lat_ur, 6)))
        else:

            form.bbox = list([lon_ll, lat_ll_new, lon_ur, lat_ur])

            print('\nProcessing band {} of {} with latitude extent of min: {}\tmax: {}'
                  .format(num_band, nsteps, round(lat_ll_new,6), round(lat_ur, 6)))
            QApplication.processEvents()

            _generate_soil_files(form, num_band)  # does actual work

        # check to see if the last band is completed
        if lat_ll_aoi > lat_ll_new or num_band == nsteps:
            print('Finished processing after {} bands of latitude extents'.format(num_band))
            for ichan in range(len(form.fstudy)):
                form.fstudy[ichan].close()
            break

        lat_ur = lat_ll_new

    return

# ===============================================================
#
def simplify_soil_recs(soil_recs, use_dom_soil_flag):
    """
    compress soil records if duplicates are present
    simplify soil records if requested
    each mu_global points to a group of soils
    a soil group can have up to ten soils
    """
    func_name =  __prog__ + ' _simplify_soil_recs'

    num_raw = 0 # total number of sub-soils
    num_compress = 0 # total number of sub-soils after compressions

    new_soil_recs = {}
    for mu_global in soil_recs:

        # no processing necessary
        # =======================
        num_sub_soils = len(soil_recs[mu_global])
        num_raw += num_sub_soils
        if num_sub_soils == 1:
            num_compress += 1
            new_soil_recs[mu_global] = soil_recs[mu_global]
            continue

        # check each soil for duplicates
        # ==============================
        new_soil_group = []
        soil_group = sorted(soil_recs[mu_global])

        # skip empty groups
        # =================
        if len(soil_group) == 0:
            continue

        first_soil = soil_group[0]
        metrics1 = first_soil[:-1]
        share1   = first_soil[-1]
        for soil in soil_group[1:]:
            metrics2 = soil[:-1]
            share2 =   soil[-1]
            if metrics1 == metrics2:
                share1 += share2
            else:
                new_soil_group.append(metrics1 + [share1])
                metrics1 = metrics2
                share1 = share2

        new_soil_group.append(metrics1 + [share1])
        num_sub_soils = len(new_soil_group)
        num_compress += num_sub_soils
        if num_sub_soils == 1:
            new_soil_recs[mu_global] = new_soil_group
            continue

        if use_dom_soil_flag:
            # assign 100% to the first entry of sorted list
            # =============================================
            dom_soil = copy(sorted(new_soil_group, reverse = True, key=itemgetter(-1))[0])
            dom_soil[-1] = 100.0
            new_soil_recs[mu_global] = list([dom_soil])

    mess = 'Leaving {}\trecords in: {} out: {}'.format(func_name, len(soil_recs),len(new_soil_recs))
    print(mess + '\tnum raw sub-soils: {}\tafter compression: {}'.format(num_raw, num_compress))
    return new_soil_recs

def _simplify_aoi(aoi_res):
    """
    simplify AOI records
    """
    aoi_res_new = []
    j = 0
    for site_rec in aoi_res:
        content = site_rec[-1]
        npairs = len(content)
        if npairs == 0:
            print('No soil information for AOI cell {} - will skip'.format(site_rec))
        elif npairs == 1:
            aoi_res_new.append(site_rec)
        else:
            site_rec_list = list(site_rec)  # convert tuple to a list so we can edit last element
            new_content = sorted(content.items(), reverse = True, key = itemgetter(1))  # sort content so we can pick up most dominant mu_global
            total_proportion = sum(content.values())    # add up proportions
            site_rec_list[-1] = {new_content[0][0]: total_proportion}       # create a new single mu global with summed proportions

            aoi_res_new.append(tuple(site_rec_list)) # convert list to tuple

    return aoi_res_new
