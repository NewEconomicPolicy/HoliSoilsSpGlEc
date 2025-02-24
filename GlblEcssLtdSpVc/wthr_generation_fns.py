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
__prog__ = 'wthr_generation_fns'
__version__ = '0.0.1'
__author__ = 's03mm5'

from time import time
from os.path import join, normpath, isdir, split
from os import listdir, walk, makedirs
from PyQt5.QtWidgets import QApplication

from getClimGenNC_ltd import ClimGenNC
from getClimGenFns_ss import (genLocalGrid, open_wthr_NC_sets, fetch_wthr_dset_overlap, join_hist_fut_to_sim_wthr)
from glbl_ecsse_low_level_fns_sv import update_wthr_progress, update_avemet_progress
from prepare_ecosse_low_level_ss import fetch_long_term_ave_wthr_recs, make_met_files
# from make_site_spec_files_classes import MakeSiteFiles
from hwsd_soil_class import _gran_coords_from_lat_lon

from thornthwaite import thornthwaite

ERROR_STR = '*** Error *** '
WARNING_STR = '*** Warning *** '
QUICK_FLAG = False       # forces break from loops after max cells reached in first GCM and SSP

GRANULARITY = 120
NEXPCTD_MET_FILES = 302
MAX_BANDS = 500000
LTA_RECS_FN = 'lta_ave.txt'

SPACER_LEN = 12

def generate_all_weather(form):
    """
    C
    """
    lat_ll_aoi = form.hwsd_mu_globals.lat_ll_aoi
    lon_ll_aoi = form.hwsd_mu_globals.lon_ll_aoi
    lat_ur_aoi = form.hwsd_mu_globals.lat_ur_aoi
    lon_ur_aoi = form.hwsd_mu_globals.lon_ur_aoi
    bbox_aoi = list([lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi])

    max_cells = int(form.w_max_sims.text())

    resol_deg = 0.5
    resol_d2 = resol_deg/2

    sims_dir = form.sims_dir
    sim_strt_year = 1901

    wthr_set_nm = form.weather_set_linkages['EFISCEN-ISIMIP'][1]
    fut_wthr_set = form.weather_sets[wthr_set_nm]

    wthr_set_nm = form.weather_set_linkages['EFISCEN-ISIMIP'][0]
    hist_wthr_set = form.weather_sets[wthr_set_nm]

    # sim_end_year = form.weather_sets[fut_wthr_set]['year_end']

    # generate weather dataset indices which enclose the AOI
    # ======================================================
    climgen = ClimGenNC(form)
    bbox_wthr = fetch_wthr_dset_overlap(hist_wthr_set, fut_wthr_set)
    bbox_wthr = (12.0, 47.0, 13.0, 48.0)
    aoi_indices_fut = genLocalGrid(fut_wthr_set, bbox_wthr, bbox_aoi)
    aoi_indices_hist = genLocalGrid(hist_wthr_set, bbox_wthr, bbox_aoi)

    # for each GCM and SSP dataset group e.g. UKESM1-0-LL 585
    # =======================================================
    print('')
    ntotal_wrttn, no_data, ncmpltd, nalrdys = 4*[0]
    last_time = time()
    for wthr_set in form.weather_set_linkages['EFISCEN-ISIMIP']:
        this_gcm, scnr = wthr_set.split('_')
        if scnr == 'hist':
            print('Getting historic weather data from weather set: ' + hist_wthr_set['ds_precip'])
            QApplication.processEvents()

            pettmp_hist = climgen.fetch_cru_historic_NC_data(aoi_indices_hist, num_band=-999)

        print('\nGetting future data from weather set: ' + this_gcm + '\tScenario: ' + scnr)
        QApplication.processEvents()
        #      =============================

        pettmp_fut = climgen.fetch_isimap_NC_data(aoi_indices_fut)
        # pettmp_fut = -999

        if pettmp_fut is None or pettmp_hist is None:
            pettmp_sim = None
            no_data += 1
        else:
            pettmp_sim = join_hist_fut_to_sim_wthr(climgen, pettmp_hist, pettmp_fut)

        # create weather
        # ==============
        lat, lon = 2*[None]
        # site_obj = MakeSiteFiles(form, climgen)
        site_obj = None
        make_wthr_files(site_obj, lat, lon, climgen, pettmp_hist, pettmp_sim)
        ncmpltd += 1
        ntotal_wrttn += 1

        last_time = update_wthr_progress(last_time, ncmpltd)
        if ncmpltd >= max_cells:
            break

        # finished this latitude band - report progress
        # =============================================
        mess = 'already existing: {}\tskipped: {}'.format(nalrdys, no_data)
        form.lgr.info(mess)
        print(mess)

        if ncmpltd >= max_cells:
            print('\nFinished checking after {} cells completed'.format(ncmpltd))
            break

        # close NC files
        # ==============
        """
        for metric in list(['precip', 'tas']):
            hist_wthr_dsets[metric].close()
            fut_wthr_dsets[metric].close()
        """
        if QUICK_FLAG:
            break

        print('Completed weather set: ' + this_gcm + '\tScenario: ' + scnr + '\n')

    print('Finished weather generation - total number of sets written: {}'.format(ntotal_wrttn))

    return

def make_avemet_file(clim_dir, lta_precip, lta_pet, lta_tmean):
    """
    will be copied
    """
    avemet_dat = join(clim_dir, 'AVEMET.DAT')
    with open(avemet_dat, 'w') as fobj:
        for imnth, (precip, pet, tmean) in enumerate(zip(lta_precip, lta_pet, lta_tmean)):
            fobj.write('{} {} {} {}\n'.format(imnth + 1, precip, pet, tmean))

    return

def make_wthr_files(site, lat, lon, climgen, pettmp_hist, pettmp_sim):
    """
    generate ECOSSE historic and future weather data
    """
    gran_lat, gran_lon = _gran_coords_from_lat_lon(lat, lon)
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    clim_dir = normpath(join(site.wthr_prj_dir, climgen.region_wthr_dir, gran_coord))
    if not isdir(clim_dir):
        makedirs(clim_dir)  # always create even if no weather data

    if pettmp_hist is None:
        return

    # calculate historic average weather
    # ==================================
    hist_lta_precip, hist_lta_tmean, hist_weather_recs = fetch_long_term_ave_wthr_recs(climgen, pettmp_hist)

    # write a single set of met files for all simulations for this grid cell
    # ======================================================================
    met_fnames = make_met_files(clim_dir, lat, climgen, pettmp_sim)  # future weather

    # create additional weather related files from already existing met files
    # =======================================================================
    irc = climgen.create_FutureAverages(clim_dir, lat, site, hist_lta_precip, hist_lta_tmean)
    lta_ave_fn = _make_lta_file(site, clim_dir)

    return

def fetch_hist_lta_from_lat_lon(sims_dir, climgen, lat, lon):
    """
    check existence of weather cell
    """
    read_lta_flag = True
    integrity_flag, hist_lta_recs, met_fnames = _check_wthr_cell_exstnc(sims_dir, climgen, lat, lon, read_lta_flag)

    return integrity_flag, hist_lta_recs, met_fnames

def _check_wthr_cell_exstnc(sims_dir, climgen, lat, lon, read_lta_flag=False):
    """
    check existence and integrity of weather cell
    allowable criteria are 1) a full set of weather files, namely 300 met files e.g. met2014s.txt, lta_ave.txt and AVEMET.DAT
                           2) an empty directory
    """
    integrity_flag = False
    hist_lta_recs = None
    met_fnames = None
    gran_lat, gran_lon = _gran_coords_from_lat_lon(lat, lon)
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    clim_dir = normpath(join(sims_dir, climgen.region_wthr_dir, gran_coord))
    if isdir(clim_dir):
        fns = listdir(clim_dir)
        nfiles = len(fns)
        if nfiles == 0 or nfiles >= 302:
            if nfiles == 0:
                integrity_flag = True
                hist_lta_recs, met_fnames = None, None
            else:
                if 'lta_ave.txt' in fns:
                    if read_lta_flag:
                        lta_ave_fn = join(clim_dir, 'lta_ave.txt')
                        hist_lta_recs = []
                        with open(lta_ave_fn, 'r') as fave:
                            for line in fave:
                                line = line.rstrip()  # strip out all tailing whitespace
                                hist_lta_recs.append(line)

                    integrity_flag = True
                    met_fnames = fns[2:]

    return integrity_flag, hist_lta_recs, met_fnames

def write_avemet_files(form):
    """
    traverse each GCM and SSP dataset group e.g. UKESM1-0-LL 585
    """
    print('')
    max_cells = int(form.w_max_cells.text())
    sims_dir = form.setup['sims_dir']

    nwrote = 0
    for wthr_set in form.weather_set_linkages['EFISCEN-ISIMIP']:
        wthr_rsrce, scnr = wthr_set.split('_')
        if scnr == 'hist':  # mod
            continue

        # for each region
        # ===============
        for irow, region in enumerate(form.regions['Region']):
            lon_ll, lon_ur, lat_ll, lat_ur, wthr_dir_abbrv = form.regions.iloc[irow][1:]

            # main traversal loop
            # ===================
            region_wthr_dir = wthr_dir_abbrv + wthr_rsrce + '_' + scnr
            clim_dir = normpath(join(sims_dir, region_wthr_dir))

            mess = '\nProcessing weather set: ' + wthr_rsrce + '\tScenario: ' + scnr + '\tRegion: ' + region
            mess += '\t\tabbrev: ' + wthr_dir_abbrv + '\n\tclim_dir: ' + clim_dir
            print(mess)

            if not isdir(clim_dir):
                print(clim_dir + ' *** does not exist ***')
                break

            # step through each directory comprising ECOSSE met files
            # =======================================================
            last_time = time()
            nwrote = 0
            for drctry, subdirs, files in walk(clim_dir):
                nfiles = len(files)
                nsubdirs = len(subdirs)
                if nsubdirs > 0:  # first directory has four scenarios e.g.  Y:\GlblEcssOutputsSv\EcosseSims\AfUKESM1-0-LL_126
                    continue

                # there should be 300 met files plus lta_ave.txt and AVEMET.DAT
                # =============================================================
                last_time = update_avemet_progress(last_time, wthr_rsrce, scnr, region, nwrote)
                if nfiles >= NEXPCTD_MET_FILES:
                    continue

                # if lta_ave.txt is not present then something is wrong
                # =====================================================
                if LTA_RECS_FN in files:
                    lta_ave_fn = join(drctry, LTA_RECS_FN)
                    with open(lta_ave_fn, 'r') as flta_ave:

                        lta_recs = flta_ave.readlines()
                        vals = [float(rec.split('#')[0]) for rec in lta_recs]
                        lta_precip, lta_tmean = vals[:12], vals[12:]

                        gran_coord = split(drctry)[1]
                        gran_lat = int(gran_coord.split('_')[0])
                        cell_lat = 90.0 - gran_lat / GRANULARITY
                        lta_pet = thornthwaite(lta_tmean, cell_lat)

                        make_avemet_file(drctry, lta_precip, lta_pet, lta_tmean)
                        nwrote += 1
                else:
                    print(WARNING_STR + LTA_RECS_FN + ' file should be present in ' + drctry)

            if nwrote >= max_cells:
                print('\nFinished checking having written {} AVEMET.DAT files'.format(nwrote))
                break

            print('Completed Region: ' + region)

        print('Completed weather set: ' + wthr_rsrce + '\tScenario: ' + scnr + '\n')

    print('Finished AVEMET creation - checked: {} cells'.format(nwrote))
    return

def _make_lta_file(site, clim_dir):
    """
    write long term average climate section of site.txt file
    """

    lines = []
    lta_precip, lta_tmean = site.lta_precip, site.lta_tmean
    if lta_precip is None or lta_tmean is None:
        return

    for precip, month in zip(lta_precip, site.months):
        lines.append(_make_line('{}'.format(precip), '{} long term average monthly precipitation [mm]'.format(month)))

    for tmean, month in zip(lta_tmean, site.months):
        lines.append(_make_line('{}'.format(tmean), '{} long term average monthly temperature [mm]'.format(month)))

    lta_ave_fn = join(clim_dir, 'lta_ave.txt')
    with open(lta_ave_fn, 'w') as fhand:
        fhand.writelines(lines)

    # will be copied
    # ==============
    make_avemet_file(clim_dir, site.lta_precip, site.lta_pet, site.lta_tmean)

    return lta_ave_fn

def _make_line(data, comment):
    """

    """
    spacer_len = max(SPACER_LEN - len(data), 2)
    spacer = ' ' * spacer_len

    return '{}{}# {}\n'.format(data, spacer, comment)