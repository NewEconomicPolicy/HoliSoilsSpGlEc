"""
#-------------------------------------------------------------------------------
# Name:
# Purpose:
# Author:      s03mm5
# Created:     08/12/2015
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""
__version__ = '1.0.00'
__prog__ = 'prepare_ecosse_files.py'

import csv
from os.path import join, exists, split, isdir
from os import makedirs, fsync, remove
from time import time
import sys
from calendar import isleap
import json
from glob import glob
from PyQt5.QtWidgets import QApplication

from thornthwaite import thornthwaite

_monthdays = [31,28,31,30,31,30,31,31,30,31,30,31]
_leap_monthdays = [31,29,31,30,31,30,31,31,30,31,30,31]
set_spacer_len = 12

ERROR_STR = '*** Error *** '

def write_study_definition_file(form):
    """
    write current selections to config file
    """
    # necessary if no simulations have been performed
    # ===============================================
    if 'region_study' not in form.setup:
        return

    # prepare the bounding box
    # ========================
    try:
        ll_lon = float(form.w_ll_lon.text())
        ll_lat = float(form.w_ll_lat.text())
        ur_lon = float(form.w_ur_lon.text())
        ur_lat = float(form.w_ur_lat.text())
    except ValueError:
        ll_lon = 0.0
        ll_lat = 0.0
        ur_lon = 0.0
        ur_lat = 0.0
    bbox = list([ll_lon, ll_lat, ur_lon, ur_lat])

    study = form.setup['region_study']

    weather_resource = form.combo10w.currentText()
    if weather_resource == 'CRU':
        fut_clim_scen = form.combo10.currentText()
    else:
        fut_clim_scen = weather_resource

    # construct land_use change - not elegant but adequate
    # =========================
    land_use = ''
    if hasattr(form, 'lu_pi_content'):
        for indx in form.lu_pi_content['LandusePI']:
            lu, pi = form.lu_pi_content['LandusePI'][indx]
            land_use += form.lu_type_abbrevs[lu] + '2'
    else:
        land_use = 'ara'

    land_use = land_use.rstrip('2')

    # TODO: replace "luCsvFname": "" with 'luPiJsonFname': form.fertiliser_fname
    # ==========================================================================
    if hasattr(form, 'combo00b'):
        crop_name = form.combo00b.currentText()
    else:
        crop_name = 'unknown'

    study_defn = {
        'studyDefn': {
            'bbox': bbox,
            "luCsvFname": "",
            'hwsdCsvFname': '',
            'study': study,
            'land_use': land_use,
            'histStrtYr': form.combo09s.currentText(),
            'histEndYr': form.combo09e.currentText(),
            'climScnr': fut_clim_scen,
            'futStrtYr': form.combo11s.currentText(),
            'futEndYr': form.combo11e.currentText(),
            'cropName': crop_name,
            'province': 'xxxx',
            'resolution': form.req_resol_deg,
            'shpe_file': 'xxxx',
            'version': form.version
        }
    }

    # copy to sims area
    if study != '':
        study_defn_file = join(form.setup['sims_dir'], study + '_study_definition.txt')
        with open(study_defn_file, 'w') as fstudy:
            json.dump(study_defn, fstudy, indent=2, sort_keys=True)
            print('\nWrote study definition file ' + study_defn_file)

    return

def input_txt_line_layout(data, comment):

        spacer_len = max(set_spacer_len - len(data), 2)
        spacer = ' ' * spacer_len
        return '{}{}# {}\n'.format(data, spacer, comment)

def write_line_summary(form, coord_frst, coord_last, max_cells_in_line, max_cells_in_cluster):
    """
    write line summary; function not used
    """
    gran_lat_last, latitude_last, gran_lon_last, longitude_last = coord_last
    gran_lat_frst, latitude_frst, gran_lon_frst, longitude_frst = coord_frst
    form.fstudy[1].write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.
                      format(gran_lat_frst, gran_lon_frst, round(latitude_frst,6), round(longitude_frst,6),
                                gran_lat_last, gran_lon_last, round(latitude_last,6), round(longitude_last,6),
                                                    max_cells_in_line, max_cells_in_cluster))
    form.fstudy[1].flush()
    fsync(form.fstudy[1].fileno())
    return

def write_study_manifest_files(form, lon_lats):
    """
    create study manifest file <study>_manifest.csv and summary <study>_summary_manifest.csv
    write first line and leave file object open for subsequent records
    """
    ngranularity = 120
    NoData = -999

    # set up mainfest file for the study
    # ==================================
    form.fstudy = []
    for frag_name in list(['_','_summary_']):
        study_fname = join(form.setup['sims_dir'], form.setup['study'] + frag_name + 'manifest.csv')
        if exists(study_fname):
            try:
                remove(study_fname)
            except PermissionError as err:
                print('*** Warning *** ' + str(err))
                return

        form.fstudy.append(open(study_fname,'w',100))

    # write first two lines so they are distinguishable for sorting purposes
    # ======================================================================
    for lon_lat_pair in lon_lats:
        longitude, latitude = lon_lat_pair
        gran_lon = round((180.0 + longitude)*ngranularity)
        gran_lat = round((90.0 - latitude)*ngranularity)
        form.fstudy[1].write('{}\t{}\t{}\t{}\t'.format(gran_lat, gran_lon, latitude, longitude))
        form.fstudy[0].write('{}\t{}\t{}\t{}\t{}\n'.format(gran_lat, gran_lon, latitude, longitude, NoData))

    form.fstudy[1].write('{}\t{}\n'.format(form.req_resol_deg, form.req_resol_granul))

    return

def write_manifest_file(study, fut_clim_scen, sim_dir, soil_list, mu_global, latitude, longitude, area):
    """
    write json consisting of mu_global and soil shares
    """

    # location etc.
    # =============
    manifest = {
        'location': {
            'longitude' : float(round(longitude,6)),
            'latitude'  : float(round(latitude,6)),
            'area' : round(area,8),
            'area_description' : study,
            'province' : 'province',
            'scenario' : fut_clim_scen
        }
    }
    # soil shares
    # ===========
    smu_global = str(mu_global)
    manifest[smu_global] =  {}
    for soil_num, soil in enumerate(soil_list):
        manifest[smu_global][soil_num + 1] =  soil[-1]

    # deprecated
    # ==========
    manifest['longitudes'] =  {}
    manifest['granular_longs'] = {}

    # construct file name and write
    # =============================
    manif_dir, fname_part2 = split(sim_dir)
    manifest_fname = join(manif_dir, 'manifest_' + fname_part2[:-4] + '.txt')
    with open(manifest_fname, 'w') as fmanif:
        json.dump(manifest, fmanif, indent=2, sort_keys=True)

    return

def write_signature_file(lgr, sim_dir, mu_global, soil, latitude, longitude, province = '', bad_val = 0):
    """
    write json consisting of mu_global and soil details
    """

    config = {
        'location': {
            'province' : province,
            'longitude' : float(round(longitude,6)),
            'latitude'  : float(round(latitude,6)),
            'share' : soil[-1]
        },
        'soil_lyr1': {
            'C_content': soil[0],
            'Bulk_dens': soil[1],
            'pH'       : soil[2],
            '%_clay': soil[3],
            '%_silt': soil[4],
            '%_sand': soil[5]
        }
    }
    # add subsoil layer, if it exists
    if len(soil) >= 12:
       config['soil_lyr2'] =  {
                'C_content': soil[6],
                'Bulk_dens': soil[7],
                'pH'       : soil[8],
                '%_clay': soil[9],
                '%_silt': soil[10],
                '%_sand': soil[11]
            }

    signature_fname = join(sim_dir, str(mu_global) + '.txt')
    with open(signature_fname, 'w') as fsig:
        try:
            json.dump(config, fsig, indent=2, sort_keys=True)
        except TypeError as err:
            print(str(err))

    return

def update_progress(last_time, ncompleted, nskipped, ntotal_grow, ngrowing, nno_grow, hwsd = None):
    """
    Update progress bar
    """
    new_time = time()
    if new_time - last_time > 5:
        # used to be: Estimated remaining
        mess = '\rCompleted: {:=6d} Growing cells: {:=6d} No grow cells: {:=6d}'.format(ncompleted, ngrowing, nno_grow)

        if hwsd is None:
            bad_muglobals = []
        else:
            bad_muglobals = hwsd.bad_muglobals
        mess += ' Skipped: {:=5d} Bad mu globals: {:=5d}'.format(nskipped, len(bad_muglobals))
        mess += ' Remaining: {:=6d}'.format(ntotal_grow - ncompleted)
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time

def _get_thornthwaite(temp_mean, latitude, year):
    """
    feed daily annual temperatures to Thornthwaite equations to estimate Potential Evapotranspiration [mm/month]
    """
    func_name =  __prog__ + ' _get_thornthwaite'

    ntime_steps = len(temp_mean)

    if ntime_steps == 365:
        month_days = _monthdays
    else:
        month_days =  _leap_monthdays

    indx1 = 0
    monthly_t = []
    for ndays in month_days:
        indx2 = indx1 + ndays
        accum_temp = 0.0
        for temp in temp_mean[indx1:indx2]:
            accum_temp +=temp
        monthly_t.append(accum_temp/ndays)
        indx1 = indx2

    pet_mnthly = thornthwaite(monthly_t, latitude, year)

    # now inflate pet to daily
    # ========================
    pet_daily = []
    for pet_month, ndays in zip(pet_mnthly, month_days):
        pet_day = pet_month/ndays
        for ic in range(ndays):
            pet_daily.append(pet_day)

    # print('{} {} {} {}'.format(year, ntime_steps, len(pet_daily),func_name))
    return pet_daily

def make_met_files(clim_dir, latitude, climgen, pettmp_fut_grid_cell):
    """
    feed annual temperatures to Thornthwaite equations to estimate Potential Evapotranspiration [mm/month]
    """
    if not isdir(clim_dir):
        return

    if pettmp_fut_grid_cell is None:
        return

    precip = pettmp_fut_grid_cell['precip']
    temp   = pettmp_fut_grid_cell['tas']
    start_year = climgen.sim_start_year
    end_year = climgen.sim_end_year
    met_fnames = []

    # check if met files already exist
    # ================================

    nyears = end_year - start_year + 1
    met_files = glob(clim_dir + '\\met*s.txt')
    if len(met_files) == nyears:
        for met_file in met_files:
            dummy, short_name = split(met_file)
            met_fnames.append(short_name)
        return met_fnames


    if pettmp_fut_grid_cell is None:
        return

    # create met files
    # ================
    indx1 = 0
    for year in range(start_year, end_year + 1):
        fname = 'met{}s.txt'.format(year)
        met_fnames.append(fname)
        met_path =join(clim_dir, fname)

        if climgen.sim_mnthly_flag:
            ntime_incrs = 12
        else:
            if isleap(year):
                ntime_incrs = 366
            else:
                ntime_incrs = 365

        indx2 = indx1 + ntime_incrs

        # precipitation and temperature
        precipitation = precip[indx1:indx2]
        temp_mean     = temp[indx1:indx2]

        # pet
        # ===
        if climgen.sim_mnthly_flag:
            pet = thornthwaite(temp_mean, latitude, year)
        else:
            pet = _get_thornthwaite(temp_mean, latitude, year)

        # TODO: do something about occasional runtime warning...
        pot_evapotrans = [round(p, 2) for p in pet]
        precip_out     = [round(p, 2) for p in precipitation]
        tmean_out      = [round(t, 2) for t in temp_mean]

        # write file
        # ==========
        output = []
        for tstep, mean_temp in enumerate(tmean_out):
            output.append([tstep+1, precip_out[tstep], pot_evapotrans[tstep], mean_temp])

        with open(met_path, 'w', newline='') as fpout:
            writer = csv.writer(fpout, delimiter='\t')
            writer.writerows(output)
            fpout.close()

        indx1 += ntime_incrs

    return met_fnames

def fetch_long_term_ave_wthr_recs(climgen, pettmp_hist):
    """
    generate long term average weather records
    """
    func_name = __prog__ + ' fetch_long_term_ave_wthr_recs'

    # calculate historic average weather
    # ==================================
    dset_start_year = climgen.hist_wthr_set_defn['year_start']

    hist_start_year = climgen.hist_start_year
    indx_start = 12*(hist_start_year - dset_start_year)

    hist_end_year = climgen.hist_end_year
    indx_end   = 12*(hist_end_year - dset_start_year + 1) # end year includes all 12 months - TODO: check

    # use dict-comprehension to initialise precip. and temperature dictionaries
    # =========================================================================
    hist_precip = {mnth: 0.0 for mnth in climgen.months}
    hist_tmean  = {mnth: 0.0 for mnth in climgen.months}

    for indx in range(indx_start, indx_end, 12):

        for imnth, month in enumerate(climgen.months):
            hist_precip[month] += pettmp_hist['precip'][indx + imnth]
            hist_tmean[month]  += pettmp_hist['tas'][indx + imnth]

    # write stanza for input.txt file consisting of long term average climate
    # =======================================================================
    hist_weather_recs = []
    num_hist_years = hist_end_year - hist_start_year + 1
    hist_lta_precip = []
    for month in climgen.months:
        ave_precip = hist_precip[month]/num_hist_years
        hist_weather_recs.append(input_txt_line_layout('{}'.format(round(ave_precip,1)), \
                                            '{} long term average monthly precipitation [mm]'.format(month)))
        hist_lta_precip.append(ave_precip)

    hist_lta_tmean = []
    for month in climgen.months:
        ave_tmean = hist_tmean[month]/num_hist_years
        hist_weather_recs.append(input_txt_line_layout('{}'.format(round(ave_tmean,2)), \
                                            '{} long term average monthly temperature [degC]'.format(month)))
        hist_lta_tmean.append(ave_tmean)

    return hist_lta_precip, hist_lta_tmean, hist_weather_recs
