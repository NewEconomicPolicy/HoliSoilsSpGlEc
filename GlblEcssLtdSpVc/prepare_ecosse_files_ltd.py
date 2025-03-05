"""
#-------------------------------------------------------------------------------
# Name:        prepareEcosseFiles.py
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

# Version history
# ---------------
#
from os.path import exists, normpath, join, lexists, basename
from os import remove, makedirs
import csv
import time
import sys
from time import sleep
import shutil
import json
from PyQt5.QtWidgets import QApplication

from thornthwaite import thornthwaite
from glbl_ecss_cmmn_funcs import write_kml_file, write_manifest_file, write_signature_file, input_txt_line_layout
from weather_datasets import write_csv_wthr_file

set_spacer_len = 12
sleepTime = 5

def _weather_for_simulation(amma_2050_allowed_gcms, weather_sets, climgen, pettmp_hist, pettmp_fut):
    """
    return spliced weather for simulation period
    """
    sim_start_year = climgen.sim_start_year
    sim_end_year = climgen.sim_end_year

    # must improve on this TODO
    # =========================
    wthr_rsrc = climgen.weather_resource
    if wthr_rsrc == 'HARMONIE':
        hist_start_year = weather_sets['HARMONIE_V2']['year_start']
        hist_end_year = weather_sets['HARMONIE_V2']['year_end']
        fut_start_year = weather_sets['HARMONIE_V2']['year_start']
    elif wthr_rsrc == 'NCAR_CCSM4':
        hist_start_year = weather_sets['NCAR_CCSM4']['year_start']
        hist_end_year = weather_sets['NCAR_CCSM4']['year_end']
        fut_start_year = weather_sets['NCAR_CCSM4']['year_start']
    elif wthr_rsrc in amma_2050_allowed_gcms:
        hist_start_year = weather_sets[wthr_rsrc + '_historical']['year_start']
        hist_end_year   = weather_sets[wthr_rsrc + '_historical']['year_end']
        fut_start_year  = weather_sets[wthr_rsrc + '_rcp26']['year_start']
    elif wthr_rsrc == 'EObs':
        hist_start_year = weather_sets['EObs_Mnth']['year_start']
        hist_end_year = weather_sets['EObs_Mnth']['year_end']
        fut_start_year = weather_sets['EObs_Mnth']['year_start']
    elif wthr_rsrc == 'EFISCEN-ISIMIP':
        hist_start_year = weather_sets['CRU_hist']['year_start']
        hist_end_year   = weather_sets['CRU_hist']['year_end']
        fut_start_year  = weather_sets[wthr_rsrc + '_ssp126']['year_start']
    else:
        # CRU is default
        hist_start_year = weather_sets['CRU_hist']['year_start']
        hist_end_year = weather_sets['CRU_hist']['year_end']
        fut_start_year = weather_sets['ClimGen_A1B']['year_start']

    pettmp_sim = {}
    if sim_start_year >= fut_start_year:
        indx_strt = 12*(sim_start_year - fut_start_year)
        for metric in pettmp_fut:
            pettmp_sim[metric] = pettmp_fut[metric][indx_strt:]
    else:
        # historic takes priority over future weather
        # ===========================================
        indx_hist_strt = 12*(sim_start_year - hist_start_year)
        indx_fut_strt  = 12*(hist_end_year + 1 - fut_start_year)
        for metric in pettmp_fut:
            pettmp_sim[metric] = pettmp_hist[metric][indx_hist_strt:] + pettmp_fut[metric][indx_fut_strt:]

    return pettmp_sim

def update_progress(last_time, start_time, completed, est_num_sims, skipped, warning_count):

    """Update progress bar."""
    new_time = time.time()
    if new_time - last_time > 5:
        # used to be: Estimated remaining
        mess = '\rCompleted: {:=6d} Skipped: {:=5d} Warnings: {:=5d} Remaining: {:=6d}'\
                                    .format(completed, skipped, warning_count, est_num_sims - completed)
        QApplication.processEvents()
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time

def _make_met_files(clim_dir, latitude, climgen, pettmp_grid_cell):
    '''
    feed annual temperatures to Thornthwaite equations to estimate Potential Evapotranspiration [mm/month]
    '''
    func_name =  '_make_met_files'

    if not lexists(clim_dir):
        try:
            makedirs(clim_dir)
        except FileNotFoundError as e:
            print('Error ' + str(e) + '\n\tin module: ' + __prog__ + ' function: ' + func_name)
            sleep(sleepTime)
            exit(0)

    start_year = climgen.sim_start_year
    end_year   = climgen.sim_end_year
    precip = pettmp_grid_cell['precipitation']         #
    temper = pettmp_grid_cell['temperature']
    nmnths = len(temper)

    indx1 = 0
    for year in range(start_year, end_year + 1):
        fname = 'met{0}s.txt'.format(year)
        met_path = join(clim_dir, fname)

        indx2 = indx1 + 12
        if indx2 > nmnths:
            break

        # precipitation and temperature
        precipitation = precip[indx1:indx2]            #
        tmean =         temper[indx1:indx2]

        # pet
        if max(tmean) > 0.0:
            pet = thornthwaite(tmean, latitude, year)
        else:
            pet = [0.0]*12
            mess = '*** Warning *** monthly temperatures are all below zero for latitude: {}\tclimate directory: {}'\
                                                                                            .format(latitude, clim_dir)
            print(mess)

        # TODO: do something about occasional runtime warning...
        pot_evapotrans = [round(p, 2) for p in pet]
        precip_out     = [round(p, 2) for p in precipitation]
        tmean_out      = [round(t, 2) for t in tmean]

        # write file
        output = []
        for tstep, mean_temp in enumerate(tmean_out):
            output.append([tstep+1, precip_out[tstep], pot_evapotrans[tstep], mean_temp])

        with open(met_path, 'w', newline='') as fpout:
            writer = csv.writer(fpout, delimiter='\t')
            writer.writerows(output)
            fpout.close()

        indx1 += 12

    return

def make_ecosse_file(form, climgen, ltd_data, site_rec, province, pettmp_grid_cell):
    """
    generate sets of Ecosse files for each site
    where each site has one or more soils and each soil can have one or more dominant soils
    pettmp_grid_cell is climate data for this soil grid point
    """
    func_name = 'make_ecosse_file'

    gran_lat, gran_lon, latitude, longitude, area, mu_globals_props = site_rec
    sims_dir = form.sims_dir
    fut_clim_scen = climgen.fut_clim_scen

    # unpack weather
    # ==============
    pettmp_hist = {}
    pettmp_fut  = {}
    for metric in list(['precipitation','temperature']):
        pettmp_hist[metric] = pettmp_grid_cell[metric][0]
        pettmp_fut[metric]  = pettmp_grid_cell[metric][1]

    # calculate historic average weather
    # ==================================
    wthr_rsrc = climgen.weather_resource
    if wthr_rsrc == 'NCAR_CCSM4':
        dset_start_year = form.weather_sets['NCAR_CCSM4']['year_start']
    elif wthr_rsrc == 'HARMONIE':
        dset_start_year = form.weather_sets['HARMONIE_V2']['year_start']
    elif wthr_rsrc in form.amma_2050_allowed_gcms:
        dset_start_year = form.weather_sets[wthr_rsrc + '_historical']['year_start']
    elif wthr_rsrc == 'EObs':
        dset_start_year = form.weather_sets['EObs_Mnth']['year_start']
    else:
        dset_start_year = form.weather_sets['CRU_hist']['year_start']

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
            try:
                hist_precip[month] += pettmp_hist['precipitation'][indx + imnth]
                hist_tmean[month]  += pettmp_hist['temperature'][indx + imnth]
            except IndexError as err:
                mess = str(err) + ' in ' + func_name
                print(mess); form.lgr.info(mess)
                return

    # write stanza for input.txt file consisting of long term average climate
    # =======================================================================
    hist_weather_recs = []
    num_hist_years = hist_end_year - hist_start_year + 1
    for month in climgen.months:
        ave_precip = hist_precip[month]/num_hist_years
        hist_weather_recs.append(input_txt_line_layout('{}'.format(round(ave_precip,1)), \
                                            '{} long term average monthly precipitation [mm]'.format(month)))

    for month in climgen.months:
        ave_tmean = hist_tmean[month]/num_hist_years
        hist_weather_recs.append(input_txt_line_layout('{}'.format(round(ave_tmean,2)), \
                                            '{} long term average monthly temperature [degC]'.format(month)))

    # write met files set for all simulations for this grid cell
    # ==========================================================
    if hasattr(form, 'study'):
        study = form.study
    else:
        study = form.w_study.text()
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)

    met_rel_path = join('..', '..', '..', 'Wthr', climgen.fut_clim_scen, gran_coord, '')

    clim_dir = normpath(join(climgen.wthr_out_dir, gran_coord))
    simulation_weather = _weather_for_simulation(form.amma_2050_allowed_gcms, form.weather_sets, climgen,
                                                                                            pettmp_hist, pettmp_fut)
    _make_met_files(clim_dir, latitude, climgen, simulation_weather)

    # create additional weather related files from already existing met files
    irc = climgen.create_FutureAverages(clim_dir, latitude, gran_lat, longitude, gran_lon)

    # TODO: improve
    write_csv_wthr_file(form.lgr, study, wthr_rsrc, fut_clim_scen, latitude, longitude,
                                climgen.sim_start_year, climgen.sim_end_year,
                                pettmp_fut['precipitation'], pettmp_fut['temperature'], clim_dir)

    #------------------------------------------------------------------
    # Create a set of simulation input files for each dominant
    # soil-land use type combination
    #------------------------------------------------------------------
    # construct directory name with all dominant soils

    for pair in mu_globals_props.items():
        mu_global, proportion = pair
        area_for_soil = area*proportion
        soil_list = form.hwsd_mu_globals.soil_recs[mu_global]

        for soil_num, soil in enumerate(soil_list):
            identifer = 'lat{0:0=7d}_lon{1:0=7d}_mu{2:0=5d}_s{3:0=2d}'.format(gran_lat, gran_lon, mu_global, soil_num + 1)

            sim_dir = join(sims_dir, study, identifer)

            if not lexists(sim_dir):
                makedirs(sim_dir)

            ltd_data.write(sim_dir, soil, latitude, hist_weather_recs, met_rel_path)

            # write kml file if requested and signature file
            # ==============================================
            if form.kml_flag and soil_num == 0:
                write_kml_file(sim_dir,  str(mu_global), mu_global, latitude, longitude)

            write_signature_file(sim_dir, mu_global, soil, latitude, longitude, province)

            # copy across Model_Switches.dat file
            # ===================================
            outMdlSwtchs = join(sim_dir, basename(form.default_model_switches))
            shutil.copyfile(form.default_model_switches, outMdlSwtchs)

        # TODO: is manifest file is essential for subsequent processing?
        # ====================================================
        write_kml_file(sim_dir, str(mu_global), mu_global, latitude, longitude,)
        write_manifest_file(study, fut_clim_scen, sim_dir, soil_list, mu_global, latitude, longitude, area_for_soil)

    # end of Soil loop
    # ================

    return
