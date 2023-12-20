#-------------------------------------------------------------------------------
# Name:        weather_datasets.py
# Purpose:     script to create weather object and other functions
# Author:      Mike Martin
# Created:     31/07/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'weather_datasets.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import normpath, join, isdir, isfile
from netCDF4 import Dataset, num2date
from glob import glob
from copy import copy
from pandas import read_csv
from thornthwaite import thornthwaite

CHESS_LOOKUP = 'GBR_hwsd_lkup_tble.csv'
EXSTNG_WTHR_RSRCS = ['EObs', 'HARMONIE', 'NCAR_CCSM4', 'CRU', 'CHESS']
REALISATIONS = list(['01', '04', '06', '15'])
RCPS = list(['rcp26', 'rcp45', 'rcp60', 'rcp85'])

NGRANULARITY = 120
MONTH_NAMES_SHORT = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
WARN_STR = '*** Warning *** '

sleepTime = 5

def _fetch_chess_wthr_nc_parms(nc_fname, lkup_tble_df, wthr_rsrce, resol_time, scenario):
    """
    create object describing weather dataset characteristics
    """

    # standard names
    # ==============
    time_var_name = 'time'
    lat = 'lat'
    lon = 'lon'

    nc_fname = normpath(nc_fname)
    nc_dset = Dataset(nc_fname, 'r')
    time_var = nc_dset.variables[time_var_name]
    if 'calendar' in time_var.ncattrs():
        calendar_attr = time_var.calendar
    else:
        calendar_attr = 'standard'

    lat_var = nc_dset.variables[lat]
    lon_var = nc_dset.variables[lon]

    lats = lkup_tble_df['lat'].unique()
    lons = lkup_tble_df['lon'].unique()

    lat_ll = lats.min()
    lon_ll = lons.min()
    lat_ur = lats.max()
    lon_ur = lons.max()

    # resolutions
    # ===========
    resol_lon = round((lon_ur - lon_ll)/(len(lons) - 1), 8)
    resol_lat = round((lat_ur - lat_ll)/(len(lats) - 1), 8)
    if abs(resol_lat) != abs(resol_lon):
        print(WARN_STR + 'weather resource {}\tlat/lon resolutions differ: {} {}'
                                                        .format(wthr_rsrce, resol_lat, resol_lon))

    # Get the start and end date of the time series (as datetime objects):
    # ====================================================================
    time_var_units = time_var.units
    start_day = time_var[0]
    try:
        start_date = num2date(start_day, units = time_var_units, calendar = calendar_attr)
    except (TypeError) as e:
        print('Error deriving start and end year for dataset: ' + nc_fname)
        return None

    end_day = int(time_var[-1])
    end_date = num2date(end_day, units = time_var_units, calendar = calendar_attr)
    start_year = start_date.year
    end_year = end_date.year

    nc_dset.close()

    # construct weather resource
    # ==========================
    wthr_rsrc = {'year_start': start_year,  'year_end': end_year,
            'resol_lat': resol_lat, 'lat_frst': lat_ll, 'lat_last': lat_ur, 'lat_ll': lat_ll, 'lat_ur': lat_ur,
            'resol_lon': resol_lon, 'lon_frst': lon_ll, 'lon_last': lon_ur, 'lon_ll': lon_ll, 'lon_ur': lon_ur,
            'longitudes': lons, 'latitudes': lats,
            'resol_time': resol_time,  'scenario': scenario}

    print('{} start and end year: {} {}\tresolution: {} degrees'
            .format(wthr_rsrce, wthr_rsrc['year_start'],  wthr_rsrc['year_end'], abs(wthr_rsrc['resol_lat'])))

    return wthr_rsrc

def _fetch_weather_nc_parms(nc_fname, wthr_rsrce, resol_time, scenario):
    """
    create object describing weather dataset characteristics
    """

    # standard names
    # ==============
    time_var_name = 'time'
    if wthr_rsrce == 'NASA' or wthr_rsrce[0:5] == 'EObs_' or wthr_rsrce[0:8] == 'ClimGen_':
        lat = 'latitude'
        lon = 'longitude'
    else:
        lat = 'lat'
        lon = 'lon'

    nc_fname = normpath(nc_fname)
    nc_dset = Dataset(nc_fname, 'r')
    time_var = nc_dset.variables[time_var_name]
    if 'calendar' in time_var.ncattrs():
        calendar_attr = time_var.calendar
    else:
        calendar_attr = 'standard'

    lat_var = nc_dset.variables[lat]
    lon_var = nc_dset.variables[lon]

    if wthr_rsrce.find('EObs_') == 0:
        lats = [round(float(lat), 2) for lat in list(lat_var)]  # rounding introduced for EObs
        lons = [round(float(lon), 2) for lon in list(lon_var)]
    else:
        lats = [round(float(lat), 8) for lat in list(lat_var)]  # rounding introduced for NCAR_CCSM4
        lons = [round(float(lon), 8) for lon in list(lon_var)]

    lat_frst = lats[0]
    lon_frst = lons[0]
    lat_last = lats[-1]
    lon_last = lons[-1]

    if lat_last > lat_frst:
        lat_ll = lat_frst; lat_ur = lat_last
    else:
        lat_ll = lat_last; lat_ur = lat_frst

    if lon_last > lon_frst:
        lon_ll = lon_frst; lon_ur = lon_last
    else:
        lon_ll = lon_last; lon_ur = lon_frst

    # resolutions
    # ===========
    resol_lon = (lons[-1] - lons[0])/(len(lons) - 1)
    resol_lat = (lats[-1] - lats[0])/(len(lats) - 1)
    if abs(resol_lat) != abs(resol_lon):
        print('Warning - weather resource {} has different lat/lon resolutions: {} {}'
                                                        .format(wthr_rsrce, resol_lat, resol_lon))

    # Get the start and end date of the time series (as datetime objects):
    # ====================================================================
    if wthr_rsrce[0:8] == 'ClimGen_':
        # print(wthr_rsrce + ' future time units attribute: ' + time_var.units)
        start_year = int(time_var.units.split(' ')[-1])
        end_year = start_year + int(len(time_var)/12) - 1
    else:
        time_var_units = time_var.units
        start_day = int(time_var[0])
        try:
            start_date = num2date(start_day, units = time_var_units, calendar = calendar_attr)
        except (TypeError) as err:
            print('Error deriving start and end year for dataset: ' + nc_fname)
            return None

        end_day = int(time_var[-1])
        end_date = num2date(end_day, units = time_var_units, calendar = calendar_attr)
        start_year = start_date.year
        end_year = end_date.year

    nc_dset.close()

    # construct weather resource
    # ==========================
    wthr_rsrc = {'year_start': start_year,  'year_end': end_year,
            'resol_lat': resol_lat, 'lat_frst': lat_frst, 'lat_last': lat_last, 'lat_ll': lat_ll, 'lat_ur': lat_ur,
            'resol_lon': resol_lon, 'lon_frst': lon_frst, 'lon_last': lon_last, 'lon_ll': lon_ll, 'lon_ur': lon_ur,
            'longitudes': lons, 'latitudes': lats,
            'resol_time': resol_time,  'scenario': scenario}

    print('{} start and end year: {} {}\tresolution: {} degrees'
            .format(wthr_rsrce, wthr_rsrc['year_start'],  wthr_rsrc['year_end'], abs(wthr_rsrc['resol_lat'])))

    return wthr_rsrc

def read_weather_dsets_detail(form, rqrd_rsces = EXSTNG_WTHR_RSRCS):
    """
    ascertain the year span for historic datasets
    TODO: replace with approach adopted for Site Specific version of Global Ecosse
    """

    # weather set linkages
    # ====================
    form.amma_2050_allowed_gcms = {}
    wthr_rsrces_generic = list([])
    weather_set_linkages = {}
    weather_sets = {}

    form.weather_resources_generic = wthr_rsrces_generic
    form.weather_set_linkages = weather_set_linkages
    form.weather_sets = weather_sets
    form.wthr_settings_prev = {}

    if form.settings['weather_dir'] is None:
        return

    weather_dir = form.settings['weather_dir']

    # check CHESS historic
    # ====================
    lkup_tble_flag = False
    gbr_lkup_tble_df = None
    gnrc_rsrce = 'CHESS'
    if gnrc_rsrce in rqrd_rsces:

        # lookup table is mandatory
        # =========================
        wthr_rsrce_base = 'CHESS_historic'
        chss_lkup_tble = join(weather_dir, wthr_rsrce_base, 'lookup_table', CHESS_LOOKUP)
        if isfile(chss_lkup_tble):
            print('Found CHESS lookup table: ' + chss_lkup_tble)
            gbr_lkup_tble_df = read_csv(chss_lkup_tble, sep = ',')
            lkup_tble_flag = True
        else:
            print('Could not find CHESS lookup table: ' + chss_lkup_tble)

    if lkup_tble_flag:

        chss_hist_flag = False
        valid_wthr_dset_rsrces = []
        chss_dir = join(weather_dir, wthr_rsrce_base, 'Monthly')
        if isdir(chss_dir):
            chss_fns = sorted(glob(chss_dir + '/*.nc'))
            chss_pr_fns = [fn for fn in chss_fns if fn.rfind('precip') >= 0]
            if len(chss_pr_fns) > 0:
                chss_tas_fns = [fn for fn in chss_fns if fn.rfind('tas_') >= 0]
                if len(chss_tas_fns) > 0:
                    wthr_rsrce = copy(wthr_rsrce_base)
                    weather_sets[wthr_rsrce] = \
                        _fetch_chess_wthr_nc_parms(chss_fns[0], gbr_lkup_tble_df, wthr_rsrce, 'Monthly', 'historic')
                    weather_sets[wthr_rsrce]['base_dir'] = chss_dir
                    weather_sets[wthr_rsrce]['ds_precip'] = chss_pr_fns[0]
                    weather_sets[wthr_rsrce]['ds_tas'] = chss_tas_fns[0]
                    valid_wthr_dset_rsrces.append(wthr_rsrce)
                    chss_hist_flag = True
            else:
                print('No CHESS historic datasets present in ' + chss_dir)

        chss_hist_dir = copy(chss_dir)

        # check CHESS RCPs
        # ================
        chss_rcp_flag = False
        wthr_rsrce_base = 'CHESS_RCPs'
        chss_rcp_dir = join(weather_dir, wthr_rsrce_base)
        for scenario in RCPS:
            for realis in REALISATIONS:
                chss_dir = join(chss_rcp_dir, 'Monthly', scenario + '_bias-corrected', realis)
                wthr_rsrce = 'chess_' + scenario + '_' + realis
                if isdir(chss_dir):
                    chss_pr_fns = glob(chss_dir + '/*_pr_*.nc')
                    if len(chss_pr_fns) > 0:
                        weather_sets[wthr_rsrce] = \
                            _fetch_chess_wthr_nc_parms(chss_pr_fns[0], gbr_lkup_tble_df, wthr_rsrce, 'Monthly', scenario)
                        weather_sets[wthr_rsrce]['base_dir'] = chss_dir
                        weather_sets[wthr_rsrce]['ds_precip'] = chss_pr_fns[0]
                        weather_sets[wthr_rsrce]['ds_tas'] = glob(chss_dir + '/*_tas_*.nc')[0]
                        valid_wthr_dset_rsrces.append(wthr_rsrce)
                        chss_rcp_flag = True
                else:
                    print('No CHESS RCP datasets present in ' + chss_dir)

        if chss_hist_flag and chss_rcp_flag:
            wthr_rsrces_generic.append(gnrc_rsrce)
            weather_set_linkages[gnrc_rsrce] = valid_wthr_dset_rsrces
        else:
            if not chss_hist_flag :
                print('CHESS historic datasets incomplete in ' + chss_hist_dir)
            if not chss_rcp_flag:
                print('CHESS RCP datasets incomplete in ' + chss_rcp_dir)

    # check EObs monthly: rr_ and tg_
    # ===============================
    gnrc_rsrce = 'EObs'
    if gnrc_rsrce in rqrd_rsces:
        eobs_mnthly_dir  = weather_dir + '\\EObs_v23\\Monthly'
        if isdir(eobs_mnthly_dir):
            wthr_rsrce = 'EObs_Mnth'
            eobs_fnames = glob(eobs_mnthly_dir + '/[rr-tg]*Monthly.nc')
            if len(eobs_fnames) > 0:
                wthr_nc_parms =  _fetch_weather_nc_parms(eobs_fnames[0], wthr_rsrce, 'Monthly', 'historic')
                weather_sets[wthr_rsrce] = wthr_nc_parms
                if  wthr_nc_parms is None:
                    print('Problem reading EObs monthly datasets in ' + eobs_mnthly_dir)
                else:
                    weather_sets[wthr_rsrce]['base_dir']   = eobs_mnthly_dir
                    weather_sets[wthr_rsrce]['ds_precip']  = eobs_fnames[0]
                    weather_sets[wthr_rsrce]['ds_tas']     = eobs_fnames[1]
                    wthr_rsrces_generic.append(gnrc_rsrce)
                    weather_set_linkages[gnrc_rsrce] = list([wthr_rsrce, wthr_rsrce])
            else:
                print('No EObs monthly datasets present in ' + eobs_mnthly_dir)

    # check HARMONIE_V2 monthly: Tair and Precip
    # ==========================================
    gnrc_rsrce = 'HARMONIE'
    if gnrc_rsrce in rqrd_rsces:
        harmonie_dir = weather_dir + '\\HARMONIE_V2\\Monthly'
        if isdir(harmonie_dir):
            wthr_rsrce = 'HARMONIE_V2'
            harmonie_fnames = glob(harmonie_dir + '/cruhar*.nc')
            if len(harmonie_fnames) > 0:
                weather_sets[wthr_rsrce] = _fetch_weather_nc_parms(harmonie_fnames[0], wthr_rsrce, 'Monthly', 'historic')
                weather_sets[wthr_rsrce]['base_dir']   = harmonie_dir
                weather_sets[wthr_rsrce]['ds_precip']  = harmonie_fnames[0]
                weather_sets[wthr_rsrce]['ds_tas']     = harmonie_fnames[1]
                wthr_rsrces_generic.append(gnrc_rsrce)
                weather_set_linkages[gnrc_rsrce] = list([wthr_rsrce, wthr_rsrce])
            else:
                print('No HARMONIE datasets present in ' + harmonie_dir)

    # check NASA monthly
    # ==================
    gnrc_rsrce = 'NCAR_CCSM4'
    if gnrc_rsrce in rqrd_rsces:
        ncar_mnthly_dir  = weather_dir + '\\NCAR_CCSM4\\Monthly'
        if isdir(ncar_mnthly_dir):
            wthr_rsrce = 'NCAR_CCSM4'
            ncar_fnames = glob(ncar_mnthly_dir + '\\rcp26\\*_Amon*.nc')
            if len(ncar_fnames) > 0:
                weather_sets[wthr_rsrce] = _fetch_weather_nc_parms(ncar_fnames[0], wthr_rsrce, 'Monthly', 'historic')
                weather_sets[wthr_rsrce]['base_dir']   = ncar_mnthly_dir
                weather_sets[wthr_rsrce]['ds_precip']  = ncar_fnames[0]
                weather_sets[wthr_rsrce]['ds_tas']     = ncar_fnames[1]
                wthr_rsrces_generic.append(gnrc_rsrce)
                weather_set_linkages[gnrc_rsrce] = list([wthr_rsrce, wthr_rsrce])
            else:
                print('No ' + wthr_rsrce + ' monthly datasets present in ' + ncar_mnthly_dir)

    # check CRU historic
    # ==================
    gnrc_rsrce = 'CRU'
    if gnrc_rsrce in rqrd_rsces:
        cru_flag = False
        valid_wthr_dset_rsrces = []
        cru_dir  = weather_dir + '/CRU_Data'
        if isdir(cru_dir):
            wthr_rsrce = 'CRU_hist'
            cru_fnames = sorted(glob(cru_dir + '/cru*dat.nc'))
            if len(cru_fnames) > 0:
                weather_sets[wthr_rsrce] = _fetch_weather_nc_parms(cru_fnames[0], wthr_rsrce, 'Monthly', 'historic')
                weather_sets[wthr_rsrce]['base_dir']   = cru_dir
                weather_sets[wthr_rsrce]['ds_precip']  = cru_fnames[0]
                weather_sets[wthr_rsrce]['ds_tas']     = cru_fnames[1]
                valid_wthr_dset_rsrces.append(wthr_rsrce)
                cru_flag = True
            else:
                print('No CRU historic datasets present in ' + cru_dir)

        # check ClimGen
        # =============
        climgen_flag = False
        for dset_scenario in list(['A1B','A2','B1','B2']):
            climgen_dir = join(weather_dir, 'ClimGen', dset_scenario)
            wthr_rsrce = 'ClimGen_' + dset_scenario
            if isdir(climgen_dir):
                climgen_fnames = glob(climgen_dir + '/*.nc')
                if len(climgen_fnames) > 0:
                    weather_sets[wthr_rsrce] = _fetch_weather_nc_parms(climgen_fnames[0], wthr_rsrce, 'Monthly', dset_scenario)
                    weather_sets[wthr_rsrce]['base_dir']   = climgen_dir
                    weather_sets[wthr_rsrce]['ds_precip']  = climgen_fnames[0]
                    weather_sets[wthr_rsrce]['ds_tas']     = climgen_fnames[1]
                    valid_wthr_dset_rsrces.append(wthr_rsrce)
                    climgen_flag = True
            else:
                print('ClimGen datasets not present in ' + climgen_dir)

        if cru_flag and climgen_flag:
            wthr_rsrces_generic.append(gnrc_rsrce)
            weather_set_linkages[gnrc_rsrce] = valid_wthr_dset_rsrces
        else:
            if not cru_flag:
                print('CRU historic datasets incomplete in ' + cru_dir)
            if not climgen_flag:
                print('ClimGen future datasets incomplete in ' + climgen_dir)

    form.weather_resources_generic = wthr_rsrces_generic
    form.weather_set_linkages = weather_set_linkages
    form.weather_sets = weather_sets
    form.gbr_lkup_tble_df = gbr_lkup_tble_df

    print('')
    return

def report_aoi_size(form, lon_ll, lat_ll, lon_ur, lat_ur):
    """
    write ASCII climate files
    """
    func_name =  __prog__ + ' report_aoi_size'

    # this will be initially only NASA
    # ================================
    resource = form.combo10w.currentText()
    weather_set = form.weather_sets[resource]
    resol_lat = weather_set['resol_lat']
    lat0 = weather_set['lat0']
    resol_lon = weather_set['resol_lon']
    lon0 = weather_set['lon0']

    lat_indx_ll = int(round((lat_ll - lat0)/resol_lat))
    lon_indx_ll = int(round((lon_ll - lon0)/resol_lon))

    lat_indx_ur = int(round((lat_ur - lat0)/resol_lat))
    lon_indx_ur = int(round((lon_ur - lon0)/resol_lon))

    lat_indx_min = min(lat_indx_ll, lat_indx_ur)
    lat_indx_max = max(lat_indx_ll, lat_indx_ur)
    nlats = lat_indx_max - lat_indx_min + 1

    lon_indx_min = min(lon_indx_ll, lon_indx_ur)
    lon_indx_max = max(lon_indx_ll, lon_indx_ur)
    nlons = lon_indx_max - lon_indx_min + 1

    # get slice for each dataset metric
    # =================================
    mess = 'will retrieve weather for {} locations - nlats/nlons: {} x {} '.format(nlats*nlons, nlats, nlons)

    print(mess)

    return

def write_csv_wthr_file(lgr, country, gcm_name, scenario, latitude, longitude,
                                                        start_year, end_year, pettmp_pr, pettmp_tas, out_dir):
    """
    write to file, simulation weather for the given time period
    """
    func_name =  __prog__ + ' _write_csv_wthr_file'

    metric_list = list(['Precipitation', 'Temperature','Potentional Evapotranspiration'])

    # file comprising rain and temperature
    # ====================================
    short_fname =  country + '_' + gcm_name + '_' + scenario + '.txt'
    metrics_fname = join(out_dir, short_fname)
    try:
        fhand_out = open(metrics_fname, 'w')
    except PermissionError as e:
        print(str(e))
        return

    # stanza for Potentional Evapotranspiration [mm/month]
    # ===================================================
    indx1 = 0
    pettmp_pet = []
    for year in range(start_year, end_year + 1):
        indx2 = indx1 + 12

        # temperature only is required
        # ============================
        tmean = pettmp_tas[indx1:indx2]

        # pet
        if max(tmean) > 0.0:
            pet = thornthwaite(tmean, latitude, year)
        else:
            pet = [0.0]*12
            mess = '*** Warning *** monthly temperatures are all below zero for latitude: {}\tclimate directory: {}'\
                                                                                            .format(latitude, out_dir)
            print(mess)

        pettmp_pet += [round(val,2) for val in pet]

    # identify file
    # =============
    header = 'Area of interest: {}\tLatitude: {}\tLongitude : {}\n'.format(country, latitude, longitude)
    fhand_out.write(header)

    # header for each metric
    # =======================
    header_sub = 'Month'
    for year in range(start_year, end_year + 1):
        header_sub += '\t' + str(year)

    # main loop
    # =========
    for metric, pettmp in zip(metric_list, list([pettmp_pr, pettmp_tas, pettmp_pet])):
        fhand_out.write('\n')
        fhand_out.write(metric + '\n')
        fhand_out.write(header_sub + '\n')

        recs = {}
        for mnth_nme in MONTH_NAMES_SHORT:
            recs[mnth_nme] = mnth_nme

        imnth = 0
        year = start_year
        for val in pettmp:
            mnth_nme = MONTH_NAMES_SHORT[imnth]
            recs[mnth_nme] += '\t' + str(val)
            imnth += 1
            if imnth == 12:
                imnth = 0
                year += 1
                if year > end_year:
                    break

        # write records to file
        # =====================
        for mnth_nme in MONTH_NAMES_SHORT:
            fhand_out.write(recs[mnth_nme] + '\n')

    # end of writing; close file and inform user
    # ==========================================
    nyears =  end_year - start_year + 1
    mess = 'Wrote {} years of weather data for area: {}\tto file: {}\n\tpath: {}'.format(country, nyears, short_fname,
                                                                                         out_dir)
    lgr.info(mess)
    fhand_out.close()

    return

def write_csv_wthr_file_v1(lgr, country, gcm_name, scenario, latitude, longitude,
                        start_year, end_year, pettmp_pr, pettmp_tas, out_dir):
    """
    write to file, simulation weather for the given time period
    """
    func_name =  __prog__ + ' _write_csv_wthr_file'

    # file comprising rain and temperature
    # ====================================
    short_fname =  country + '_' + gcm_name + '_' + scenario + '.txt'
    metrics_fname = join(out_dir, short_fname)
    try:
        fhand_out = open(metrics_fname, 'w')
    except PermissionError as e:
        print(str(e))
        return

    header = 'AOI\tgran_lat\tgran_lon\tlatitude\tlongitude\tyear\tmonth\tprecipitation\ttemperature\n'
    fhand_out.write(header)

    gran_lat = int(round((90.0 - latitude)*NGRANULARITY))
    gran_lon = int(round((180.0 + longitude)*NGRANULARITY))
    prefix = '{}\t{}\t{}\t{}\t{}\t'.format(country, gran_lat, gran_lon, latitude, longitude)

    # write the two metrics to file
    # =============================
    iyear = start_year
    imnth = 0
    irecs = 0
    for rain, temperature in zip(pettmp_pr, pettmp_tas):
        mnth_nme = MONTH_NAMES_SHORT[imnth]
        record = prefix + '{}\t{}\t{:.1f}\t{:.1f}\n'.format(iyear, mnth_nme, rain, temperature)
        imnth += 1
        if imnth == 12:
            imnth = 0
            iyear += 1
            if iyear > end_year:
                fhand_out.write(record)
                irecs += 1
                break

        fhand_out.write(record)
        irecs += 1

    # end of writing; close file and inform user
    # ==========================================
    mess = 'Wrote {} lines of weather data for area: {}\tto file: {}\n\tpath: {}'.format(country, irecs, short_fname,
                                                                                         out_dir)
    lgr.info(mess)
    fhand_out.close()

    return

def record_weather_settings(scenario, hist_strt_year, hist_end_year, fut_strt_year, fut_end_year):
    """
    record weather settings
    """
    previous_settings = {'scenario': scenario, 'hist_strt_year': hist_strt_year, 'hist_end_year': hist_end_year,
                                                    'fut_strt_year': fut_strt_year, 'fut_end_year': fut_end_year}
    return previous_settings

def change_weather_resource(form, weather_resource = None):
    """
    during initialisation, weather_resource will be specified
    otherwise get it from GUI
    """
    if weather_resource == '':
        return
    if weather_resource is None:
        weather_resource = form.combo10w.currentText()
        '''
        form.wthr_settings_prev[weather_resource] = record_weather_settings(scenario, hist_strt_year, hist_end_year,
                                                                                        fut_strt_year, fut_end_year)
        '''

    # invoked when setting up the GUI or when there has been a change in weather resource
    # ===================================================================================
    if weather_resource not in form.weather_set_linkages:
        print('weather resource ' + weather_resource + ' not in weather_set_linkages, cannot proceed')
        return

    weather_set_linkage = form.weather_set_linkages[weather_resource]
    wthr_set_hist = weather_set_linkage[0]
    wthr_set_fut  = weather_set_linkage[1]
    start_year = form.weather_sets[wthr_set_hist]['year_start']
    end_year   = form.weather_sets[wthr_set_hist]['year_end']
    hist_syears = list(range(start_year, end_year))
    hist_eyears = list(range(start_year + 1, end_year + 1))

    # simulation years can extend back into historical period
    # =======================================================
    start_year = min(1970, end_year)
    end_year   = form.weather_sets[wthr_set_fut]['year_end']
    fut_syears = list(range(start_year, end_year))
    fut_eyears = list(range(start_year + 1, end_year + 1))

    # get scenarios from the future weather sets for this resource
    # ============================================================
    scenarios = []
    for wthr_set in weather_set_linkage[1:]:
        scenarios.append(form.weather_sets[wthr_set]['scenario'])

    form.combo10.clear()
    for scenario in scenarios:
        form.combo10.addItem(str(scenario))

    form.combo09s.clear()
    for year in hist_syears:
        form.combo09s.addItem(str(year))

    form.combo09e.clear()
    for year in hist_eyears:
        form.combo09e.addItem(str(year))

    form.combo11s.clear()
    for year in fut_syears:
        form.combo11s.addItem(str(year))

    form.combo11e.clear()
    for year in fut_eyears:
        form.combo11e.addItem(str(year))

    if weather_resource in form.wthr_settings_prev:
        wthr_settings_prev = form.wthr_settings_prev[weather_resource]
        form.combo09s.setCurrentText(wthr_settings_prev['hist_strt_year'])
        form.combo09e.setCurrentText(wthr_settings_prev['hist_end_year'])
        form.combo10.setCurrentText(wthr_settings_prev['scenario'])
        form.combo11s.setCurrentText(wthr_settings_prev['fut_strt_year'])
        form.combo11e.setCurrentText(wthr_settings_prev['fut_end_year'])

    return
