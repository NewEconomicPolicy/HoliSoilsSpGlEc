#-------------------------------------------------------------------------------
# Name:        getClimGenNC.py
# Purpose:     read netCDF files comprising ClimGen data
# Author:      s03mm5
# Created:     08/12/2015
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'getClimGenNC.py'
__author__ = 's03mm5'

# Version history
# ---------------
#
from os.path import normpath, isfile, join
from calendar import monthrange
from netCDF4 import Dataset
import math
from numpy import arange, seterr, ma
import warnings
import csv

from thornthwaite import thornthwaite

null_value = -9999
set_spacer_len = 12

numSecsDay = 3600*24
ngranularity = 120
weather_resource_permitted = list(['CRU', 'EObs', 'EWEMBI'])

def _consistency_check(pettmp, varnams_mapped):
    '''
    make sure for a each key if one metric is zero length then all other metrics for that key are also blank
    TODO: this function only works for two metrics and is unpythonic!
    '''
    metric_list = list(varnams_mapped.values())
    metric0 = metric_list[0]
    metric1 = metric_list[1]
    for key in pettmp[metric0]:
        len_key0 = len(pettmp[metric0][key])

        if len_key0 == 0:
            pettmp[metric1][key] = []

        len_key1 = len(pettmp[metric1][key])
        if len_key1 == 0:
            pettmp[metric0][key] = []

    return pettmp

def _check_list_for_none(metric_list):
    '''
    if a None is found then return an empty list
    '''
    for indx, val in enumerate(metric_list):
        if val is None:
            return []

    return metric_list

def _input_txt_line_layout(data, comment):

        spacer_len = max(set_spacer_len - len(data), 2)
        spacer = ' ' * spacer_len
        return '{}{}# {}\n'.format(data, spacer, comment)

class ClimGenNC(object,):

    def __init__(self, form):
        """
        # typically form.inpnc_dir = r'E:\mark2mike\climgenNC'  (get climgen future climate netCDF4 data from here)
        #           form.inp_hist_dir = r'E:\mark2mike\fut_data'  (get CRU historic climate netCDF4 data from here)
        """
        func_name =  __prog__ +  ' ClimGenNC __init__'

        # determine user choices
        # ======================
        if hasattr(form, 'combo10w'):
            weather_resource = form.combo10w.currentText()
            fut_clim_scen = form.combo10.currentText()
            hist_start_year = int(form.combo09s.currentText())
            hist_end_year = int(form.combo09e.currentText())
            ave_weather_flag = form.w_ave_weather.isChecked()
            sim_start_year = int(form.combo11s.currentText())   # might need these later
            sim_end_year = int(form.combo11e.currentText())
        else:
            weather_resource = form.weather_resource
            fut_clim_scen = form.scenario
            hist_start_year = form.hist_strt_year
            hist_end_year = form.hist_end_year
            ave_weather_flag = form.ave_weather_flag
            sim_start_year = form.sim_strt_year
            sim_end_year = form.sim_end_year

        self.weather_resource = weather_resource

        # African Monsoon Multidisciplinary Analysis (AMMA) 2050 datasets
        # ===============================================================
        if weather_resource in form.amma_2050_allowed_gcms:
            wthr_set_key = weather_resource + '_' + fut_clim_scen
            if wthr_set_key not in form.weather_sets:
                print('key {} not in weather sets in function {} - cannot continue'.format(wthr_set_key, func_name))
                return
            hist_weather_set = form.weather_sets[weather_resource + '_historical']
            fut_weather_set  = form.weather_sets[weather_resource + '_' + fut_clim_scen]
            lat = 'lat'
            lon = 'lon'
        elif weather_resource == 'HARMONIE':
            hist_weather_set = form.weather_sets['HARMONIE_V2']
            fut_weather_set  = form.weather_sets['HARMONIE_V2']
            lat = 'lat'
            lon = 'lon'
        elif weather_resource == 'EObs':
            hist_weather_set = form.weather_sets['EObs_Mnth']
            fut_weather_set  = form.weather_sets['EObs_Mnth']
            lat = 'latitude'
            lon = 'longitude'
        elif weather_resource == 'CRU':
            wthr_set_key = 'ClimGen_' + fut_clim_scen
            if wthr_set_key not in form.weather_sets:
                print('key {} not in weather sets in function {} - cannot continue'.format(wthr_set_key, func_name))
                return
            hist_weather_set = form.weather_sets['CRU_hist']
            fut_weather_set  = form.weather_sets['ClimGen_' + fut_clim_scen]
            lat = 'latitude'
            lon = 'longitude'
        elif weather_resource == 'NCAR_CCSM4':
            hist_weather_set = form.weather_sets['NCAR_CCSM4']
            fut_weather_set  = form.weather_sets['NCAR_CCSM4']
            lat = 'lat'
            lon = 'lon'
        else:
            print('weather resource ' + weather_resource + ' not recognised in ' + func_name + ' - cannot continue')
            return

        # make sure start and end years are within dataset limits
        # =======================================================
        hist_start_year = max(hist_weather_set['year_start'], hist_start_year)
        hist_end_year   = min(hist_weather_set['year_end'], hist_end_year)

        self.ave_weather_flag = ave_weather_flag
        num_hist_years = hist_end_year - hist_start_year + 1
        self.num_hist_years = num_hist_years
        self.hist_start_year = hist_start_year
        self.hist_end_year   = hist_end_year
        self.months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        self.fut_clim_scen = fut_clim_scen

        # future data
        # ===========
        self.fut_precip_fname = fut_weather_set['ds_precip']
        self.fut_tas_fname = fut_weather_set['ds_tas']
        self.resolution_lon = fut_weather_set['resol_lon']
        self.resolution_lat = fut_weather_set['resol_lat']

        # granularity
        # ===========
        self.lon = lon
        self.lat = lat
        self.lon_min = fut_weather_set['lon_ll']
        self.lat_min = fut_weather_set['lat_ll']
        self.longitudes = fut_weather_set['longitudes']
        self.latitudes =  fut_weather_set['latitudes']
        self.pettmp = {}        # dictionary whose keys will reference the climate grid, pt_grid
        self.lgr = form.lgr

        # past (monthly) data
        # ===================
        self.hist_precip_fname = hist_weather_set['ds_precip']
        self.hist_tas_fname    = hist_weather_set['ds_tas']
        self.latitudes_hist    = hist_weather_set['latitudes']

        # New stanza to facilitate option when user selects "use average weather"
        # =======================================================================
        if num_hist_years == 1:
            self.met_ave_file = 'met' + str(hist_start_year) + 'a.txt'
        else:
            self.met_ave_file = 'met' + str(hist_start_year) + '_' + str(hist_end_year) + 'a.txt'

        # replaces:
        # num_years_str = '{:0=3d}'.format(num_hist_years)
        # self.met_ave_file = 'met' + num_years_str + 'a.txt'

        # only the years for which we have historic data will be taken into account
        self.num_ave_wthr_years = num_hist_years

        self.sim_start_year = sim_start_year
        self.sim_end_year   = sim_end_year
        self.num_fut_years = sim_end_year - sim_start_year + 1
        self.fut_ave_file = 'met{}_to_{}_ave.txt'.format(sim_start_year, sim_end_year)

    def genLocalGrid(self, bbox, hwsd, snglPntFlag, num_band = None):
        """
        # return the weather indices for the area which encloses the supplied bounding box
        # this function does not alter the ClimGenNC (self) object
        """
        func_name =  __prog__ +  ' genLocalGrid'
        junk = seterr(all='ignore') # switch off warning messages

        bbLonMin, bbLatMin, bbLonMax, bbLatMax =  bbox
        if snglPntFlag:
            bbLonMax = bbLonMin
            bbLatMax = bbLatMin

        # determine bounds for climate grid which will enclose the supplied bounding box
        # ==============================================================================
        resol_lat = self.resolution_lat   # negative for future CRU data
        lat_indices = []
        lat_indices_hist = []
        clim_lat_min = self.lat_min
        num_lats = math.ceil( abs((bbLatMax - clim_lat_min)/resol_lat) )
        latMax = round(abs(num_lats*resol_lat) + clim_lat_min, 8)   # rounding introduced for NCAR_CCSM4
        lat_indices.append(self.latitudes.index(latMax))
        lat_indices_hist.append(self.latitudes_hist.index(latMax))

        num_lats = int( abs((bbLatMin - clim_lat_min)/resol_lat) )
        latMin = round(abs(num_lats*resol_lat) + clim_lat_min, 8)   # rounding introduced for NCAR_CCSM4
        lat_indices.append(self.latitudes.index(latMin))
        lat_indices_hist.append(self.latitudes_hist.index(latMin))

        # longitudes
        # ==========
        lon_indices = []
        resol_lon = self.resolution_lon
        clim_lon_min = self.lon_min
        num_lons = math.ceil((bbLonMax - clim_lon_min)/resol_lon)
        lonMax = round(num_lons*resol_lon + clim_lon_min, 8)
        lon_indices.append(self.longitudes.index(lonMax))

        num_lons = int((bbLonMin - clim_lon_min)/resol_lon)
        lonMin = round(num_lons*resol_lon + clim_lon_min, 8)
        lon_indices.append(self.longitudes.index(lonMin))

        # generate ClimGen grid    NB need to add one when slicing!!!
        # =====================    ==================================
        alons = arange(lonMin, lonMax, resol_lon)
        alats = arange(latMin, latMax, resol_lat)
        nlats = len(alats)
        nlons = len(alons)

        granlons = arange(nlons)
        for ic, lon in enumerate(alons):
            granlons[ic] = (180.0 + lon)*hwsd.granularity
        granlons.sort()

        granlats = arange(nlats)
        for ic, lat in enumerate(alats):
            granlats[ic] = (90.0 - lat)*hwsd.granularity
        granlats.sort()

        # must be in correct order
        # ========================
        lat_indices.sort()
        lat_indices_hist.sort()
        lon_indices.sort()

        aoi_indices_fut = lat_indices + lon_indices
        aoi_indices_hist = lat_indices_hist + lon_indices
        return aoi_indices_fut, aoi_indices_hist

    def fetch_ewembi_NC_data(self, aoi_indices, num_band, future_flag = True):
        '''
        get precipitation or temperature data for a given variable and lat/long index for all times
        CRU uses NETCDF4 format
        '''
        func_name = __prog__ +  ' fetch_fut_future_NC_data'
        warnings.simplefilter('default')

        num_key_masked = 0
        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}

        weather_resource =  self.weather_resource
        if future_flag:
            precip_fname = self.fut_precip_fname
            tas_fname    = self.fut_tas_fname
            start_year   = self.sim_start_year
        else:
            precip_fname = self.hist_precip_fname
            tas_fname    = self.hist_tas_fname
            start_year   = self.hist_start_year

        varnams_mapped = {'pr':'precipitation','tas':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            ncfile = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            slice = ncfile.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]

            if ma.isMaskedArray(slice):
                slice_is_masked_flag = True
                self.lgr.info('Future slice is a masked array in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # generate days per month
            # ======================
            if varname == 'pr':
                days_per_month = []
                nmonths = len(ncfile.variables[varname])
                for year in range(start_year, start_year + int(nmonths/12) + 1):
                    for imnth in range(12):
                        dummy, ndays = monthrange(year, imnth + 1)
                        days_per_month.append(ndays)

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*ngranularity)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*ngranularity)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None
                            num_key_masked += 1

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        if varname == 'pr':
                            pettmp[varnam_map][key] = [round(val*ndays*numSecsDay, 1) for val, ndays in
                                                                                zip(slice[:,ilat,ilon], days_per_month)]
                        elif varname == 'tas':

                            # risk of masked constant
                            # =======================
                            try:
                                pettmp[varnam_map][key] = [round(float(val) - 273.15, 1) for val in slice[:,ilat,ilon]]
                            except UserWarning as e:
                                print('{}\n\tLatitude:{}\tLongitude {}'
                                      .format(e, self.latitudes[lat_indx], self.longitudes[lon_indx]))

            # close netCDF file
            ncfile.close()
            print('# masked weather keys: {}'.format(num_key_masked))

        return pettmp

    def fetch_eobs_NC_data(self, aoi_indices, num_band, future_flag = True):
        '''
        get precipitation or temperature data for a given variable and lat/long index for all times
        EObs uses NETCDF format
        '''
        func_name = __prog__ +  ' fetch_eobs_NC_data'
        warnings.simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            temper_fname = self.fut_tas_fname
        else:
            precip_fname = self.hist_precip_fname
            temper_fname = self.hist_tas_fname

        # process future climate
        # ======================
        varnams_mapped = {'rr':'precipitation','tg':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, temper_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            ncfile = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            try:
                slice = ncfile.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as e:
                print(e)

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                self.lgr.info('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                lat = self.latitudes[lat_indx]
                gran_lat = round((90.0 - lat)*ngranularity)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    long = self.longitudes[lon_indx]
                    gran_lon = round((180.0 + long)*ngranularity)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key {} lat {} long {}'.format(key, lat, long))
                            pettmp[varnam_map][key] = []

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        pettmp[varnam_map][key] = _check_list_for_none(slice[:,ilat,ilon].tolist())

            # close netCDF file
            ncfile.close()

        pettmp = _consistency_check(pettmp, varnams_mapped)
        return pettmp

    def fetch_ncar_ccsm4_NC_data(self, aoi_indices, num_band, future_flag = True):
        '''
        get precipitation or temperature data for a given variable and lat/long index for all times
        CORDEX uses NETCDF3_64BIT format
        '''
        func_name = __prog__ +  ' fetch_ncar_ccsm4_NC_data'
        warnings.simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            temper_fname = self.fut_tas_fname
            start_year   = self.sim_start_year
        else:
            precip_fname = self.hist_precip_fname
            temper_fname = self.hist_tas_fname
            start_year   = self.hist_start_year

        # process future climate
        # ======================
        varnams_mapped = {'pr':'precipitation','tas':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, temper_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            ncfile = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            try:
                slice = ncfile.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as e:
                print(e)

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                print('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # generate days per month
            # ======================
            if varname == 'pr':
                days_per_month = []
                nmonths = len(ncfile.variables[varname])
                for year in range(start_year, start_year + int(nmonths/12) + 1):
                    for imnth in range(12):
                        dummy, ndays = monthrange(year, imnth + 1)
                        days_per_month.append(ndays)

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*ngranularity)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*ngranularity)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        if varname == 'pr':
                            pettmp[varnam_map][key] = [round(val*ndays*numSecsDay, 1) for val in slice[:,ilat,ilon]]
                        elif varname == 'tas':
                            pettmp[varnam_map][key] = [round(val - 273.15, 1) for val in slice[:,ilat,ilon]]

            # close netCDF file
            ncfile.close()

        return pettmp

    def fetch_harmonie_NC_data(self, aoi_indices, num_band, future_flag = True):
        '''
        get precipitation or temperature data for a given variable and lat/long index for all times
        '''
        func_name = __prog__ +  ' fetch_harmonie_NC_data'
        warnings.simplefilter('default')

        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}
        if future_flag == True:
            precip_fname = self.fut_precip_fname
            temper_fname = self.fut_tas_fname
            start_year   = self.sim_start_year
        else:
            precip_fname = self.hist_precip_fname
            temper_fname = self.hist_tas_fname
            start_year   = self.hist_start_year

        # process future climate
        # ======================
        varnams_mapped = {'Precipalign':'precipitation','Tairalign':'temperature'}
        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([precip_fname, temper_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            ncfile = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            try:
                slice = ncfile.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
            except RuntimeWarning as e:
                print(e)

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                print('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # generate days per month
            # ======================
            if varname == 'Precipalign':
                days_per_month = []
                nmonths = len(ncfile.variables[varname])
                for year in range(start_year, start_year + int(nmonths/12) + 1):
                    for imnth in range(12):
                        dummy, ndays = monthrange(year, imnth + 1)
                        days_per_month.append(ndays)

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*ngranularity)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*ngranularity)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        if varname == 'Precipalign':
                            pettmp[varnam_map][key] = [round(val, 2) for val in slice[:,ilat,ilon]]
                        elif varname == 'Tairalign':
                            pettmp[varnam_map][key] = [round(val - 273.15, 1) for val in slice[:,ilat,ilon]]

            # close netCDF file
            ncfile.close()

        return pettmp

    def fetch_cru_future_NC_data(self, aoi_indices, num_band, fut_start_indx = 0):
        '''
        get precipitation or temperature data for a given variable and lat/long index for all times
        CRU uses NETCDF4 format
        '''
        func_name = __prog__ +  ' fetch_fut_future_NC_data'
        warnings.simplefilter('default')

        num_key_masked = 0
        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}

        # process future climate
        # ======================
        varnams_mapped = {'precipitation':'precipitation','temperature':'temperature'}

        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([self.fut_precip_fname, self.fut_tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            ncfile = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================
            slice = ncfile.variables[varname][lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1, :]

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                self.lgr.info('Future slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes[lat_indx])*ngranularity)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*ngranularity)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[ilat,ilon,0]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None
                            num_key_masked += 1

                    # add data for this coordinate
                    # ============================
                    if pettmp[varnam_map][key] == null_value:
                        # remove overlap with historic data - for CRU data only
                        record = [round(val, 1) for val in slice[ilat,ilon,:]]
                        pettmp[varnam_map][key] = record[fut_start_indx:]

            # close netCDF file
            ncfile.close()
            if num_key_masked > 0:
                print('# masked weather keys: {}'.format(num_key_masked))

        return pettmp

    def fetch_cru_historic_NC_data(self, aoi_indices, num_band):
        '''
        get precipitation or temperature data for a given variable and lat/long index for all times
        CRU uses NETCDF4 format
        '''
        func_name = __prog__ +  ' fetch_historic_NC_data'
        warnings.simplefilter('default')

        num_key_masked = 0
        lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
        pettmp = {}

        # process historic climate
        # ========================
        varnams_mapped = {'pre':'precipitation','tmp':'temperature'}

        varnams = sorted(varnams_mapped.keys())

        for varname, fname in zip(varnams, list([self.hist_precip_fname, self.hist_tas_fname])):
            varnam_map = varnams_mapped[varname]
            pettmp[varnam_map] = {}
            ncfile = Dataset(fname, mode='r')

            # collect readings for all time values
            # ====================================

            slice = ncfile.variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]

            if ma.is_masked(slice):
                slice_is_masked_flag = True
                self.lgr.info('Historic weather slice is masked in band {}'.format(num_band))
            else:
                slice_is_masked_flag = False

            # reform slice
            # ============
            for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
                gran_lat = round((90.0 - self.latitudes_hist[lat_indx])*ngranularity)

                for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                    gran_lon = round((180.0 + self.longitudes[lon_indx])*ngranularity)
                    key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))

                    # validate values
                    # ===============
                    pettmp[varnam_map][key] = null_value
                    if slice_is_masked_flag :
                        val = slice[0,ilat,ilon]
                        if val is ma.masked:
                            self.lgr.info('val is ma.masked for key ' + key)
                            pettmp[varnam_map][key] = None
                            num_key_masked += 1

                    # add data for this coordinate
                    if pettmp[varnam_map][key] == null_value:
                        pettmp[varnam_map][key] = [round(val, 1) for val in slice[:, ilat,ilon]]

            # close netCDF file
            ncfile.close()
            if num_key_masked > 0:
                print('# masked weather keys: {}'.format(num_key_masked))

        return pettmp

    def create_FutureAverages(self, clim_dir, lat_inp, granLat, long_inp, granLon):
        '''
        use prexisting metyyyys.txt files to generate a text file of average weather which will subsequently
        be included in the input.txt file
        also create a climate file for each of the simulation years based on average weather from the CRU year range
        '''
        func_name =  ' create_FutureAverages'
        full_func_name =  __prog__ +  func_name

        sim_start_year = self.sim_start_year
        sim_end_year = self.sim_end_year
        months = self.months

        # skip if already exists
        ave_met_file = join(normpath(clim_dir), self.fut_ave_file)
        met_ave_file = join(normpath(clim_dir), self.met_ave_file)
        if isfile(ave_met_file) and isfile(met_ave_file):
            return 0

        # read  precipitation and temperature
        fut_precip = {}
        fut_tmean = {}
        for month in months:
            fut_precip[month] = 0.0
            fut_tmean[month] = 0.0

        for year in range(sim_start_year, sim_end_year):
            fname = 'met{0}s.txt'.format(year)
            met_fpath = join(clim_dir, fname)

            if not isfile(met_fpath):
                print('File ' + met_fpath + ' does not exist - will abandon average weather creation')
                return -1

            with open(met_fpath, 'r', newline='') as fpmet:
                lines = fpmet.readlines()

            for line, month in zip(lines, months):
                tlst = line.split('\t')
                fut_precip[month] += float(tlst[1])
                fut_tmean[month]  += float(tlst[3].rstrip('\r\n'))

        # write stanza for input.txt file consisting of long term average climate
        # =======================================================================
        output = []
        num_fut_years = self.num_fut_years
        for month in self.months:
            ave_precip = fut_precip[month]/num_fut_years
            output.append(_input_txt_line_layout('{}'.format(round(ave_precip,1)), \
                                                '{} long term average monthly precipitation [mm]'.format(month)))

        for month in self.months:
            ave_tmean = fut_tmean[month]/num_fut_years
            output.append(_input_txt_line_layout('{}'.format(round(ave_tmean,2)), \
                                                '{} long term average monthly temperature [degC]'.format(month)))

        # write text file of average weather which will subsequently be included in the input.txt file
        try:
            fhand = open(ave_met_file, 'w')
        except IOError:
            raise IOError('Unable to open file 0}'.format(ave_met_file))
        else:
            fhand.writelines(output)
            fhand.close()

        self.lgr.info('Successfully wrote average weather file {} in function {}'.format(ave_met_file, func_name))

        # write long term average climate file
        # ====================================
        ave_precip = [round(fut_precip[month]/num_fut_years, 1) for month in months]
        ave_tmean  = [round(fut_tmean[month]/num_fut_years, 1) for month in months]

        # pet
        if max(ave_tmean) > 0.0:
            pet = thornthwaite(ave_tmean, lat_inp, year)
        else:
            pet = [0.0]*12
            mess = '*** Warning *** all monthly average temperatures are below zero in ' + full_func_name + \
                   ' for lat/lon: {}/{}\tgranular lat/lon: {}/{}'.format(lat_inp, long_inp, granLat, granLon)
            print(mess)

        pot_evapotrans = [round(p, 1) for p in pet]

        # write file
        output = []
        for tstep, mean_temp in enumerate(ave_tmean):
            output.append([tstep+1, ave_precip[tstep], pot_evapotrans[tstep], mean_temp])

        with open(met_ave_file, 'w', newline='') as fpout:
            writer = csv.writer(fpout, delimiter='\t')
            writer.writerows(output)
            fpout.close()

        self.lgr.info('Successfully wrote average weather file {} in function {}'.format(met_ave_file, func_name))

        return 0