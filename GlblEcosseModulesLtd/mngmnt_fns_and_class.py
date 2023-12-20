#-------------------------------------------------------------------------------
# Name:        mngmnt_fns_and_class.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------

__prog__ = 'mngmnt_fns_and_class.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import normpath, split, splitext, join, exists, isfile
from pandas import read_excel
from netCDF4 import Dataset, num2date
from numpy import arange, seterr, ma, zeros
from glob import glob

ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '
CNVRSN_FACT = 10000 * 365 / 1000  # gC/m**2/day to kgC/ha/yr
OCHIDEE_MANDAT_VARS = ('lat', 'lon', 'veget', 'TOTAL_BM_LITTER_c')

def check_ochidee_dset(nc_fname):
    """
    C
    """
    nc_fname = normpath(nc_fname)

    nc_dset = Dataset(nc_fname)
    vars_not_prsnt = []
    for var in OCHIDEE_MANDAT_VARS:
        if var not in nc_dset.variables:
            vars_not_prsnt.append(var)
    nc_dset.close()

    if len(vars_not_prsnt) > 0:
        print(ERROR_STR + 'invalid OCHIDEE dataset - variables ' + str(vars_not_prsnt) + ' must be present')
        retcode = False
    else:
        retcode = True

    return retcode

class OchideeSet(object, ):
    """
    Create object from an OCHIDEE NC file
    """
    def __init__(self, nc_fname, pfts):
        """
        assumption is that dataset has been pre-checked using check_ochidee_dset function
        """
        lat_var = 'lat'
        lon_var = 'lon'

        nc_fname = normpath(nc_fname)

        nc_dset = Dataset(nc_fname)
        print('\nReading OCHIDEE file ' + nc_fname)

        lats = nc_dset.variables[lat_var][:]
        nlats = len(lats)
        lons = nc_dset.variables[lon_var][:]
        nlons = len(lons)
        lookup_table = zeros((nlats, nlons), dtype=bool)

        nvegets = len(nc_dset.variables['veget'])

        # record var names
        # ================
        exclude_vars = list([lat_var, lon_var, 'time'])
        start_year = None
        end_year = None
        var_names = []
        for var in nc_dset.variables:
            if var not in exclude_vars:
                var_names.append(var)

            # stanza to get start of time series
            # ==================================
            if var == 'time_centered':

                time_var = nc_dset.variables[var]
                nyears = len(time_var)
                start_date = num2date(int(time_var[0]), units=time_var.units, calendar=time_var.calendar)
                start_year = start_date.year
                end_year = start_year + nyears - 1

        lat_frst = float(lats[0])
        lon_frst = float(lons[0])
        lat_last = float(lats[-1])
        lon_last = float(lons[-1])

        # required for bounding box
        # =========================
        if lat_last > lat_frst:
            lat_ll = lat_frst
            lat_ur = lat_last
        else:
            lat_ll = lat_last
            lat_ur = lat_frst

        if lon_last > lon_frst:
            lon_ll = lon_frst
            lon_ur = lon_last
        else:
            lon_ll = lon_last
            lon_ur = lon_frst

        self.lat_frst = float(lats[0])
        self.lon_frst = float(lons[0])
        self.lat_last = float(lats[-1])
        self.lon_last = float(lons[-1])

        self.lat_var = lat_var
        self.lon_var = lon_var
        self.bbox = lon_ll, lat_ll, lon_ur, lat_ur

        self.nc_fname = nc_fname
        self.var_names = var_names
        self.nc_dset = None

        # resolutions
        # ===========
        self.resol_lon = (lons[-1] - lons[0])/(nlons - 1)
        self.resol_lat = (lats[-1] - lats[0])/(nlats - 1)
        self.max_lat_indx = nlats - 1
        self.max_lon_indx = nlons - 1

        #
        self.lats = list(lats)
        self.lons = list(lons)

        extent_lats = 'N latitudes: {}   extent: {} {}\t'.format(nlats, lat_frst, lat_last)
        extent_lons = 'N longitudes: {}   extent: {} {}\t'.format(nlons, lon_frst, lon_last)
        grid_resol =  'grid resolution: {}'.format(self.resol_lat)
        self.nc_extnt = extent_lats + extent_lons + grid_resol

        # Create a boolean table of cells with and without data
        # =====================================================
        vals = {}
        aves = {}
        for pft_indx in range(nvegets):
            pft_key = '{0:0=2d}'.format(pft_indx + 1)
            tmp_vals = nc_dset.variables['TOTAL_BM_LITTER_c'][:, pft_indx, :, :]

            if tmp_vals.sum() == 0.0:
                print('\t' + WARN_STR + 'No data for PFT ' + pfts[pft_key])
                continue

            vals[pft_key] = tmp_vals
            aves[pft_key] = CNVRSN_FACT * tmp_vals.mean()

            for lat_indx in range(nlats):
                for lon_indx in range(nlons):
                    n_masked = ma.count_masked(vals[pft_key][:, lat_indx, lon_indx])
                    if n_masked == nyears:
                        lookup_table[lat_indx][lon_indx] = False
                    else:
                        lookup_table[lat_indx][lon_indx] = True

        nc_dset.close()

        self.lookup_table = lookup_table
        self.vals = vals
        self.aves = aves
        self.nyears = nyears
        self.start_year = start_year
        self.end_year = end_year

    def get_ochidee_nc_data(self, pft_key, lat, long):
        """
        retrieve data on condition that lat, long are within bounds and data is present
        """
        wthn_bnds = True
        data = None

        lat_indx = int(round((lat - self.lat_frst)/self.resol_lat))
        lon_indx = int(round((long - self.lon_frst)/self.resol_lon))

        if lat_indx < 0 or lat_indx > self.max_lat_indx:
            wthn_bnds = False
            print(WARN_STR + 'latitude index {} out of bounds for latitude {}\tmax indx: {}'.format(lat_indx,
                                                                            round(lat, 4), self.max_lat_indx))

        if lon_indx < 0 or lon_indx > self.max_lon_indx:
            wthn_bnds = False
            print(WARN_STR + 'longitude index {} out of bounds for longitude {}\tmax indx: {}'.format(lon_indx,
                                                                            round(long, 4), self.max_lon_indx))

        if wthn_bnds and self.lookup_table[lat_indx][lon_indx]:
            plnt_inpts = {'yrs': [], 'pis': []}
            strt_yr = self.start_year
            nyears = self.nyears
            vals = list(self.vals[pft_key][:, lat_indx, lon_indx])

            plnt_inpts['yrs'] = [yr for yr in range(strt_yr, strt_yr + nyears)]
            plnt_inpts['pis'] = [CNVRSN_FACT * val for val in vals]
        else:
            plnt_inpts = None

        return plnt_inpts

def check_xls_coords_fname(fname, w_lbl31, data_flag=False):
    """
    C
    """
    results = None
    if not exists(fname):
        print('CSV file ' + fname + ' does not exist')
        return None

    nrecs = 0
    if isfile(fname):
        results = read_excel(fname)
        nrecs = len(results)

    w_lbl31.setText('N recs: {}'.format(nrecs))

    if data_flag:
        return results
    else:
        return nrecs

def get_hilda_land_uses(w_hilda_lus):
    """
    C
    """
    land_uses = []
    for lu in w_hilda_lus:
        if w_hilda_lus[lu].isChecked():
            land_uses.append(lu)

    return land_uses

def check_mask_location(mask_defn, site_rec, land_uses, resol_deg):
    """
    resol_deg is size of cell in degrees
    """
    gran_lat, gran_lon, lat, lon, dummy, dummy = site_rec

    res_d2 = resol_deg/2

    lat_ur = lat + res_d2
    lon_ur = lon + res_d2
    lat_ur_indx, lon_ur_indx, ret_code = mask_defn.get_nc_coords(lat_ur, lon_ur)

    lat_ll = lat - res_d2
    lon_ll = lon - res_d2
    lat_ll_indx, lon_ll_indx, ret_code = mask_defn.get_nc_coords(lat_ll, lon_ll)

    for land_use in land_uses:
        vals = mask_defn.nc_dset.variables[land_use][lat_ll_indx:lat_ur_indx + 1, lon_ll_indx:lon_ur_indx + 1]
        val_mean = vals.mean()
        if val_mean > 0:
            return True

    return False

def create_proj_data_defns(project_path, crop_name, req_resol_deg):
    """
    C
    """
    crop = crop_name.lower()
    resol = str(req_resol_deg)

    resource = 'cropmasks'
    crop_mask_path = join(project_path, resource, crop)
    mask_fnames = glob(crop_mask_path + '\\*' + resol + '*.nc')
    if len(mask_fnames) == 0:
        print('*** Error *** ' + crop_mask_path + ' no mask files found')
        return None

    mask_defn = ManagementSet(mask_fnames[0], resource)

    # for now we just use the mean of the yields from 2000 to 2014;  NB yield is a Python reserved word
    # =================================================================================================
    resource = 'yields'
    obs_path = join(project_path, 'obs_' + crop, resol + 'deg', 'nc')
    yield_fnames = glob(obs_path + '\\*_20*.nc')
    if len(yield_fnames) == 0:
        print('*** Error *** ' + obs_path + ' no yield files found')
        return None

    yield_defn = ManagementSet(yield_fnames[0], resource)

    # fertiliser
    # ==========
    resource = 'fertilizer'
    fert_path = join(project_path, 'fertilizer')
    fert_fnames = glob(fert_path + '\\*.nc4')
    if len(fert_fnames) == 0:
        print('*** Error *** ' + fert_path + ' no ' + resource + ' files found')
        return None

    fert_defns = {}
    for fert_fname in fert_fnames:
        dummy, short_fn = split(fert_fname)
        root_name, dummy = splitext(short_fn)
        # skip not needed dataset
        if root_name == 'Ninput_date_random_ver1':
            continue
        fert_defns[root_name] = ManagementSet(fert_fname, resource)

    fert_defns = fert_defns

    # sowing_harvest
    # ==============
    resource = 'sowing_harvest'
    dates_path = join(project_path, 'sowing_harvest')
    dates_fname = glob(dates_path + '\\' + crop_name + '*.nc')
    if len(dates_fname) == 0:
        print('*** Error *** ' + dates_path + ' no sowing harvest dates file found')
        return None

    dates_defn = ManagementSet(dates_fname[0], resource)

    return mask_defn, yield_defn, dates_defn, fert_defns

def open_proj_NC_sets(mask_defn, yield_defn, dates_defn, fert_defns):
    """
    C
    """
    mask_defn.nc_dset = Dataset(mask_defn.nc_fname)   # mode='r' is default
    yield_defn.nc_dset = Dataset(yield_defn.nc_fname)
    dates_defn.nc_dset = Dataset(dates_defn.nc_fname)
    for metric in fert_defns:
        fert_defns[metric].nc_dset = Dataset(fert_defns[metric].nc_fname)

    return

def close_proj_NC_sets(mask_defn, yield_defn, dates_defn, fert_defns):
    """
    C
    """
    mask_defn.nc_dset.close()
    yield_defn.nc_dset.close()
    dates_defn.nc_dset.close()
    for metric in fert_defns:
        fert_defns[metric].nc_dset.close()

class ManagementSet(object, ):
    """
    C
    """
    def __init__(self, nc_fname, resource):
        """
        C
        """
        nc_fname = normpath(nc_fname)

        nc_dset = Dataset(nc_fname)
        if 'lat' in nc_dset.variables:
            lat_var = 'lat'
            lon_var = 'lon'
        elif 'LAT' in nc_dset.variables:
            lat_var = 'LAT'
            lon_var = 'LON'
        else:
            lat_var = 'latitude'
            lon_var = 'longitude'

        lats = nc_dset.variables[lat_var][:]
        nlats = len(lats)
        lons = nc_dset.variables[lon_var][:]
        nlons = len(lons)

        # record var names
        # ================
        exclude_vars = list([lat_var, lon_var, 'time'])
        start_year = None
        end_year = None
        var_names = []
        for var in nc_dset.variables:
            if var not in exclude_vars:
                var_names.append(var)
            '''
            if var == 'time':
                time_var = nc_dset.variables[var]
                time_units = time_var.units
                since_time = time_units.split('since')[1]
                start_year = int(since_time.split('-')[0]) + time_var[0]    # messy way to get to 1961
                end_year =   start_year + len(time_var) - 1
            '''
        nc_dset.close()

        lat_frst = float(lats[0])
        lon_frst = float(lons[0])
        lat_last = float(lats[-1])
        lon_last = float(lons[-1])

        # required for bounding box
        # =========================
        if lat_last > lat_frst:
            lat_ll = lat_frst
            lat_ur = lat_last
        else:
            lat_ll = lat_last
            lat_ur = lat_frst

        if lon_last > lon_frst:
            lon_ll = lon_frst
            lon_ur = lon_last
        else:
            lon_ll = lon_last
            lon_ur = lon_frst

        self.lat_frst = float(lats[0])
        self.lon_frst = float(lons[0])
        self.lat_last = float(lats[-1])
        self.lon_last = float(lons[-1])

        self.lat_var = lat_var
        self.lon_var = lon_var
        self.bbox = lon_ll, lat_ll, lon_ur, lat_ur

        self.nc_fname = nc_fname
        self.var_names = var_names
        if resource == 'yields':
            self.var_name = None
            for var_name in list(['yield', 'yield_national']):
                if var_name in var_names:
                    self.var_name = var_name
        self.nc_dset = None

        # resolutions
        # ===========
        self.resol_lon = (lons[-1] - lons[0])/(nlons - 1)
        self.resol_lat = (lats[-1] - lats[0])/(nlats - 1)
        self.max_lat_indx = nlats - 1
        self.max_lon_indx = nlons - 1

        #
        self.lats = list(lats)
        self.lons = list(lons)
        self.start_year = start_year
        self.end_year = end_year

    def get_nc_coords(self, latitude, longitude):
        """
        C
        """
        ret_code = 'OK'

        lat_indx = int(round((latitude - self.lat_frst)/self.resol_lat))
        lon_indx = int(round((longitude - self.lon_frst)/self.resol_lon))

        if lat_indx < 0 or lat_indx > self.max_lat_indx:
            ret_code = WARN_STR + 'latitude index {} out of bounds for latitude {}\tmax indx: {}'.format(lat_indx,
                            round(latitude, 4), self.max_lat_indx)
        if lon_indx < 0 or lon_indx > self.max_lon_indx:
            ret_code = WARN_STR + 'longitude index {} out of bounds for longitude {}\tmax indx: {}'.format(lon_indx,
                                                                                round(longitude, 4), self.max_lon_indx)
        return lat_indx, lon_indx, ret_code
