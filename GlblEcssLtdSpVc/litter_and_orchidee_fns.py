#-------------------------------------------------------------------------------
# Name:        litter_and_orchidee_fns.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------

__prog__ = 'litter_and_orchidee_fns.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import exists, normpath, split, splitext, join, exists, isfile
from numpy import arange, seterr, ma, zeros
from netCDF4 import Dataset, num2date
from PyQt5.QtWidgets import QApplication

ERROR_STR = '*** Error *** '
WARN_STR = '*** Warning *** '

CNVRSN_FACT_BM_LITTER = 10000 * 365 / 1000  # gC/m**2/day to kgC/ha/yr
CNVRSN_FACT_LITTER_SOIL = 10000 / 1000      # gC/m**2 to kgC/ha

OCHIDEE_MANDAT_VARS = ('lat', 'lon', 'veget', 'TOTAL_BM_LITTER_c', 'TOTAL_LITTER_SOIL_c', 'TOTAL_SOIL_c')

EFISCEN_CARBON_VAR = 'TOTAL_BM_LITTER_c'
EFISCEN_MANDAT_VARS = ('lat', 'lon', EFISCEN_CARBON_VAR)       # EFISCEN | European Forest Institute
OCHIDEE_OVERRIDE_YEAR = 1970

def fetch_nc_litter(form, fname):
    """
    currently permit only a single cell
    """
    if not exists(fname):
        mess = WARN_STR + 'ORCHIDEE or EFISCEN NetCDF litter file name ' + fname
        if fname.isspace() or fname == '':
            mess += ' must not be blank'
        else:
            mess += '  does not exist'

        print(mess)
        return None

    if not check_ochidee_dset(fname):
        return None

    pfts = form.pfts

    carbon_var = form.combo08.currentText()
    lttr_defn = OchideeSet(fname, pfts, carbon_var)
    form.litter_defn = lttr_defn
    form.w_nc_extnt.setText(lttr_defn.nc_extnt)

    # reset PFT combo widget
    # ======================
    form.w_combo_pfts.clear()
    if lttr_defn.dset_type == 'EFISCEN':
        form.combo08.setCurrentText(EFISCEN_CARBON_VAR)
        form.combo08.setEnabled(False)
        form.w_var_desc.setText(form.carbon_vars[carbon_var])
        form.w_combo_pfts.addItem('total biomass litter C')
    else:
        form.combo08.setEnabled(True)
        for pft in lttr_defn.aves:
            form.w_combo_pfts.addItem(form.pfts[pft])

    # report average value
    # ====================
    pft_name = form.w_combo_pfts.currentText()
    if pft_name == '':
        mess = 'No data for ' + carbon_var
        form.w_ave_val.setText(mess)
    else:
        if form.litter_defn.dset_type == 'EFISCEN':
            pft_key = '00'
            ave_val = form.litter_defn.aves[pft_key]
            mess = 'average value: ' + str(round(float(ave_val), 2))
        else:
            pft_key = list({elem for elem in pfts if pfts[elem] == pft_name})[0]
            if pft_key in lttr_defn.aves:
                ave_val = lttr_defn.aves[pft_key]
            else:
                ave_val = 0.0

            mess = 'veget type: ' + pft_key + '  average value: ' + str(round(float(ave_val), 2))

    form.w_ave_val.setText(mess)

    return None

def check_ochidee_dset(nc_fname):
    """
    C
    """
    nc_fname = normpath(nc_fname)
    nc_dset = Dataset(nc_fname)

    # check for OCHIDEE dataset
    # =========================
    ochidee_flag = True
    vars_not_prsnt = []
    for var in OCHIDEE_MANDAT_VARS:
        if var not in nc_dset.variables:
            vars_not_prsnt.append(var)

    if len(vars_not_prsnt) == 0:
        efiscen_flag = False
    else:
        ochidee_flag = False

        # check for EFISCEN dataset
        # =========================
        efiscen_flag = True
        vars_not_prsnt = []
        for var in EFISCEN_MANDAT_VARS:
            if var not in nc_dset.variables:
                vars_not_prsnt.append(var)

        if len(vars_not_prsnt) > 0:
            efiscen_flag = False

    nc_dset.close()

    if not ochidee_flag and not efiscen_flag:
        print(ERROR_STR + 'invalid dataset, could not be identified as OCHIDEE or EFISCEN')
    elif ochidee_flag:
        print('Dataset ' + nc_fname + ' is identified as OCHIDEE')
    else:
        print('Dataset ' + nc_fname + ' is identified as EFISCEN')

    return ochidee_flag, efiscen_flag

class OchideeSet(object, ):
    """
    Create object from an OCHIDEE NC file
    """
    def __init__(self, nc_fname, pfts, carbon_var):
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

        if 'veget' in nc_dset.variables:
            dset_type = 'OCHIDEE'
            nvegets = len(nc_dset.variables['veget'])
        else:
            dset_type = 'EFISCEN'
            nvegets = 0

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

                start_year = OCHIDEE_OVERRIDE_YEAR      # e.g. 1970
                print('\t' + WARN_STR + 'Start year set to: {}'.format(start_year))
                QApplication.processEvents()

                end_year = start_year + nyears - 1

            if var == 'time' and dset_type == 'EFISCEN':
                time_var = nc_dset.variables[var]
                nyears = len(time_var)
                start_year = time_var[0]
                end_year = time_var[-1]

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

        self.dset_type = dset_type
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
        if dset_type == 'EFISCEN':
            vals, aves, lookup_table = fetch_efiscen_vals(nc_dset, nyears, nlats, nlons, lookup_table)
        else:
            vals, aves, lookup_table = fetch_orchidee_vals(nc_dset, carbon_var, nvegets,
                                                           pfts, nyears, nlats, nlons, lookup_table)
        nc_dset.close()

        self.lookup_table = lookup_table
        self.vals = vals
        self.aves = aves
        self.nyears = nyears
        self.start_year = start_year
        self.end_year = end_year

    def get_ochidee_nc_data(self, pft_key, lat, long, baseline_flag):
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
            if pft_key is None:
                vals = nyears * [0]
            else:
                slice = self.vals[pft_key][:, lat_indx, lon_indx]
                is_masked = ma.is_masked(slice)
                if is_masked:
                    vals = ma.getdata(slice)
                else:
                    vals = list(slice)

            plnt_inpts['yrs'] = [yr for yr in range(strt_yr, strt_yr + nyears)]
            if baseline_flag:
                plnt_inpts['pis'] = nyears*[0]
            else:
                plnt_inpts['pis'] = [val for val in vals]
        else:
            plnt_inpts = None

        return plnt_inpts

def resize_yrs_pi(sim_strt_yr, sim_end_yr, yrs_pi):
    """
    patch to enable adjust yrs_pi to correspond to user specified simulation period
    """
    if yrs_pi is None:
        return None

    yr_frst = yrs_pi['yrs'][0]
    yr_last = yrs_pi['yrs'][-1]

    pi_frst = yrs_pi['pis'][0]
    pi_last = yrs_pi['pis'][-1]

    sim_yrs = list(range(sim_strt_yr, sim_end_yr + 1))

    sim_pis = []
    for iyr, yr in enumerate(sim_yrs):
        if yr < yr_frst:
            sim_pis.append(pi_frst)
        elif yr > yr_last:
            sim_pis.append(pi_last)
        else:
            sim_pis.append(yrs_pi['pis'][iyr])

    new_yrs_pi = {'yrs': sim_yrs, 'pis': sim_pis}

    return new_yrs_pi

def fetch_efiscen_vals(nc_dset, nyears, nlats, nlons, lookup_table):
    """
    Vegetation and corresponding plant functional type as defined in OCHIDEE model
    """
    vals = {}
    aves = {}
    pft_key = '00'
    tmp_vals = nc_dset.variables[EFISCEN_CARBON_VAR][:, :, :]
    vals[pft_key] = tmp_vals
    aves[pft_key] = tmp_vals.mean()

    for lat_indx in range(nlats):
        for lon_indx in range(nlons):
            n_masked = ma.count_masked(vals[pft_key][lat_indx, lon_indx, :])
            if n_masked == nyears:
                lookup_table[lat_indx][lon_indx] = False
            else:
                lookup_table[lat_indx][lon_indx] = True

        return vals, aves, lookup_table

def fetch_orchidee_vals(nc_dset, carbon_var, nvegets, pfts, nyears, nlats, nlons, lookup_table):
    """
    Vegetation and corresponding plant functional type as defined in OCHIDEE model
    """
    vals = {}
    aves = {}
    for pft_indx in range(nvegets):
        pft_key = '{0:0=2d}'.format(pft_indx + 1)
        tmp_vals = nc_dset.variables[carbon_var][:, pft_indx, :, :]  # TOTAL_BM_LITTER_c or TOTAL_LITTER_SOIL_c

        if tmp_vals.sum() == 0.0:
            print('\t' + WARN_STR + 'No data for vegetation type: ' + pft_key + ' PFT: ' + pfts[pft_key])

        if carbon_var == 'TOTAL_LITTER_SOIL_c':
            tmp_vals2 = nc_dset.variables['TOTAL_SOIL_c'][:, pft_indx, :, :]
            tmp_vals3 = tmp_vals - tmp_vals2
            tmp_vals = tmp_vals3  # TODO: check if necessary
            CNVRSN_FACT = CNVRSN_FACT_LITTER_SOIL  # gC/m**2 to kgC/ha
        else:
            CNVRSN_FACT = CNVRSN_FACT_BM_LITTER  # gC/m**2/day to kgC/ha/yr

        vals[pft_key] = CNVRSN_FACT * tmp_vals
        aves[pft_key] = CNVRSN_FACT * tmp_vals.mean()

        for lat_indx in range(nlats):
            for lon_indx in range(nlons):
                n_masked = ma.count_masked(vals[pft_key][:, lat_indx, lon_indx])
                if n_masked == nyears:
                    lookup_table[lat_indx][lon_indx] = False
                else:
                    lookup_table[lat_indx][lon_indx] = True

        return vals, aves, lookup_table

def orchidee_pfts():
    """
    Vegetation and corresponding plant functional type as defined in OCHIDEE model
    """
    pfts = {'01': 'SoilBareGlobal',
            '02': 'BroadLeavedEvergreenTropical',
            '03': 'BroadLeavedRaingreenTropical',
            '04': 'NeedleleafEvergreenTemperate',
            '05': 'BroadLeavedEvergreenTemperate',
            '06': 'BroadLeavedSummergreenTemperate',
            '07': 'NeedleleafEvergreenBoreal',
            '08': 'BroadLeavedSummergreenBoreal',
            '09': 'LarixSpBoreal',
            '10': 'C3GrassTemperate',
            '11': 'C4GrassTemperate',
            '12': 'C3AgricultureTemperate',
            '13': 'C4AgricultureTemperate',
            '14': 'C3GrassTropical',
            '15': 'C3GrassBoreal'}

    return pfts
