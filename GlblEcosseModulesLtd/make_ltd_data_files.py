#-------------------------------------------------------------------------------
# Name:        make_ltd_data_files.py
# Purpose:     based on Mark Richard's original but with some some mods
# Author:      soi698
# Created:     30/09/2011
# Copyright:   (c) soi698 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------

__version__ = '1.0.00'

from os.path import normpath, join
from calendar import month_abbr
from copy import copy
import validate

NoData = -999
MONTH_ABBRS = [mnth for mnth in month_abbr[1:]]
MaxNumFutureYears = 150
ERROR_STR = '*** Error *** '

class SoilLyr(object, ):

    def __init__(self, c, bulk_dens, ph, clay_pc, silt_pc, sand_pc, no_data=NoData):
        self.bulk_dens = bulk_dens
        self.ph = ph
        self.c = c     # C content
        self.clay_pc = clay_pc
        self.silt_pc = silt_pc
        self.sand_pc = sand_pc
        self.no_data = no_data

        return

    def validate(self):
        if self.c != self.no_data: validate.total_soil_carbon([self.c], depth='1 m')
        if self.bulk_dens != self.no_data: validate.bulk_density([self.bulk_dens])
        if self.ph != self.no_data: validate.soil_ph([self.ph])
        if self.clay_pc != self.no_data: validate.percent([self.clay_pc])
        if self.silt_pc != self.no_data: validate.percent([self.silt_pc])
        if self.sand_pc != self.no_data: validate.percent([self.sand_pc])
        total = 0
        for val in [self.clay_pc, self.silt_pc, self.sand_pc]:
            if val != self.no_data:
                total += val
        validate.percent([total])
##        assert(total <= 100.0)

        return

class MakeLtdDataFiles(object):
    """
    C
    """
    def __init__(self, form,  climgen, yrs_pi, comments=True, spacer_len=12, no_data = -999):
        """
        C
        """
        if hasattr(form, 'combo11s'):
            equil_mode = form.w_equimode.text()
        else:
            equil_mode = form.equimode

        if hasattr(form, 'combo11s'):
            sim_start_year = int(form.combo11s.currentText())
            sim_end_year = int(form.combo11e.currentText())
        else:
            sim_start_year = form.sim_strt_year
            sim_end_year = form.sim_end_year

        num_years = len(range(sim_start_year, sim_end_year + 1))
        if num_years > MaxNumFutureYears:
            print('Cannot specify more than {} of future years')
            return

        self.num_years = num_years

        # General options/settings
        # ========================
        self.comments = comments     # True = write comments, False = leave them out
        self.spacer_len = spacer_len  # Number of spaces between data and comment
        self.no_data = no_data
        self._luts = ['ara', 'gra', 'for', 'nat', 'mis', 'src']   # mksims values
        self.num_luts = len(self._luts)
        # self._elumluts =  ['ara', 'gra', 'for', 'nat', 'mis', 'src', 'sug', 'osr', 'srf', 'whe']
        self.months = MONTH_ABBRS

        # Parameters
        self.equil_mode = equil_mode           # Mode of equilibrium run
        self.num_lyrs = no_data
        self.lyr_depths = [] # list
        self.soil_lyrs = {}  # dict
        for lut in self._luts:
            self.soil_lyrs[lut] = []
        self.del_lyrs()   # is this necessary?
        self.plant_inputs = {}
        for lut in self._luts:
            self.plant_inputs[lut] = 0
        self.latitude = no_data
        self.wt_at_start = no_data
        if form.ecosse_exe is None or form.ecosse_exe == 'ecossev6_2c' or form.ecosse_exe == 'ecossev6_2b':
            self.wt_max_stand = no_data
        else:
            self.wt_max_stand = 0
        self.drain_class = no_data
        self.c_accum_b4_change = no_data
        self.ch4_b4_change = no_data
        self.co2_b4_change = no_data
        self.doc_loss_b4_change = no_data
        self.num_grow_seasons = no_data
        self.future_lu = []
        self.future_pi = []

        self.start_year = sim_start_year
        self.end_year = sim_end_year

        # create list of met file names
        # =============================
        self.met_fnames = []
        if climgen.ave_weather_flag:
            met_ave_file = climgen.met_ave_file
            for year in range(sim_start_year, sim_end_year + 1):
                self.met_fnames.append(met_ave_file.format(year))
        else:
            for year in range(sim_start_year, sim_end_year + 1):
                self.met_fnames.append('met{0}s.txt'.format(year))

        depths = form.depths
        ndepths = len(depths)
        if ndepths > 0:
            self.num_lyrs = ndepths
            for dep in depths:
                self.lyr_depths.append(dep)

        self.wt_at_start = 300
        self.drain_class = 2
        self.num_grow_seasons = num_years

        # generate landuse changes to be used in each lnput.txt file
        # ==========================================================
        # sorted(form.lu_pi_content.keys())
        #
        if yrs_pi is None:
            landuse_pi = copy(form.lu_pi_content['LandusePI'])
            landUses = []
            plantInputs = []
            num_years = self.end_year - self.start_year + 1
            for year_num in range(num_years):
                yr_str = str(year_num)
                if yr_str in landuse_pi.keys():
                    land_use, plant_input = landuse_pi[yr_str]
                    del landuse_pi[yr_str]

                land_use_num = form.land_use_types[land_use]
                landUses.append(land_use_num)
                plantInputs.append(plant_input)

            self.plantInput = plantInputs
            self.landUses = landUses
        else:
            # stanza for HoliSoils project
            # ============================
            lu_type = form.land_use_types['Forestry']
            self.landUses = len(yrs_pi['yrs'])*[lu_type]      # forest
            self.plantInput = yrs_pi['pis']
        pass

# ========================= end of init =================================================

    def add_lyr(self, lut_name, c_content, bulk_density, ph, clay_pc,
                silt_pc, sand_pc):
        self.soil_lyrs[lut_name.lower()].append(SoilLyr(c_content, bulk_density, ph,
                                        clay_pc, silt_pc, sand_pc, self.no_data))
        return

    def del_lyrs(self):
        for lut in self._luts:
            self.soil_lyrs[lut] = []

        return

    def line(self, data, comment):
        spacer_len = max(self.spacer_len - len(data), 2)
        spacer = ' ' * spacer_len
        return '{0}{1}# {2}\n'.format(data, spacer, comment)

    def validate(self):
        # Misc
        assert(0 < self.equil_mode < 10)
        validate.latitude([self.latitude])
        assert(0 <= self.wt_at_start <= 300)
        assert(self.num_grow_seasons > 0)

        # Soil layers
        assert(0 < self.num_lyrs < 11)
        for i in range(len(self.lyr_depths)):
            assert(0 < self.lyr_depths[i] < 300)
            assert(self.lyr_depths[i] % 5 == 0)  # Depth should be a multiple of 5
            if i > 0:   # If not the top layer check that its depth is greater than the layer above
               assert(self.lyr_depths[i] > self.lyr_depths[i-1])
        for lut in self._luts:
            for lyr_num, lyr in enumerate(self.soil_lyrs[lut]):
                lyr.validate()
                assert(lyr_num < self.num_lyrs)

        # Plant inputs
        if self.equil_mode in [1,3]:
            for key in self._luts:
                assert(0 <= self.plant_inputs[key] < 20000) # Upper limit?
        validate.plant_c_inputs(self.future_pi, 'annual')

        # Land use
        for lu in self.future_lu:
            assert(1 <= lu <= len(self._elumluts))

        # Climate'
        # validate.rainfall(self.lta_precip, 'monthly')
        # validate.self.lta_precip (self.lta_tmean)

        assert(len(self.future_lu) == self.num_grow_seasons)
        assert(len(self.future_pi) == self.num_grow_seasons)

        return

    def write(self, sim_dir, soil, latitude, hist_weather_recs, met_rel_path, input_fname='input.txt'):
        """
        MJM: this function has been hacked around from mksims original

        previously:
            # Read the file comprising historic precipitation and temperature
            with open(hist_historic_fname, 'r') as finp:
                hist_historic_lines = finp.readlines()
        """
        type_soil = type(soil)
        if type_soil is not list:
            print(ERROR_STR + 'Problem writing {} for simulation directory {} - soil has type {} should have list type'
                  .format(input_fname, sim_dir, type_soil))
            return

        if self.landUses is None or self.plant_inputs is None:
            print(ERROR_STR + 'Land uses and/or plant inputs are None - cannot write meaningful input file')
            return

        self.latitude = latitude

        if len(soil) == 7:
            num_lyrs = 1
            lyr_depths = self.lyr_depths[0:1]
        else:
            num_lyrs = self.num_lyrs
            lyr_depths = self.lyr_depths

        # Soil parameters
        # ===============
        output_buff = []
        output_buff.append(self.line('{0}'.format(self.equil_mode), 'Mode of equilibrium run'))
        output_buff.append(self.line('{0}'.format(num_lyrs), 'Number of soil layers (max 10)'))
        for lyr_num, lyr_depth in enumerate(lyr_depths):
            output_buff.append(self.line('{0}'.format(lyr_depth), 'Depth of bottom of SOM layer {0} [cm]'.format(lyr_num+1)))

        # MJM: modified to be cruder than MR original
        # for each LUT write soil characteristics for each soil layer
        # ===========================================================
        for key in self._luts:      # typically: ['ara', 'gra', 'for', 'nat', 'mis', 'src']
            for lyr_num in range(num_lyrs):    # typically: 2
                strt_indx = 6*lyr_num
                output_buff.append(self.line('{}'.format(soil[strt_indx]), 'C content [kgC/ha] for this soil under {0} in SOM layer {1}'.format(key, lyr_num+1)))
                output_buff.append(self.line('{}'.format(soil[strt_indx + 1]), 'Bulk density [g/cm3] for this soil under {0} in SOM layer {1}'.format(key, lyr_num+1)))
                output_buff.append(self.line('{}'.format(soil[strt_indx + 2]), 'pH for this soil under {0} in SOM layer {1}'.format(key, lyr_num+1)))
                output_buff.append(self.line('{}'.format(soil[strt_indx + 3]), '% clay by weight for this soil under {0} in SOM layer {1}'.format(key, lyr_num+1)))
                output_buff.append(self.line('{}'.format(soil[strt_indx + 4]), '% silt by weight for this soil under {0} in SOM layer {1}'.format(key, lyr_num+1)))
                output_buff.append(self.line('{}'.format(soil[strt_indx + 5]), '% sand by weight for this soil under {0} in SOM layer {1}'.format(key, lyr_num+1)))

        # Long term average plant C input
        # ===============================
        for key in self._luts:
            output_buff.append(self.line('{}'.format(self.plant_inputs[key]),
                        '{} long term average plant C input [kgC/ha/yr] (obsolete - use a dummy value)'.format(key)))

        # Long term average climate - 24 records
        # ======================================
        for weather_rec in hist_weather_recs:
            output_buff.append(weather_rec)

        # Other parameters
        # ================
        output_buff.append(self.line('{0}'.format(round(latitude,3)), 'Latitude [decimal deg]'))
        output_buff.append(self.line('{0}'.format(self.wt_at_start), 'Water table depth at start [cm]'))
        if self.wt_max_stand >= 0:  # ensures backward compatibility
            output_buff.append(self.line('{0}'.format(self.wt_max_stand), 'Max standing water [cm]'))
        output_buff.append(self.line('{0}'.format(self.drain_class), 'Drainage class (not yet used)'))
        output_buff.append(self.line('{0}'.format(self.c_accum_b4_change), 'C accumulated before change [kgC/ha/yr] (obsolete - use a dummy value)'))
        output_buff.append(self.line('{0}'.format(self.ch4_b4_change), 'CH4 emission before change [kgC/ha/yr] (not used yet)'))
        output_buff.append(self.line('{0}'.format(self.co2_b4_change), 'CO2 emission before change [kgC/ha/yr] (not used yet)'))
        output_buff.append(self.line('{0}'.format(self.doc_loss_b4_change), 'DOC loss before change [kgC/ha/yr] (not used yet)'))
        output_buff.append(self.line('{0}'.format(self.num_grow_seasons), 'Number of growing seasons to simulate'))

        # Future land use and plant inputs
        # ================================
        # NB "(if plant input set to zero it is obtained from RothC instead)" has been omitted to reduce size of file
        num_years = self.end_year - self.start_year  + 1
        for year_num in range(num_years):
            lu = self.landUses[year_num]
            plnt_inpt = round(self.plantInput[year_num], 3)
            rec = '{}, {}'.format(lu, plnt_inpt)
            comment = 'Year {} land use code and plant C input [kgC/ha/yr]'.format(year_num)
            output_buff.append(self.line(rec, comment))

        # Climate file names
        # ==================
        for year_num, fname in enumerate(self.met_fnames):
            output_buff.append(self.line('{}{}'.format(met_rel_path, fname), 'Year {0} climate file'.format(year_num+1)))

        path_input_txt = join(normpath(sim_dir), input_fname)
        try:
            fhand = open(path_input_txt, 'w')
        except IOError:
            raise IOError('Unable to open file {0}'.format(path_input_txt))
        else:
            fhand.writelines(output_buff)
            fhand.close()

        return
