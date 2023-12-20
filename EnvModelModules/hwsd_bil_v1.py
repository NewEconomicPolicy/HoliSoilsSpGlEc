#-------------------------------------------------------------------------------
# Name:        hwsd_bil_v1.py
# Purpose:     Class to read HWSD .hdr and .bil files and return a 2D array comprising MU_GLOBAL values for a given
#              bounding box
#               input files are: hwsd.hdr and hwsd.bil
#               notes:
#                   a) lat longs are stored in unites of 30 arc seconds for sake of accuracy and simplicity - so we
#                       convert lat/longs read from shapefiles to nearest seconds using round()
#
#                   b) output file names take the form: mu_global_NNNNN.nc where NNNNN is the region name??
#
#                   c) modified on 03-09-2016 to include soils which have only a top layer
#
# Author:      Mike Martin
# Created:     16/09/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

__prog__ = 'hwsd_bil_v1.py'
__version__ = '0.0.0'

#
from csv import reader, DictReader
from os.path import join, exists
from sys import exit, stdout
from struct import unpack
from locale import LC_ALL, setlocale, format_string
from numpy import arange, dtype, zeros, int32
from time import sleep
import math

sleepTime = 5

def validate_hwsd_rec (lgr, mu_global, data_rec):
    """
    # function to make sure we are supplying a valid HWSD data record for a given mu_global

    # TODO: improve this function by using a named tuple:
                 data_rec = list([
                 0 row['SHARE'],
                 1 row['T_SAND'],
                 2 row['T_SILT'],
                 3 row['T_CLAY'],
                 4 row['T_BULK_DENSITY'],
                 5 row['T_OC'],
                 6 row['T_PH_H2O'],
                 7 row['S_SAND'],
                 8 row['S_SILT'],
                 9 row['S_CLAY'],
                10 row['S_BULK_DENSITY'],
                11 row['S_OC'],
                12 row['S_PH_H2O']
    """
    share, t_sand, t_silt, t_clay, t_bulk_density, t_oc, t_ph_h2o, \
        s_sand, s_silt, s_clay, s_bulk_density, s_oc, s_ph_h2o = data_rec

    # topsoil is mandatory
    ts_rec = data_rec[1:7]
    for val in ts_rec:
        if val == '':
            # log a message
            mess =  'Incomplete topsoil data - mu_global: ' + str(mu_global) + '\trecord: ' + ', '.join(data_rec)
            # print(mess)
            lgr.info(mess)
            return False

    '''
    MJM: borrowed from mksims.py:
    Convert top and sub soil Organic Carbon percentage weight to kgC/ha as follows:
                      / 100: % to proportion,
                     * 100000000: cm3 to hectares,
                      /1000: g to kg,
                     *30: cm to lyr depth
    '''
    # lgr.info('Soil data check for mu_global {}\tshare: {}%'.format(mu_global,share))
    t_bulk = float(t_bulk_density)  # T_BULK_DENSITY
    t_oc = float(t_oc)    # T_OC
    t_c = round(t_bulk * t_oc / 100.0 * 30 * 100000000 / 1000.0, 1)

    # sub soil
    ss_rec = data_rec[7:]
    for val in ss_rec:
        if val == '':
            # log a message
            lgr.info('Incomplete subsoil data - mu_global: ' + str(mu_global) + '\trecord: ' + ', '.join(data_rec))
            # modified share to from int to float
            return list([t_c, t_bulk, float(t_ph_h2o), int(t_clay), int(t_silt), int(t_sand), float(share)])

    s_bulk = float(s_bulk_density) # S_BULK_DENSITY
    s_oc = float(s_oc)   # S_OC
    s_c = round(s_bulk * s_oc / 100.0 * 70 * 100000000 / 1000.0, 1)

    # Two layers: C content, bulk, PH, clay, silt, sand
    try:
        # modified share to from int to float
        ret_list = list([t_c, t_bulk, float(t_ph_h2o), int(t_clay), int(t_silt), int(t_sand),
                 s_c, s_bulk, float(s_ph_h2o), int(s_clay), int(s_silt), int(s_sand), float(share)])
    except ValueError as e:
        ret_list = []
        print('problem {} with mu_global {}\tdata_rec: {}'.format(e, mu_global, data_rec))
        sleep(sleepTime)
        exit(0)

    return ret_list

class HWSD_bil(object,):

    def __init__(self,form):

        # open header file and read first 9 lines
        # ========================
        inpfname = 'hwsd.hdr'
        inp_file = join(form.hwsd_dir, inpfname)
        finp = open(inp_file, 'r')
        finp_reader = reader(finp)

        row = next(finp_reader); dummy, byteorder = row[0].split()
        row = next(finp_reader)   # skip: LAYOUT       BIL

        row = next(finp_reader); dummy, f1 = row[0].split(); nrows = int(f1)
        row = next(finp_reader); dummy, f1 = row[0].split(); ncols = int(f1)
        row = next(finp_reader); dummy, f1 = row[0].split(); nbands = int(f1)
        row = next(finp_reader); dummy, f1 = row[0].split(); nbits = int(f1)


        row = next(finp_reader); dummy, f1 = row[0].split(); bandrowbytes = int(f1)
        row = next(finp_reader); dummy, f1 = row[0].split(); totalrowbytes = int(f1)
        row = next(finp_reader); dummy, f1 = row[0].split(); bandgapbytes = int(f1)

        finp.close()

        # the HWSD grid covers the globe's land area with 30 arc-sec grid-cells
        ngranularity = 120

        self.hwsd_dir = form.hwsd_dir
        self.lgr = form.lgr
        self.nrows = nrows
        self.ncols = ncols

        # wordsize is in bytes
        self.wordsize = round(nbits/8)
        self.granularity = float(ngranularity)

        # these attrubutes change according to the bounding box
        self.rows = arange(1)   # create ndarray
        self.mu_globals = []
        self.badCntr = 0
        self.zeroCntr = 0

        # these define the extent of the grid
        self.nlats = 0; self.nlons = 0
        self.nrow1 = 0; self.nrow2 = 0
        self.ncol1 = 0; self.ncol2 = 0

        # this will store dictionary of mu_global keys and corresponding soil data
        # will be refreshed each time user selects a new bounding box
        self.soil_recs = {}

        # this section is adapted from MR's hwsd_uk.py
        # ===========================================
        '''
        # self.nan = nan

        # self.Rec = namedtuple('Rec', ('id', 'lon', 'lat', 'soils'))
        self.Rec = namedtuple('Rec', ('mu_global', 'soils'))
        # define a tuple
        soil_fields = ('share', 't_sand', 't_silt', 't_clay', 't_bulk', 't_oc',
                       't_ph', 's_sand', 's_silt', 's_clay', 's_bulk', 's_oc',
                       's_ph')
        self.soil_fields = soil_fields
        self.nsoil_fields = len(soil_fields)
        self.Soil = namedtuple('Soil', soil_fields)
        '''

    def read_bbox_hwsd_mu_globals(self, bbox, hwsd_mu_globals, upscale_resol):

        # this function creates a grid of MU_GLOBAL values corresponding to a given
        # bounding box defined by two lat/lon pairs

        func_name =  __prog__ + '.read_bbox_hwsd_mu_globals'
        self.lgr.info('Running programme ' + func_name)
        setlocale(LC_ALL, '')

        # the HWSD grid covers the globe's land area with 30 arc-sec grid-cells
        # AOI is typically county sized e.g. Isle of Man
        granularity = self.granularity

        # these are lat/lon
        upper_left_coord = [bbox[3],bbox[0]]
        lower_right_coord = [bbox[1],bbox[2]]
        self.upper_left_coord = upper_left_coord
        self.lower_right_coord = lower_right_coord

        # round these so that they are divisable by the requested resolution
        nlats = round((upper_left_coord[0] - lower_right_coord[0])*granularity)
        nlats = upscale_resol*math.ceil(nlats/upscale_resol)
        nlons = round((lower_right_coord[1] - upper_left_coord[1])*granularity)
        nlons = upscale_resol*math.ceil(nlons/upscale_resol)

        # construct 2D array
        rows = zeros(nlats*nlons, dtype = int32)
        rows.shape = (nlats,nlons)
        #
        # work out first and last rows
        nrow1 = round((90.0 - upper_left_coord[0])*granularity)
        nrow2 = nrow1 + nlats
        ncol1 = round((180.0 + upper_left_coord[1])*granularity)
        ncol2 = ncol1 + nlons

        # read a chunk the CSV file
        # =========================
        aoi_chunk = hwsd_mu_globals.data_frame.loc[hwsd_mu_globals.data_frame['gran_lat'].isin(range(nrow1, nrow2))]
        nvals_read = len(aoi_chunk)
        nrejects = 0
        for index, rec in aoi_chunk.iterrows():
            icol = int(rec[1])
            irow = int(rec[0])
            mu_global = int(rec[2])

            # make sure indices are kept within limits
            # ========================================
            col_indx = icol - ncol1
            if col_indx < 0 or col_indx >= nlons:
                nrejects += 1
            else:
                rows[irow-nrow1][icol-ncol1] = mu_global

        nvals_read_str = format_string('%d', nvals_read, grouping=True)
        nrejects_str = format_string('%d', nrejects, grouping=True)
        mess = 'Read ' + nvals_read_str + ' cells from hwsd mu globals data frame, rejected ' + nrejects_str
        self.lgr.info(mess)

        self.nrow1 = nrow1
        self.nrow2 = nrow2
        self.ncol1 = ncol1
        self.ncol2 = ncol2

        self.rows = rows
        self.nlats = nlats
        self.nlons = nlons

        mess =  'Retrieved {} mu globals for AOI of {} lats and {} lons'.format(nvals_read, nlats, nlons)
        self.lgr.info(mess)
        nretrieve_str = format_string('%d', nvals_read - nrejects, grouping=True)
        print('Found mu_globals for ' + nretrieve_str + ' cells')

        return nvals_read

    def get_soil_recs(self,keys):
        """
        # get soil records associated with each MU_GLOBAL list entry
        # assumption: the keys passed are already sorted,
        # also and crucially: that in the .csv file, the mu_globals are ordered from lowest value to highest
        """
        func_name =  __prog__ + ' get_soil_recs'

        inpfname = 'HWSD_DATA.csv'
        csv_file = join(self.hwsd_dir, inpfname)
        if not exists(csv_file):
            stdout.write('Error in  get_soil_recs file does not exist: ' + csv_file + ' cannot proceed \n')
            return -1

        # build a dictionary with mu_globals as keys
        soil_recs = {}

        kiters = iter(keys)
        key = next(kiters)
        if key == 0:
            key = next(kiters)
        else:
            # initialise dictionary entry - TODO: is this necessary?
            soil_recs[key] = []

        # initialisations
        mu_glob_prev = -999

        data_recs = []

        with open(csv_file) as csvf:
            csvreader = DictReader(csvf)

            self.lgr.info('csv file of mu_global data: {0}\t dialect: {1}\t\tnum fields: {2}'
                  .format(csv_file,csvreader.dialect,len(csvreader.fieldnames)))

            # step through rows in CSV file until we pick up all data
            # soil_data_recs = []

            for row in csvreader:
                # output required columns - see hwsd_uk.txt in mksims.py
                # mu_global,share1,t_sand1,t_silt1,t_clay1,t_bulk1,t_oc1,t_ph1,s_sand1,s_silt1,s_clay1,s_bulk1,s_oc1,s_ph1
                # all floats - then fill with nans...
                mu_global = int(row['MU_GLOBAL'])

                # if mu_global changes then record accumulated soil data
                if mu_global != mu_glob_prev:

                    # if soil data has been accumulated then record it
                    if len(data_recs) > 0:
                        soil_recs[mu_glob_prev] = data_recs
                        data_recs = []

                    #
                    mu_glob_prev = mu_global

                if mu_global > key:
                    try:
                        key = next(kiters)
                    except StopIteration:
                        # once list of keys (mu_globals) is exhausted then close csv file and break out
                        csvf.close()
                        break

                if mu_global == key:
                    data_rec = list([row['SHARE'],row['T_SAND'],row['T_SILT'],row['T_CLAY'],
                                   row['T_BULK_DENSITY'],row['T_OC'],row['T_PH_H2O'],
                                   row['S_SAND'],row['S_SILT'],row['S_CLAY'],
                                   row['S_BULK_DENSITY'],row['S_OC'],row['S_PH_H2O']])
                    adjusted_datarec = validate_hwsd_rec(self.lgr, mu_global, data_rec)
                    if adjusted_datarec != False:
                        data_recs.append(adjusted_datarec)
                    self.lgr.handlers[0].flush()
                    # os.fsync(self.lgr.handlers[0].fileno())

        # End of main loop - make doubly sure that the CSV file is closed
        if not csvf.closed:
            self.lgr.info('Warning: csv file of mu_global data: {} closed'.format(csv_file))
            csvf.close()

        # if soil data has been accumulated then record it
        if len(data_recs) > 0:
             soil_recs[mu_glob_prev] = data_recs

        self.lgr.info('Exiting function {0} after generating {1} records from soil table {2}\n'
                                        .format(func_name, len(soil_recs), csv_file))

        # compare keys and create a list of mu_globals which are missing data in the HWSD
        # TODO: is this good enough?
        bad_muglobals = []
        for key in keys:
            if soil_recs.get(key) == None:
                bad_muglobals.append(key)
            elif len(soil_recs[key]) == 0:
                bad_muglobals.append(key)
        self.bad_muglobals = bad_muglobals

        return soil_recs

    def read_bbox_mu_globals(self, bbox, snglPntFlag = False):

        # this function creates a grid of MU_GLOBAL values corresponding to a given
        # bounding box defined by two lat/lon pairs

        func_name =  __prog__ + ' read_bbox_mu_globals'
        self.lgr.info('Running programme ' + func_name)

        # the HWSD grid covers the globe's land area with 30 arc-sec grid-cells
        # AOI is typically county sized e.g. Isle of Man

        granularity = self.granularity
        ncols = self.ncols
        wordsize = self.wordsize

        if snglPntFlag:
            nlats = 1
            nlons = 1
            lower_left_coord = [bbox[1],bbox[0]]
            nrow1 = round((90.0 - lower_left_coord[0])*granularity)
            nrow2 = nrow1
            ncol1 = round((180.0 + lower_left_coord[1])*granularity)
            ncol2 = ncol1
            nlats = 1
            nlons = 1
        else:
            # these are lat/lon
            upper_left_coord = [bbox[3],bbox[0]]
            lower_right_coord = [bbox[1],bbox[2]]
            self.upper_left_coord = upper_left_coord
            self.lower_right_coord = lower_right_coord

            nlats = round((upper_left_coord[0] - lower_right_coord[0])*granularity)
            nlons = round((lower_right_coord[1] - upper_left_coord[1])*granularity)
            #
            # work out first and last rows
            nrow1 = round((90.0 - upper_left_coord[0])*granularity) + 1
            nrow2 = round((90.0 - lower_right_coord[0])*granularity)
            ncol1 = round((180.0 + upper_left_coord[1])*granularity) + 1
            ncol2 = round((180.0 + lower_right_coord[1])*granularity)

        chunksize = wordsize*nlons
        form = str(int(chunksize/2)) + 'H'  # format for unpack

        # TODO - construct a 2D array from extraction
        inpfname = 'hwsd.bil'
        inp_file = join(self.hwsd_dir, inpfname)
        finp = open(inp_file, 'rb')

        nrows_read = 0
        rows = arange(nlats*nlons)
        rows.shape = (nlats,nlons)

        # read subset of HWSD .bil file
        # =============================
        for nrow in range(nrow1, nrow2 + 1):
            offst = (ncols*nrow + ncol1)*wordsize

            finp.seek(offst)
            chunk = finp.read(chunksize)
            if chunk:
                rows[nrows_read:] = unpack(form, chunk)
                # build array of mu_globals
                nrows_read += 1
                # print(' row {0:d}\tform {1:s}'.format(nrow,form))
            else:
                break

        finp.close()

        self.nrow1 = nrow1
        self.nrow2 = nrow2
        self.ncol1 = ncol1
        self.ncol2 = ncol2

        self.rows = rows
        self.nlats = nlats
        self.nlons = nlons

        self.lgr.info('Exiting function {0} after reading {1} records from .bil file\n'.format(func_name,nrows_read))
        return nrows_read

    def mu_global_list(self):

        # build a list of mu_globals from grid

        # reshape
        mu_globals = []
        nlats = self.nlats
        nlons = self.nlons

        # reshape tuple
        self.rows.shape = (nlats*nlons)
        for mu_global in self.rows:
            if mu_global not in mu_globals:
                mu_globals.append(mu_global)

        # revert shape of tuple
        self.rows.shape = (nlats,nlons)

        self.mu_globals = mu_globals

        return len(mu_globals)

    def get_mu_globals_dict(self):

        func_name =  __prog__ + ' get_mu_globals_dict'
        # function returns list of MU_GLOBALs in grid with number of occurrences of each MU_GLOBAL value
        # i.e. build dictionary of MU_GLOBALs from grid

        mu_globals = {}     # initialise dictionary

        nlats = self.nlats
        nlons = self.nlons

        for irow in range(0,nlats):

            val_rec = []
            for icol in range(0,nlons):

                mu_global = self.rows[irow,icol]

                if mu_global not in mu_globals.keys():
                    # add detail
                    mu_globals[mu_global] = 1
                else:
                    mu_globals[mu_global] += 1

        self.lgr.info('Func: {0}\tUnique number of mu_globals: {1}'.format(func_name, len(mu_globals)))
        for key in sorted(mu_globals.keys()):
            self.lgr.info('\t' + str(key) + ' :\t' + str(mu_globals[key]))

        # filter condition where there is only sea, i.e. mu_global == 0
        if len(mu_globals.keys()) == 1:
            for key in mu_globals.keys():
                if key == 0:
                    mu_globals = {}
                    print('Only one mu_global with value of zero - nothing to process')

        return mu_globals

    def get_data_mapit(self):

        # function to return lists which can be plotted using mapit
        xcoords = []  # e.g. longitude or easting
        ycoords = []  # e.g. latitude or northing
        vals = []
        nrow1 = self.nrow1
        nrow2 = self.nrow2
        nlats = self.nlats

        ncol1 = self.ncol1
        ncol2 = self.ncol2
        nlons = self.nlons

        yc = self.nrow1
        for irow in range(0,nlats):
            xc = self.ncol1
            for icol in range(0,nlons):
                xcoords.append(xc)
                ycoords.append(yc)
                vals.append(self.rows[irow,icol])
                xc += 1
            yc += 1

        return xcoords, ycoords, vals