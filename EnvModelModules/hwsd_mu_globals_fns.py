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

__prog__ = 'hwsd_mu_globals_fns.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

from os.path import isdir, join, normpath, splitext, exists, lexists, isfile, split
from locale import setlocale, LC_ALL, format_string
from pandas import read_csv
from numpy import count_nonzero, zeros, int32, float64
from copy import copy
from shape_funcs import calculate_area

ERROR_STR = '*** Error *** '
null_value = -9999
big_number = 999.0
nonzero_proportion = 0.75  # threshold proportion of non-zero mu globals for a meta-cell below which centroid is calculated

def meta_cell_centroid(aggre_cell, irow, icol, num_non_zeros, resol_upscale, granularity):

    # work out centroid if there is a significant proportion of sea ie zero values
    irow_sum = 0
    icol_sum = 0
    for ir in range(resol_upscale):
        for ic in range(resol_upscale):
            if aggre_cell[ir,ic] > 0:
                irow_sum += ir
                icol_sum += ic

    irow_ave = irow_sum/num_non_zeros
    icol_ave = icol_sum/num_non_zeros
    # print('irow/icol ave: {} {}'.format(irow_ave, icol_ave))

    gran_lat = irow + irow_ave
    cell_lat = 90.0 - gran_lat/granularity
    gran_lon = icol + icol_ave
    cell_lon = gran_lon/granularity - 180.0

    return gran_lat, gran_lon, cell_lat, cell_lon

def process_aggregated_cell(irow, icol, hwsd, resol_d2, resol_upscale, aggre_cell):
    '''
    process aggregated cell (array aggre_cell) into a single representative cell (the meta-cell)
    filter out zero mu globals and those without data
    '''
    func_name =  __prog__ + ' process_aggregated_cell'

    granularity = hwsd.granularity

    # build dictionary of MU_GLOBALs with number of occurrences of each MU_GLOBAL value
    # =================================================================================
    mu_globals = {}     # initialise dictionary
    for mu_global in aggre_cell.ravel():

        # exclude elements with zero or bad mu_globals
        # =============================================
        if mu_global in hwsd.bad_muglobals:
            if mu_global == 0:
                # zero mu_globals correspond to sea and areas clipped out
                hwsd.zeroCntr += 1
            else:
                hwsd.badCntr += 1
        else:
            #  create dictionary entry or increment existing entry
            if mu_global not in mu_globals.keys():
                mu_globals[mu_global] = 1
            else:
                mu_globals[mu_global] += 1


    mu_globals_props = {int(key):round(val/aggre_cell.size, 5) for (key,val) in mu_globals.items()}

    # work out percentage coverage for each mu_global
    # ===============================================
    num_non_zeros = count_nonzero(aggre_cell)
    coverage_ratio = num_non_zeros/aggre_cell.size
    if coverage_ratio < nonzero_proportion:
        gran_lat, gran_lon, cell_lat, cell_lon = \
                            meta_cell_centroid(aggre_cell, irow, icol, num_non_zeros, resol_upscale, granularity)
    else:
        gran_lat = irow + resol_d2
        cell_lat = 90.0 - gran_lat/granularity
        gran_lon = icol + resol_d2
        cell_lon = gran_lon/granularity - 180.0

    # calculate cell area
    # ===================
    lat_ur = 90.0 - irow/granularity
    lat_ll = 90.0 - (irow + resol_upscale)/granularity
    lon_ll = icol/granularity - 180.0
    lon_ur = (icol + resol_upscale)/granularity - 180.0
    bbox = list([lon_ll, lat_ll, lon_ur, lat_ur])
    area = calculate_area(bbox)
    coverage = area*coverage_ratio

    return  int(gran_lat), int(gran_lon), cell_lat, cell_lon, coverage, mu_globals_props

def gen_grid_cells_for_band(hwsd, resol_upscale):
    """
    # resolution is an integer multiple where 1 corresponds to the HWSD i.e. 30 arc seconds, so a resolution of 12 would
    # correspond to 12x30 arc seconds i.e. 6 arc minutes
    # this function associates each soil grid point with the proximate gridded climate data
    # at the time of writing (Dec 2015) soil data is on 30 second grid whereas climate data is half degree.
    """
    func_name =  __prog__ + ' gen_grid_cells_for_band'

    # these are lat/lon
    # determine the largest bbox enclosing the province of interest
    lon_ll_aoi = big_number
    lat_ll_aoi = big_number
    lon_ur_aoi = -big_number
    lat_ur_aoi = -big_number

    # now take each soil coord (integers) and store key to pettmp data
    nrow1 = hwsd.nrow1
    nrow2 = hwsd.nrow2
    nlats = hwsd.nlats

    ncol1 = hwsd.ncol1
    ncol2 = hwsd.ncol2
    nlons = hwsd.nlons

    # move north to south and west to east (decreasing latitude and longitude)
    # build list of mu_globals with their coordinates (granular and lat lon)
    # when mu_globals for the aggrecell (aggregated cell) have been collected then append aggrecell to soil_res
    soilres = []
    aggregated_cell = zeros(resol_upscale*resol_upscale, dtype = int32)
    aggregated_cell.shape = (resol_upscale, resol_upscale)
    resol_d2 = int(resol_upscale/2)
    ixhws1 = 0
    iyhws1 = 0
    for irow in range(nrow1, nrow2, resol_upscale):
        iyhws2 = iyhws1 + resol_upscale

        for icol in range(ncol1, ncol2, resol_upscale):

            ixhws2 = ixhws1 + resol_upscale
            # NB. numbers, strings, tuples, frozensets, booleans, None, are immutable
            # lists, dictionaries, sets, bytearrays, are mutable
            aggregated_cell = copy(hwsd.rows[iyhws1:iyhws2,ixhws1:ixhws2])

            # skip processing if mu_global is 0 i.e. comprises water only
            # =================================
            if aggregated_cell.max() > 0:
                # record this ac using middle point
                grid_cell = process_aggregated_cell(irow, icol, hwsd, resol_d2, resol_upscale, aggregated_cell)
                gran_lat, gran_lon, cell_lat, cell_lon, area, mu_global_props = grid_cell
                soilres.append(grid_cell)

                # is this useful?
                # ===============
                lon_ll_aoi = min(lon_ll_aoi, cell_lon)
                lat_ll_aoi = min(lat_ll_aoi, cell_lat)
                lon_ur_aoi = max(lon_ur_aoi, cell_lon)
                lat_ur_aoi = max(lat_ur_aoi, cell_lat)

            # end of inner loop
            # ================
            ixhws1 = ixhws2
        # end of outer loop
        # ================
        ixhws1 = 0
        iyhws1 = iyhws2

    hwsd.lgr.info('Finished in function {}\tgenerated {} soil elements\tdiscarded {} zero or bad mu_globals'
                  .format(func_name,len(soilres),(hwsd.zeroCntr + hwsd.badCntr)))
    hwsd.lgr.info('\tgrid points summary - zero mu_global: {}\tbad non-zero: {}\tbad mu_globals: {}'
                  .format(hwsd.zeroCntr, hwsd.badCntr, hwsd.bad_muglobals))

    bbox = list([round(lon_ll_aoi,4), round(lat_ll_aoi,4), round(lon_ur_aoi,4), round(lat_ur_aoi,4)])
    return soilres, bbox

class HWSD_mu_globals_csv(object,):

    def __init__(self, form, csv_fname):

        '''
        this object uses pandas to read a file comprising granular lat, granular lon, mu_global, latitude, longitude
        and ascertains its bounding box
        '''
        class_name =  'HWSD_aoi_mu_globals_csv'
        mess = 'Creating class ' + class_name + ' from CSV file ' + csv_fname
        form.lgr.info(mess); print(mess + '\n\t- this may take some time...')

        # read the CSV file
        # =================
        headers = ['gran_lat', 'gran_lon', 'mu_global', 'latitude', 'longitude', 'land_use']
        try:
            data_frame = read_csv(csv_fname, sep = ',', names = headers, dtype={'gran_lat': int32, 'gran_lon': int32,
                                                        'mu_global': int32, 'latitude': float64, 'longitude': float64})
        except ValueError as err:
            form.success_flag = False
            print(ERROR_STR + str(err) + ' reading ' + csv_fname)
        else:
            form.success_flag = True
            nlines = len(data_frame)
            setlocale(LC_ALL, '')
            nlines_str = format_string("%d", nlines, grouping=True )

            # sort values North to South and West to East
            # ===========================================
            data_frame = data_frame.sort_values(by=["gran_lat", "gran_lon"])

            # determine the largest bbox enclosing the AOI
            lon_ll_aoi = data_frame['longitude'].min()
            lat_ll_aoi = data_frame['latitude'].min()
            lon_ur_aoi = data_frame['longitude'].max()
            lat_ur_aoi = data_frame['latitude'].max()

            mess = 'read and sorted {} lines using pandas'.format(nlines_str)
            form.lgr.info(mess); print(mess)

            mu_global_dict = dict(data_frame['mu_global'].value_counts())

            aoi_label = 'LL: {:.2f} {:.2f}\tUR: {:.2f} {:.2f}\t# cells: {}' \
                                    .format(lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi, nlines_str)
            mess = 'Exiting function {} AOI: {}\tCSV file {}\n'.format(class_name, aoi_label, csv_fname)
            form.lgr.info(mess)

            # list file contents
            # ==================
            mu_global_list = list(mu_global_dict.keys())
            dummy, fname_short = split(csv_fname)
            mess = 'HWSD csv file {} has {} unique mu globals'.format(fname_short, len(mu_global_list))
            print(mess);         form.lgr.info(mess)
            # self.soil_recs = hwsd.get_soil_recs(sorted(mu_global_dict.keys())) # sorted key list (of mu_globals)

            # complete object
            # ===============
            self.lon_ll_aoi = lon_ll_aoi
            self.lat_ll_aoi = lat_ll_aoi
            self.lon_ur_aoi = lon_ur_aoi
            self.lat_ur_aoi = lat_ur_aoi
            self.nlines = nlines
            self.mu_global_list = sorted(mu_global_list)
            self.aoi_label = aoi_label
            self.csv_fname = csv_fname
            self.data_frame = data_frame
