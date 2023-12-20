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

__prog__ = 'plant_input_csv_fns.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

import os
from pandas import read_csv, DataFrame
from csv import reader

REQUIRED_FIELDS = list(['mu_global', 'gran_lat', 'gran_lon'])
WARN_STR = '*** Warning *** '
sleepTime = 5
max_lines = 10

def associate_plant_inputs(lggr_info, strt_year, gran_lat, gran_lon, pi_df, ltd_data, write_to_lggr = False):
    """
    pi_df is plant input data frame comprising monthly PI values therefore need to sum for each year
    then modify ltd_data object with plant inputs
    """
    pi_ann = {}
    locat = pi_df.loc[(pi_df['gran_lat'] == gran_lat) & (pi_df['gran_lon'] == gran_lon)]
    if len(locat.values) == 0:
        # print('No data at cordinates: gran_lat: {}\tgran_lon: {}'.format(gran_lat, gran_lon))
        return None

    # location exists - sum values
    # ============================
    len_locat = len(locat.values[0])
    yr = strt_year
    for indx in range(5, len_locat, 12):
        pi_ann[str(yr)] = locat.values[0][indx:indx + 12].sum()
        yr += 1

    if write_to_lggr:
        lat = -999.0
        lon = -999.0
        mess = 'Plant input: {} at latitude: {}\tlongitude: {}'.format(str(pi_ann), round(lat, 3), round(lon, 3))
        lggr_info(mess)

    # only overwrite if arable
    # ========================
    ltd_data.plantInput = []
    for yr_str, land_luse in zip(pi_ann, ltd_data.landUses):
        if land_luse == 1 or land_luse == 2:
            ltd_data.plantInput.append(round(float(pi_ann[yr_str]),3))
        else:
            ltd_data.plantInput.append(0)
    return

def cnvrt_joe_plant_inputs_to_df(fname):
    """
    plant inputs are for 21 years - 252 values
    start year is 2020
    monthly values - convert to annual, construct data frame
    """

    print('Reading plant inputs file: ' + fname)

    # get header fields and start year
    # ================================
    with open(fname) as csv_fobj:
        csv_reader = reader(csv_fobj, delimiter = ' ')
        for header in csv_reader:

            # check the header list for compliance
            # ====================================
            for fld in REQUIRED_FIELDS:
                if fld not in header:
                    print('Field: ' + fld +  ' must be in header: ' + ', '.join(header[:6]))
                    return None

            # start year and number of years
            # ==============================
            strt_year = int(header[5].split('.')[0].lstrip('X'))
            nfields = len(header)
            nyears = int(nfields/12)

            nflds_lctn = nfields % 12
            if nflds_lctn != 5:
                print(WARN_STR + 'should be 5 location fields, got ' + str(nflds_lctn))
            print(f'Location fields: {", ".join(header[:5])}')

            break

    # convert to int
    # ==============
    data_frame = read_csv(fname, sep = ' ', names = header, skiprows=1, usecols = range(1,nfields + 1))
    for metric in list(['mu_global', 'gran_lat', 'gran_lon']):
        data_frame[metric] = data_frame[metric].astype(str).astype(int)

    print('Created data frame of length {} from {}'.format(len(data_frame), fname))

    return (strt_year, nyears, data_frame)
