"""
-------------------------------------------------------------------------------
 Name:        glbl_ecss_cmmn_funcs.py
 Purpose:     module comprising functions pertaining to bio-physical processes
 Author:      Mike Martin
 Created:     15/06/2023
 Licence:     <your licence>
-------------------------------------------------------------------------------
"""

__prog__ = 'glbl_ecss_bio_fns.py'
__version__ = '0.0.0'

# Version history
# ---------------
#

from math import exp

WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '

def add_npp_zaks_by_month(management, pettmp, soil_water, tstep, t_grow):
    """
    This differs from the  calculation presented by Zaks et al. (2007) in that the net primary production was
    calculated monthly using the water stress index for the previous month.
    """
    gdds_scle_factr = 11500
    iws_scle_factr  = 2720
    if management.pi_props[tstep] > 0:
        wat_strss_indx = soil_water.data['wat_strss_indx'][tstep]
        tgdd = pettmp['grow_dds'][tstep]
        npp = (0.0396/(1 + exp(6.33 - 1.5*(tgdd/gdds_scle_factr))))*(39.58*wat_strss_indx - 14.52)
        npp_month = iws_scle_factr * max(0, npp/t_grow)    # (eq.3.2.1)
    else:
        npp_month = 0

    management.npp_zaks[tstep] = npp_month

    return

def _miami_dyce_growing_season(precip, tair, land_cover_type='ara'):
    """
    tair     mean annual temperature
    precip   mean total annual precipitation

    from section 3.1 Net primary production from average temperature and total rainfall during the growing season
    modification of the well-established Miami model (Leith, 1972)

    units are tonnes
    """
    nppp = 15 * (1 - exp(-0.000664 * precip))       # precipitation-limited npp
    nppt = 15 / (1 + exp(1.315 - 0.119 * tair))     # temperature-limited npp
    npp = min(nppt, nppp)

    return npp

class MiamiDyceAnnualNpps(object):

    def __init__(self, climgen, pettmp):
        """
        return list of miami dyce npp estimates based on rainfall and temperature for perennial crop
        """
        ntsteps = len(pettmp['precipitation'])
        nyrs = int(ntsteps/12)

        npps = []
        for iyr in range(nyrs):
            indx = iyr*12
            precip_cumul = sum(pettmp['precipitation'][indx:indx + 12])
            tair_ave = sum(pettmp['temperature'][indx:indx + 12])/12
            npp = _miami_dyce_growing_season(precip_cumul, tair_ave)
            npps.append(npp)

        self.npps = npps
        self.end_yr  = climgen.fut_end_year
        self.strt_yr = climgen.fut_start_year

def generate_miami_dyce_grow_season_npp(pettmp, mngmnt):
    """
    return list of miami dyce npp estimates based on rainfall and temperature for growing months only
    """
    ntsteps = mngmnt.ntsteps
    precip_cumul = 0
    tair_cumul = 0
    tgrow = 0
    strt_indx = None
    for tstep in range(ntsteps):

        # second condition covers last month of last year
        # ===============================================
        if mngmnt.pi_props[tstep] == 0 or tstep == (ntsteps - 1):
            if tgrow > 0:

                # fetch npp for this growing season and backfill monthly NPPs
                # ===========================================================
                tair_ave = tair_cumul / tgrow
                npp = _miami_dyce_growing_season(precip_cumul, tair_ave)
                mngmnt.npp_miami_grow.append(npp)

                npp_mnthly = npp/tgrow
                for indx in range(strt_indx, tstep):
                    mngmnt.npp_miami[indx] = npp_mnthly

                tgrow = 0
                precip_cumul = 0
                tair_cumul = 0
                strt_indx = None
        else:
            if strt_indx is None:
                strt_indx = tstep
            precip_cumul += pettmp['precip'][tstep]
            tair_cumul += pettmp['tair'][tstep]
            tgrow += 1

    return
