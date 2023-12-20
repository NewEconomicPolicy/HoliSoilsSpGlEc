__author__ = 'soi698'


def boelter(bd, dz):
    """
    Estimates the water retention properties of peat soils using the
    pedotransfer functions of Boelter (1969).

    Boelter, D.H. (1969) Physical properties of peats as related to degree
    of decomposition. Soil sci. soc. amer. proc. vol 33, 1969

    Args
    ----
    bd: float
        Bulk density [g/cm3]
    dz: float
        Depth (thickness) of layer for which soil properties apply [mm].

    Returns
    -------
    wp: float
        Soil water content at wilting point [mm/dz]
    aw_fc: float
        Available water content at field capacity [mm/dz]
    aw_sat: float
        Available water content at saturation [mm/dz]
    """
    # Calculate wilt point of layer (based on Boelter's equation for 15 bar)
    wp = 1.57 + (115.28 * bd) - (107.77 * bd**2)
    wp = wp / 100.0 * dz  # Convert % vol to mm/lyr

    # Calc field capacity of layer (based on Boelter's equation for 0.1 bar)
    aw_fc = 2.06 + (719.35 * bd) - (1809.68 * bd**2)
    aw_fc = aw_fc / 100.0 * dz  # Convert % vol to mm/lyr

    # Calculate saturated water content of layer
    aw_sat = 99.0 - (123.45 * bd) + (252.92 * bd**2)
    aw_sat = aw_sat / 100.0 * dz  # Convert % vol to mm/lyr

    return wp, aw_fc, aw_sat


def bss(sand, silt, clay, oc, bd, depth, top):
    """
    For field soils, use equation of British Soil Service to find wilt point,
    field capacity and saturation water content
    from sand, silt, clay, carbon and bulk density data [H92]
    Equation seems to perform well [D04, G04]
    Updated values obtained from LEACHM February 2011 revision, provided
    by John Hutson [H11]

    eqn of the form: SWC_p = a + b*clay + c*silt + d*OC + e*BD
    where SWC_p is soil water capacity at pressure p, a is a constant
    b is coefficient for % clay, c for silt, d for organic carbon and e for bulk density
    change of form for equation in subsoil at high pressures, as noted below

    At pressures: [-1500., -200., -40., -10., -5.]kPa
    For topsoil:
    a = [0.0611,0.0938,0.2668,0.4030,0.4981]
    b = [0.004,0.0047,0.0039,0.0034,0.0027]
    c = [0.0005,0.0011,0.0013,0.0013,0.0011]
    d = [0.005,0.0069,0.0046,0.004,0.003]
    e = [0.,0.,-0.0764,-0.125,-0.1778]

    For subsoil*:
    a = [0.0125,0.0431,0.2205,0.3086,0.4216]
    b = [0.0092,0.0108,0.0047,0.0040,0.0034]
    c = [-0.000062,-0.000079,0.0020,0.0021,0.0018]
    d = [0.,0.,0.0093,0.0126,0.0022]
    e = [0.,0.,-0.0956,-0.1246,-0.1697]

    *in subsoil at presssures -1500 and -200kPa, equation form changes to:
    SWC_p = a + b*clay + c*clay^2

    We are interested in p = -1500kPa and 0kPa for wilt point and saturation
    and a range of pressures for field capacity [R44]:
    -1kPa (>5% OC), -10kPa (loamy), -100kPa (heavy clay) and -5kPa for typical UK field soil
    loamy and heavy clay properties defined by soil texture triangle [H12]

    To obtain coefficients for equation at small p, simply extend line
    Perform similar for coefficients between any other values.  Similar to procedure used by Donatelli et al. [D04]

    For field soils saturation (0kPa), rather than use PTF (which doesn't apply at very low pressure), find total pore space per unit volume
    The fraction f of volume occupied by solid matter is roughly (bulk density)/(particle density), hence saturation is given by 1-f
    Particle density is assumed 2.65g/cm3 for mineral fractions, and 1.5g/cm3 for organic fractions [P70,R73]

    References:
    [D04] M. Donatelli, J.H.M. Wosten, G. Belocchi (2004). Methods to evaluate pedotransfer functions
        Developmennt of Pedotransfer Functions in Soil Hydrology, Developments in Soil Science 30, 357-411
    [G04] Givi et al. (2004), Agricultural Water Management 70, 83-96
    [H92] Hutson et al. (1992), LEACHM manual (http://www.apesimulator.org/help/utilities/pedotransfers/British_Soil_Service_(topsoil).html)
    [H11] Hutson (2011), LEACHM manual update
    [H12] HWSD Documentation 2012
    [P70] Puustjarvi (1970), Peat Plant News 3,48-52
    [R44] Richards and Weaver (1944), J Ag Res 215-235
    [R73] Reid (1973), Area 5(1) 10-12

    Args
    ----
    sand, silt, clay, oc, bd, depth:
        values for the current layer [%, %, %, %, g/cm3, mm]
    top: logical
        True if topsoil, eslse False.

    Returns
    -------
    wp, fc, sa: wilt point, field capacity and saturation for layer (mm)
    """
    nullv = 999  # null value
    # np: number of pressures of interest for British Soil Service eqns
    # wilt point, field capacity, saturation [vol/vol, then converted to mm]
    wpfcsa = [None] * 3
    # pressures of interest (wilt point, all possible field capacities,
    # saturation)
    p = [None] * 5
    # coefficients a-e (rows) at all pressure of interest (columns) for layer
    coef = [[None] * 5 for i in range(5)]
    # coefficients at original defined pressures in layer
    coef0 = [[None] * 5 for i in range(5)]
    # coefficient and constant for straight lines between each pressure for
    # each coefficient in layer
    m = [[None] * 5 for i in range(5)]
    # coefficient and constant for straight lines between each pressure for
    # each coefficient in layer
    c = [[None] * 5 for i in range(5)]
    vi = [None] * 2  # index of coefficient vectors for wp and fc

    # constants ----------------------------------------------------------------
    # field soils
    np0 = 5  # number of input pressures
    nc0 = 5  # number of input coefficients
    p0 = [-1500., -200., -40., -10., -5.]  # original defined pressures
    if top:
        coef0[0] = [0.0611, 0.0938, 0.2668, 0.4030, 0.4981]
        coef0[1] = [0.004, 0.0047, 0.0039, 0.0034, 0.0027]
        coef0[2] = [0.0005, 0.0011, 0.0013, 0.0013, 0.0011]
        coef0[3] = [0.005, 0.0069, 0.0046, 0.004, 0.003]
        coef0[4] = [0., 0., -0.0764, -0.125, -0.1778]
    else:   # Subsoil coefficients
        coef0[0] = [0.0125, 0.0431, 0.2205, 0.3086, 0.4216]
        coef0[1] = [0.0092, 0.0108, 0.0047, 0.0040, 0.0034]
        coef0[2] = [-0.000062, -0.000079, 0.0020, 0.0021, 0.0018]
        coef0[3] = [0., 0., 0.0093, 0.0126, 0.0022]
        coef0[4] = [0., 0., -0.0956, -0.1246, -0.1697]

    # define pressures of interest (n.b. if you change these, you must also'
    # change the method to define vi, which must correspond to indices)
    np = 5
    #  wilt point pressure first, then possible field capacities
    p = [-1500., -100., -10., -5., -1.]
    vi[0] = 1   # i.e. wp always -1500kPa

    # Speed optimisations
    nc0_range = list(range(nc0))
    np0_range = list(range(np0))
    np0_minus_1_range = list(range(np0 - 1))
    # end constants ------------------------------------------------------------

    # find straight line between each pressure for each coefficient and layer
    for i in nc0_range:
        for j in np0_minus_1_range:
            m[i][j] = (coef0[i][j + 1] - coef0[i][j]) / (p0[j + 1] - p0[j])
            # -ve as it's really +(0-p0(j))
            c[i][j] = coef0[i][j] - p0[j] * m[i][j]

    # account for change of form for equations in subsoil at -200kPa and below
    if not top:
        # assume equation form changes at exactly -200kPa, hence line from
        # -200 to -40 is same as -40 to -10kPa
        m[2][1] = m[2][2]
        c[2][1] = c[2][2]

    # calculate coefficients for pressures of interest
    for i in range(np):
        for j in nc0_range:
            # find which straight line to use
            for k in np0_minus_1_range:
                if ((p[i] > p0[k] and p[i] <= p0[k + 1])
                        or (p[i] <= p0[0] and k == 0)
                        or (p[i] > p0[np0 - 1] and k == np0 - 2)):
                    # use line between points, or extend line at either end if
                    # smaller/larger than quoted range
                    coef[j][i] = m[j][k] * p[i] + c[j][k]

    # determine index for p at field capacity from soil properties
    vi[1] = 3  # default index for -5kPa
    if ((clay <= 25 and clay > 20 and silt <= 50 and silt >= 28)
            or (clay <= 20 and clay >= 8 and silt <= 50 and sand <= 52)):
        vi[1] = 2  # loamy -10kPa
    if clay >= 60:
        vi[1] = 1  # heavy clay -100kPa
    if oc > 5:
        vi[1] = 4   # organic soil -1kPa

    # calculate WP, FC, sat (eqns assume BD in g/cm3=kg/dm3, all others as %
    # (not fractions); results as fractions (vol/vol))
    for j in range(3):
        if j == 2:
            # saturation found from bulk density and estimated particle density
            wpfcsa[j] = 1 - bd / (2.65 * (1 - oc / 100.0) + 1.5 * oc / 100.0)
        else:
            # use British Soil Service PTF
            wpfcsa[j] = (
                coef[0][vi[j]] + coef[1][vi[j]] * clay +
                coef[2][vi[j]] * silt + coef[3][vi[j]] * oc +
                coef[4][vi[j]] * bd
            )
            if not top and p[vi[j]] <= -200:
                # low pressures in subsoil
                # replacing silt term by clay^2 term for -1500 to -200kPa in
                # subsoil
                wpfcsa[j] = wpfcsa[j] + coef[2][vi[j]] * ((clay**2) - silt)

    # convert from vol/vol to mm/ref_depth by multiplying by the layer depth
    for j in range(3):
        wpfcsa[j] = depth * wpfcsa[j]

    # set any nulls in either layer to zero
    for j in range(3):
        if sand == nullv:
            wpfcsa[j] = 0.0

    wp = wpfcsa[0]
    fc = wpfcsa[1]
    sa = wpfcsa[2]

    return wp, fc, sa
