"""
#-------------------------------------------------------------------------------
# Name:        glbl_ecss_cmmn_funcs.py
# Purpose:     module comprising miscellaneous functions common to some Global Ecosse variations
# Author:      Mike Martin
# Created:     31/07/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
"""

__prog__ = 'glbl_ecss_cmmn_funcs.py'
__version__ = '0.0.0'

# Version history
# ---------------
#
from os import remove, mkdir
from os.path import isdir, join, normpath, splitext, exists, lexists, isfile, split

from json import dump as json_dump, load as json_load
from json import JSONDecodeError

from pandas import read_excel

from time import sleep
from glob import glob

sleepTime = 5

NGRANULARITY = 120

WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '
DUMMY_STR = 'dummy'
SET_SPACER_LEN= 12
SUPER_G_RUN_CNFGS = ('all', 'nsyn_grz', 'nsyn_mnr', 'grz_mnr', 'base_line')

def formatSize(bytes):
    """
    derived from https://www.tutorialexample.com/a-simple-guide-to-python-get-disk-or-directory-total-space-used-space-and-free-space-python-tutorial
    """
    try:
        bytes = float(bytes)
        kb = bytes / 1024
    except:
        return "Error"
    if kb >= 1024:
        M = kb / 1024
        if M >= 1024:
            G = M / 1024
            return "%.2fG" % (G)
        else:
            return "%.2fM" % (M)
    else:
        return "%.2fkb" % (kb)

def check_runsites(setup_file, settings, config_dir, runsites_config_fn, python_exe, runsites_py):
    """

    """
    run_ecosse_flag = True

    # file to run Ecosse
    # ==================
    if isfile(runsites_py):
        print('Will use ECOSSE run script : ' + runsites_py)
    else:
        print(WARN_STR + '\tECOSSE run script ' + runsites_py + ' does not exist')
        run_ecosse_flag = False

    # check that the runsites configuration file exists
    # =================================================
    runsites_cnfg_fn = join(config_dir, runsites_config_fn)

    mess = 'Run sites configuration file ' + runsites_cnfg_fn
    if exists(runsites_cnfg_fn):
        mess += ' exists'
    else:
        print(WARN_STR + '\tRun sites configuration file ' + runsites_cnfg_fn + ' does not exist')
        runsites_cnfg_fn = None
        run_ecosse_flag = False

    print(mess)

    # make sure all components of the command string exist
    # ====================================================

    if not exists(python_exe):
        print(WARN_STR + '\tPython programme ' + python_exe + ' does not exist')
        run_ecosse_flag = False

    if not run_ecosse_flag:
        print(WARN_STR + 'cannot run ECOSSE - check setup file: ' + setup_file)

    settings['runsites_py'] = runsites_py
    settings['run_ecosse_flag'] = run_ecosse_flag
    settings['runsites_cnfg_fn'] = runsites_cnfg_fn

    return

def build_and_display_studies(form, glbl_ecsse_str):
    """
    is called at start up and when user creates a new project
    """
    func_name =  __prog__ + ' build_and_display_studies'

    if hasattr(form, 'sttngs'):
        config_dir = form.sttngs['config_dir']
    elif hasattr(form, 'config_dir'):
        config_dir = form.config_dir
    else:
        mess = ERROR_STR + 'could not display studies in function: ' + func_name
        print(mess + ' - form object must have either sttngs or config_dir attribute')
        sleep(sleepTime)
        exit(0)

    config_files = glob(config_dir + '/' + glbl_ecsse_str + '*.json')
    studies = []
    for fname in config_files:
        dummy, remainder = fname.split(glbl_ecsse_str)
        study, dummy = splitext(remainder)
        if study != '':
            studies.append(study)

    # filter out dummy study
    # ======================
    for study in studies:
        if DUMMY_STR in study:
            for cnfg_fn in config_files:
                if cnfg_fn.rfind('_' + DUMMY_STR + '\.') > 0:
                    studies.remove(study)
                    config_files.remove(cnfg_fn)
                    print(WARN_STR + 'removed ' + DUMMY_STR + ' config file: ' + cnfg_fn)

    form.studies = studies

    # display studies list
    # ====================
    if hasattr(form, 'combo00s'):
        form.combo00s.clear()
        for study in studies:
            form.combo00s.addItem(study)

    return config_files

def check_lu_pi_json_fname(form):
    """
    Land use and plant input file should look like this:
    {
      "YieldMap" : "",
      "LandusePI": {
        "0": ["Arable", 3995.0],
        "1": ["Grassland", 3990.0],
        "5": ["Arable", 3210.0],
        "17": ["Forestry",3203.0]
      }
    }
    """
    form.w_create_files.setEnabled(False)
    lu_pi_json_fname = form.w_lbl13.text()
    form.lu_pi_content = {}

    if not exists(lu_pi_json_fname):
        return 'Land use and plant input file does not exist'

    try:
        with open(lu_pi_json_fname, 'r') as flu_pi:
            lu_pi_content = json_load(flu_pi)
            flu_pi.close()

            print('Read land use and plant input file ' + lu_pi_json_fname)
            form.lu_pi_content = lu_pi_content
    except (OSError, IOError) as err:
            print(err)
            return 'Could not read land use and plant input file'

    # check file contents
    # ===================
    fileOkFlag = True
    landuse_pi_key = 'LandusePI'
    transition_lu = None
    if landuse_pi_key in lu_pi_content.keys():
        for key in lu_pi_content[landuse_pi_key]:
            land_use = lu_pi_content[landuse_pi_key][key][0]
            if land_use not in form.land_use_types:
                mess = 'Invalid landuse: ' + land_use
                long_mess = mess + ' in {} - must be any of these:  {}' \
                                                    .format(lu_pi_json_fname, str(form.land_use_types.keys()))
                print(long_mess)
                fileOkFlag = False
            elif key == '1':
                transition_lu = land_use    # TODO: not sure what this does
    else:
        mess = 'Key {} must be present in  {}'.format(landuse_pi_key, lu_pi_json_fname)
        fileOkFlag = False
    form.transition_lu = transition_lu

    # check yield map
    # ===============
    yield_map_fname = None
    yield_map_key = 'YieldMap'
    if yield_map_key in lu_pi_content.keys():
        yield_map_fname = lu_pi_content[yield_map_key]
        if yield_map_fname is not None:
            mess = "\tYields map "
            if yield_map_fname == '':
                mess += 'not specified'
                yield_map_fname = None

            elif exists(yield_map_fname):
                mess += yield_map_fname + ' exists'
            else:
                mess += yield_map_fname + ' does not exist'
                yield_map_fname = None

            print(mess)
    else:
        print('Key {} not present in {} - will assume no Yield map and will set plant C input to zero for each year'
                                                                        .format(landuse_pi_key, lu_pi_json_fname))
    form.yield_map_fname = yield_map_fname

    if fileOkFlag:
        mess = land_use + ' landuse and plant input file is valid'
        form.w_create_files.setEnabled(True)

    return mess

def write_study_definition_file(form, glbl_ecss_variation = None):
    """

    """
    # calculate_grid_cell(form)   # pulls resolution from GUI otherwise value is Null

    # prepare the bounding box
    # ========================
    try:
        ll_lon = float(form.w_ll_lon.text())
        ll_lat = float(form.w_ll_lat.text())
        ur_lon = float(form.w_ur_lon.text())
        ur_lat = float(form.w_ur_lat.text())
    except ValueError:
        ll_lon = 0.0
        ll_lat = 0.0
        ur_lon = 0.0
        ur_lat = 0.0
    bbox =  list([ll_lon,ll_lat,ur_lon,ur_lat])
    study = form.w_study.text()
    if study == '':
        study = 'xxxx'

    wthr_rsrce = form.combo10w.currentText()
    if wthr_rsrce == 'CRU':
        fut_clim_scen = form.combo10.currentText()
    else:
        fut_clim_scen = wthr_rsrce

    # construct land_use change - not elegant but adequate
    # =========================
    land_use = ''
    if hasattr(form, 'lu_pi_content'):
        for indx in form.lu_pi_content['LandusePI']:
            lu, pi = form.lu_pi_content['LandusePI'][indx]
            land_use += form.lu_type_abbrevs[lu] + '2'
    else:
        land_use = 'unk2unk'

    land_use = land_use.rstrip('2')

    # TODO: replace "luCsvFname": "" with 'luPiJsonFname': form.fertiliser_fname
    if hasattr(form, 'combo12'):
        crop_name = form.combo12.currentText()
    else:
        crop_name = 'Unknown'

    if hasattr(form, 'w_daily'):
        daily_mode = form.w_daily.isChecked()
    else:
        daily_mode = False

    if hasattr(form, 'w_lbl06'):
        hwsd_csv_fname = form.w_lbl06.text()
    else:
        hwsd_csv_fname = None

    if hasattr(form, 'sttngs'):
        sims_dir = form.sttngs['sims_dir']
    elif hasattr(form, 'settings'):
        sims_dir = form.settings['sims_dir']
    else:
        sims_dir = form.sims_dir

    study_defn = {
        'studyDefn': {
            'dailyMode': daily_mode,
            'bbox'     : bbox,
            "luCsvFname": "",
            'hwsdCsvFname' : hwsd_csv_fname,
            'study'    : study,
            'land_use' : land_use,
            'histStrtYr': form.combo09s.currentText(),
            'histEndYr' : form.combo09e.currentText(),
            'climScnr' : fut_clim_scen,
            'futStrtYr': form.combo11s.currentText(),
            'futEndYr' : form.combo11e.currentText(),
            'cropName' : crop_name,
            'province' : 'xxxx',
            'resolution' : form.req_resol_deg,
            'shpe_file': 'xxxx',
            'version'  : form.version
            }
        }

    # copy to sims area
    # =================
    if glbl_ecss_variation is None:
        study_defn_file = join(sims_dir, study + '_study_definition.txt')
        with open(study_defn_file, 'w') as fstudy:
            json_dump(study_defn, fstudy, indent=2, sort_keys=True)
            print('\nWrote study definition file ' + study_defn_file)

    elif glbl_ecss_variation == 'SuperG_MA':
        for sim_nm in SUPER_G_RUN_CNFGS:
            study_run = study + '_' + sim_nm
            study_defn_file = join(sims_dir, study_run + '_study_definition.txt')
            study_defn['studyDefn']['study'] = study_run
            with open(study_defn_file, 'w') as fstudy:
                json_dump(study_defn, fstudy, indent=2, sort_keys=True)
                print('Wrote study definition file ' + study_defn_file)

    return

def write_runsites_cnfg_file(form, identifer = None):

    func_name =  __prog__ +  ' write_runsites_cnfg_fn'

    # read the runsites config file and edit one line
    # ======================================
    runsites_cnfg_fn = form.sttngs['runsites_cnfg_fn']
    try:
        with open(runsites_cnfg_fn, 'r') as fconfig:
            config = json_load(fconfig)
            print('Read config file ' + runsites_cnfg_fn)
    except (JSONDecodeError, OSError, IOError) as err:
            print(err)
            return False

    # overwrite config file
    # =====================
    study = form.w_study.text()
    if identifer is not None:
        study += '_' + identifer

    sims_dir = normpath(join(form.sttngs['sims_dir'], study))
    if not isdir(sims_dir):
        mkdir(sims_dir)
        print('Created simulations directory ' + sims_dir)

    config['Simulations']['sims_dir'] = sims_dir
    config['General']['cropName'] = form.combo12.currentText()

    with open(runsites_cnfg_fn, 'w') as fconfig:
        json_dump(config, fconfig, indent=2, sort_keys=True)
        print('Edited ' + runsites_cnfg_fn + '\n...with simulation location: ' + sims_dir)
        return True

def read_manure_types(manure_types_fn):
    """
    read manure types and their corresponding codes
    """
    if isfile(manure_types_fn):
        data = read_excel(manure_types_fn)  # does not support ODS file
        manure_types = data.to_dict('list')
    else:
        manure_types = None

    return manure_types

def read_crop_pars_codes(crop_pars_fname, required_crops, nset_lines = 9):
    """
    read crop names and their corresponding codes
    """
    with open(crop_pars_fname) as fhand:
        lines = fhand.readlines()
    nlines = len(lines)

    crop_codes = {}
    for iline in range(0, nlines, nset_lines):
        crop_name = lines[iline].strip()
        nline_nxt = iline + 1
        if nline_nxt >= nlines:
            break

        try:
            code = int(lines[nline_nxt].strip())
        except ValueError:
            print(ERROR_STR + 'reading file : ' + crop_pars_fname)
            break

        if crop_name in required_crops:
            crop_codes[crop_name] = code

    return crop_codes

def check_sims_dir(lggr, sims_dir):
    """
    called from _initiation func.
    """
    ret_flag = False    

    # make sure directory has write permissions and it exists
    if not lexists(sims_dir):
        mkdir(sims_dir)
        ret_flag = True
    else:
        if isdir(sims_dir):
            lggr.info('Directory {} already exists'.format(sims_dir))
            ret_flag = True

        else:
            ret_flag = False
            if isfile(sims_dir):
                print('{} already exists but as a file' + sims_dir)
            else:
                print('{} is not a directory or a file' + sims_dir)

    return ret_flag, sims_dir

def create_dump_files(form):
    """
    called from _initiation func.
    create dump files for grid point with mu_global 0
    """
    study = form.w_study.text()

    form.fobjs = {}
    output_fnames = list(['nodata_muglobal_cells_','grid_cell_info_', 'nitro_cell_info_'])
    if form.zeros_file:
        output_fnames.append('zero_muglobal_cells_')

    for fn in output_fnames:
        fname = fn + study + '.csv'
        long_fname = join(form.settings['log_dir'], fname)
        key = fname.split('_')[0]
        if exists(long_fname):
            try:
                remove(long_fname)
            except (PermissionError) as err:
                mess = 'Failed to delete dump file: {}\n\t{} '.format(long_fname, err)
                print(mess + '\n\t- check that there are no other instances of GlblEcosse'.format(long_fname, err))
                sleep(sleepTime)
                exit(0)

        form.fobjs[key] = open(long_fname,'w')

    return

def write_kml_file(sim_dir, fname_short, mu_global, latitude, longitude):
    """
    write kml consisting of mu_global and soil details
    """
    # set Icon size
    scaleVal = 0.5

    output = 15 * ['']

    output[0] = '<?xml version="1.0" encoding="utf-8"?>\n'
    output[1] = '<kml xmlns="http://www.opengis.net/kml/2.2">\n'
    output[2] = '\t<Document>\n'
    output[3] = '\t<Style id="defaultStyle">\n'
    output[4] = '\t\t<LabelStyle><scale>' + str(scaleVal) + '</scale></LabelStyle>\n'
    output[5] = '\t\t<IconStyle id="hoverIcon"><scale>' + str(scaleVal) + '</scale></IconStyle>\n'
    output[6] = '\t</Style>\n'

    output[7] = '\t\t<Placemark>\n'
    output[8] = '\t\t\t<name>' + str(mu_global) + '</name>\n'
    output[9] = '\t\t\t<description>\n' + 'Long: ' + str(round(longitude,4)) + \
            ',   Lat:' + str(round(latitude,4)) +'\nMu_global:' + str(mu_global) + '\n\t\t\t</description>\n'

    output[10] = '\t\t\t<Point><coordinates>' + str(longitude) + ',' + str(latitude) + '</coordinates></Point>\n'

    output[11] = '\t\t\t<styleUrl>#defaultStyle</styleUrl>\n'
    output[12] = '\t\t</Placemark>\n'
    output[13] = '\t</Document>\n'
    output[14] = '</kml>\n'

    kml_fname = join(sim_dir, fname_short + '.kml')
    with open(kml_fname, 'w') as fpout:    # context manager
        fpout.writelines(output)

    return

def input_txt_line_layout(data, comment):
    """

    """
    spacer_len = max(SET_SPACER_LEN - len(data), 2)
    spacer = ' ' * spacer_len
    return '{}{}# {}\n'.format(data, spacer, comment)

def write_signature_file(sim_dir, mu_global, soil, latitude, longitude, province = ''):
    """
    write json consisting of mu_global and soil details
    """

    config = {
        'location': {
            'province' : province,
            'longitude' : round(longitude,6),
            'latitude'  : round(latitude,6),
            'share' : soil[-1]
        },
        'soil_lyr1': {
            'C_content': soil[0],
            'Bulk_dens': soil[1],
            'pH'       : soil[2],
            '%_clay': soil[3],
            '%_silt': soil[4],
            '%_sand': soil[5]
        }
    }
    # add subsoil layer, if it exists
    # ===============================
    if len(soil) >= 12:
       config['soil_lyr2'] =  {
                'C_content': soil[6],
                'Bulk_dens': soil[7],
                'pH'       : soil[8],
                '%_clay': soil[9],
                '%_silt': soil[10],
                '%_sand': soil[11]
            }

    signature_fname = join(sim_dir, str(mu_global) + '.txt')
    with open(signature_fname, 'w') as fsig:
        json_dump(config, fsig, indent=2, sort_keys=True)

    return

def write_manifest_file(study, fut_clim_scen, sim_dir, soil_list, mu_global, latitude, longitude, area):
    """
    write json consisting of mu_global and soil shares for each grid cell
    """

    # location etc.
    # =============
    manifest = {
        'location': {
            'longitude' : round(longitude,6),
            'latitude'  : round(latitude,6),
            'area' : round(area,8),
            'area_description' : study,
            'province' : 'province',
            'scenario' : fut_clim_scen
        }
    }
    # soil shares
    # ===========
    smu_global = str(mu_global)
    manifest[smu_global] =  {}
    for soil_num, soil in enumerate(soil_list):
        manifest[smu_global][soil_num + 1] =  soil[-1]

    # deprecated
    # ==========
    manifest['longitudes'] =  {}
    manifest['granular_longs'] = {}

    # construct file name and write
    # =============================
    manif_dir, fname_part2 = split(sim_dir)
    manifest_fname = join(manif_dir, 'manifest_' + fname_part2[:-4] + '.txt')
    with open(manifest_fname, 'w') as fmanif:
        json_dump(manifest, fmanif, indent=2, sort_keys=True)

    return

def write_study_manifest_files(form, lon_lats):
    """
    create study manifest file <study>_manifest.csv and summary <study>_summary_manifest.csv
    write first line and leave file object open for subsequent records
    """
    func_name =  __prog__ + ' write_study_manifest_files'
    NoData = -999

    if hasattr(form, 'sttngs'):
        sims_dir = form.sttngs['sims_dir']
    elif hasattr(form, 'sims_dir'):
        sims_dir = form.sims_dir
    else:
        mess = ERROR_STR + 'could not determine value for sims_dir in function: ' + func_name
        print(mess + ' - form object must have either sttngs or sims_dir attribute')
        return

    # set up mainfest file for the study
    # ==================================
    form.fstudy = []
    for frag_name in list(['_','_summary_']):
        study_fname = join(sims_dir, form.study + frag_name + 'manifest.csv')
        if exists(study_fname):
            remove(study_fname)
        form.fstudy.append(open(study_fname,'w',100))

    # write first two lines so they are distinguisable for sorting purposes
    for lon_lat_pair in lon_lats:
        lon, lat = lon_lat_pair
        gran_lon = round((180.0 + lon)*NGRANULARITY)
        gran_lat = round((90.0 - lat)*NGRANULARITY)
        form.fstudy[1].write('{}\t{}\t{}\t{}\t'.format(gran_lat, gran_lon, lat, lon))
        form.fstudy[0].write('{}\t{}\t{}\t{}\t{}\n'.format(gran_lat, gran_lon, lat, lon, NoData))

    form.fstudy[1].write('{}\t{}\n'.format(form.req_resol_deg, form.req_resol_granul))

    return