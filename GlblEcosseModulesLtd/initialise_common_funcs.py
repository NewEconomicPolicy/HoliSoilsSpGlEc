"""
#-------------------------------------------------------------------------------
# Name:        initialise_common_funcs.py
# Purpose:     script of common initialisation functions for single site and spatial versions for
#                                                                                        Global Ecosse limited data mode
# Author:      Mike Martin
# Created:     31/07/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
"""
__prog__ = 'initialise_common_funcs.py'
__version__ = '0.0.0'

# Version history
# ---------------
#
from os.path import exists, normpath, split, splitext, isfile, isdir, join, lexists
from os import getcwd, remove, makedirs, mkdir, name as os_name
from json import load as json_load, dump as json_dump
from time import sleep
import sys
from glob import glob
from netCDF4 import Dataset, num2date
from PyQt5.QtWidgets import QApplication

from hwsd_mu_globals_fns import HWSD_mu_globals_csv
from glbl_ecss_cmmn_cmpntsGUI import print_resource_locations
from set_up_logging import set_up_logging
from weather_datasets import read_weather_dsets_detail
from hwsd_bil import check_hwsd_integrity

WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '
SETTINGS_LIST = ['config_dir', 'fname_png', 'hwsd_dir', 'log_dir', 'mask_fn', 'shp_dir', 'sims_dir', 'weather_dir']
BBOX_DEFAULT = [116.90045, 28.2294, 117.0, 29.0]  # bounding box default - somewhere in SE Europe
sleepTime = 5

# ===================================================

def initiation(form, variation=''):
    """
    this function is called to initiate the programme to process non-GUI settings.
    """
    if form.version == 'HWSD_grid':
        glbl_ecsse_str = 'global_ecosse_config_hwsd_'
        fname_setup = 'global_ecosse_setup_ver2' + variation + '.json'
    else:
        glbl_ecsse_str = 'global_ecosse_config_sngl_'
        fname_setup = 'global_ecosse_setup_ver2_sngl.json'

    # retrieve settings
    # =================
    _read_setup_file(form, fname_setup)

    form.glbl_ecsse_str = glbl_ecsse_str
    config_files = build_and_display_studies(form)

    # form.config_file = join(form.config_dir, 'global_ecosse_config.txt' ) # configuration file
    if len(config_files) > 0:
        form.config_file = config_files[0]
    else:
        form.config_file = form.config_dir + '/' + glbl_ecsse_str + 'dummy.txt'

    fname_model_switches = 'Model_Switches.dat'
    cwDir = getcwd()
    default_model_switches = join(cwDir, fname_model_switches)
    if isfile(default_model_switches):
        form.default_model_switches = default_model_switches
    else:
        print('{0} file does not exist in directory {1}'.format(fname_model_switches, cwDir))
        sleep(sleepTime)
        sys.exit(0)

    # set up logging
    # ==============
    form.settings['log_dir'] = form.log_dir
    set_up_logging(form, 'global_ecosse_min')

    if not check_sims_dir(form):
        sleep(sleepTime)
        sys.exit(0)

    # create dump files for grid point with mu_global 0
    form.fobjs = {}
    output_fnames = list(['nodata_muglobal_cells_v2b.csv'])
    if form.zeros_file:
        output_fnames.append('zero_muglobal_cells_v2b.csv')
    for file_name in output_fnames:
        long_fname = join(form.log_dir, file_name)
        key = file_name.split('_')[0]
        if exists(long_fname):
            try:
                remove(long_fname)
            except PermissionError as err:
                mess = 'Failed to delete mu global zeros dump file: {}\n\t{} '.format(long_fname, err)
                print(mess + '\n\t- check that there are no other instances of GlblEcosse')
                sleep(sleepTime)
                sys.exit(0)

        form.fobjs[key] = open(long_fname, 'w')

    # if specified then create pandas object read deom HWSD CSV file
    # ==============================================================
    if 'aoi_fname' in form.settings:
        form.hwsd_mu_globals = HWSD_mu_globals_csv(form, form.settings['aoi_fname'])
        # print('Reading AOI HWSD file ' + form.settings['aoi_fname'])

    return

def check_nc_start_end_dates(nc_fname, time_var_name='time', cru_future_flag=False):
    """
    ascertain the year span for historic datasets
    """
    dset_inp = Dataset(nc_fname)
    time_var = dset_inp.variables[time_var_name]
    if cru_future_flag:
        print('CRU future time units attribute: ' + time_var.units)
        start_year = int(time_var.units.split(' ')[-1])
        end_year = start_year + int(len(time_var) / 12) - 1  # year starts in January so knock a year off
    else:
        if 'calendar' in time_var.ncattrs():
            calendar_attr = time_var.calendar
        else:
            calendar_attr = 'standard'

        # Get the start and end date of the time series (as datetime objects):
        # ====================================================================
        start_day = time_var[0]
        start_date = num2date(start_day, units=time_var.units, calendar=calendar_attr)
        end_day = time_var[-1]
        end_date = num2date(end_day, units=time_var.units, calendar=calendar_attr)
        start_year = start_date.year
        end_year = end_date.year
    dset_inp.close()

    return start_year, end_year

def _read_setup_file(form, fname_setup):
    """
    read settings used for programme from the setup file, if it exists,
    or create setup file using default values if file does not exist
    """
    func_name = __prog__ + ' _read_setup_file'

    setup_file = join(getcwd(), fname_setup)
    if exists(setup_file):
        try:
            with open(setup_file, 'r') as fsetup:
                settings = json_load(fsetup)
        except (OSError, IOError) as e:
            sleep(sleepTime)
            exit(0)
    else:
        settings = write_default_setup_file(setup_file)
        print('Read default setup file ' + setup_file)

    # initialise vars
    # ===============
    form.settings = {}
    form.images_dir = ''
    form.lu_pi_content = {}  # TODO

    # validate setup file
    # ===================
    grp = 'setup'
    for key in SETTINGS_LIST:
        if key not in settings[grp]:
            print(ERROR_STR + 'setting {} is required in setup file {} '.format(key, setup_file))
            sleep(sleepTime)
            exit(0)

    config_dir = settings[grp]['config_dir']
    form.fname_png = settings[grp]['fname_png']
    hwsd_dir = settings[grp]['hwsd_dir']
    log_dir = settings[grp]['log_dir']
    mask_fn = settings[grp]['mask_fn']
    form.shp_dir = settings[grp]['shp_dir']
    sims_dir = settings[grp]['sims_dir']
    weather_dir = settings[grp]['weather_dir']

    # applies to Global Ecosse SpVc version only
    # ==========================================
    if 'aoi_fname' in settings[grp]:
        aoi_fname = settings[grp]['aoi_fname']
        if isfile(aoi_fname):
            form.settings['aoi_fname'] = aoi_fname
        else:
            print(ERROR_STR + 'AOI CSV file {} must exist' + aoi_fname)
            sleep(sleepTime)
            exit(0)

    # ==============
    if isfile(mask_fn):
        form.mask_fn = mask_fn
    else:
        if mask_fn != '':
            print(WARN_STR + 'Land use mask file ' + mask_fn + ' does not exist')
        form.mask_fn = None

    # additional settings to enable ECOSSE to be run
    # ==============================================
    ecosse_run_flag = False
    runsites_py = None
    python_exe = None
    ecosse_exe = None

    runsites_config_file = join(config_dir, 'global_ecosse_ltd_data_runsites_config.txt')
    if isfile(runsites_config_file):
        if 'python_exe' in settings[grp].keys() and 'runsites_py' in settings[grp].keys():

            runsites_py = settings[grp]['runsites_py']
            python_exe = settings[grp]['python_exe']

            if isfile(runsites_py) and isfile(python_exe):
                form.runsites_py = runsites_py
                form.python_exe = python_exe

                # ascertain which version of Ecosse is defined in the runsites file
                # =================================================================
                if type(runsites_config_file) is str:
                    with open(runsites_config_file, 'r') as fconfig:
                        config = json_load(fconfig)

                    try:
                        fn = split(config['Simulations']['exepath'])[1]
                    except KeyError as err:
                        print(
                            WARN_STR + 'could not identify Ecosse version in run sites config file: '
                            + runsites_config_file + '\n')
                    else:
                        ecosse_exe = splitext(fn)[0].lower()
                        ecosse_run_flag = True
            else:
                if not isfile(runsites_py):
                    print(WARN_STR + 'Run sites script file: ' + runsites_py + ' does not exist')
                if not isfile(python_exe):
                    print(WARN_STR + 'Python interpreter: ' + python_exe + ' does not exist')
        else:
            print('\n' + WARN_STR + 'python_exe and runsites_py should be in ' + setup_file)
    else:
        print(WARN_STR + 'Run sites configuration file: ' + runsites_config_file + ' does not exist')
        runsites_config_file = None

    if ecosse_run_flag:
        print('\nRun ECOSSE settings:')
        print('\tRun sites config file: ' + runsites_config_file)
        print('\tRun sites script file: ' + runsites_py)
        print('\tPython interpreter:    ' + python_exe)
        print('')
    else:
        print('\t\tcannot run ECOSSE\n')

    form.runsites_config_file = runsites_config_file
    form.runsites_py = runsites_py
    form.python_exe = python_exe
    form.ecosse_exe = ecosse_exe
    form.ecosse_run_flag = ecosse_run_flag

    # ============= End of runsites config section ====================

    # check for AMMA-2050: African Monsoon Multidisciplinary Analysis 2050
    # ===================
    allowed_gcms = []
    amma_2050_dir = ''
    gcm_list = []
    if 'amma_2050_dir' in settings[grp].keys() and 'gcm_list' in settings[grp].keys():
        amma_2050_dir = settings[grp]['amma_2050_dir']
        if isdir(amma_2050_dir):
            gcm_list = settings[grp]['gcm_list']
            for fname in glob(amma_2050_dir + '\\*'):
                if isdir(fname):
                    dummy, gcm = split(fname)
                    allowed_gcms.append(gcm)
        else:
            print(ERROR_STR + 'path in setup file for AMMA_2050 data {} does not exist'.format(amma_2050_dir))
            amma_2050_dir = ''

    form.amma_2050_allowed_gcms = allowed_gcms
    form.amma_2050_dir = amma_2050_dir
    form.gcm_list = gcm_list

    # make sure directories exist for configuration and log files
    # ===========================================================
    if not lexists(log_dir):
        makedirs(log_dir)
    form.log_dir = log_dir

    if not lexists(config_dir):
        makedirs(config_dir)
    form.config_dir = config_dir

    # HWSD is crucial
    # ===============
    if lexists(hwsd_dir):
        check_hwsd_integrity(settings[grp]['hwsd_dir'])
        form.hwsd_dir = hwsd_dir
    else:
        print('Error reading {}\tHWSD directory {} must exist'.format(setup_file, hwsd_dir))
        sleep(sleepTime)
        exit(0)

    # weather is crucial
    # ===================
    if lexists(weather_dir):
        form.weather_dir = weather_dir
        form.settings['weather_dir'] = weather_dir
    else:
        print('Error reading {}\tClimate directory {} must exist'.format(setup_file, weather_dir))
        sleep(sleepTime)
        exit(0)

    if lexists(weather_dir):
        form.weather_dir = weather_dir
        form.settings['weather_dir'] = weather_dir

        # sims dir checked later
    # ======================
    form.sims_dir = sims_dir

    # check weather data
    # ==================
    if form.version == 'HWSD_grid':
        rqurd_wthr_rsrcs = ['CRU', 'CHESS']
        rqurd_wthr_rsrcs = ['CRU', 'EObs']  # required weather resources
    else:
        rqurd_wthr_rsrcs = ['CRU', 'EObs', 'HARMONIE']

    form.wthr_settings_prev = {}
    read_weather_dsets_detail(form, rqurd_wthr_rsrcs)

    # TODO: most of these are not used
    # ================================
    grp = 'run_settings'
    try:
        form.completed_max = settings[grp]['completed_max']
        form.start_at_band = settings[grp]['start_at_band']
        form.space_remaining_limit = settings[grp]['space_remaining_limit']
        form.kml_flag = settings[grp]['kml_flag']
        form.soilTestFlag = settings[grp]['soil_test_flag']
        form.zeros_file = settings[grp]['zeros_file']
    except KeyError:
        print('{}\tError in group: {}'.format(func_name, grp))
        sleep(sleepTime)
        exit(0)

    # report settings
    # ===============
    lta_nc_fname = None
    print_resource_locations(setup_file, config_dir, hwsd_dir, weather_dir, lta_nc_fname, sims_dir, log_dir)

    return True

def write_default_setup_file(setup_file):
    """
    stanza if setup_file needs to be created
    """
    # Windows only for now
    # =====================
    os_system = os_name
    if os_system != 'nt':
        print('Operating system is ' + os_system + 'should be nt - cannot proceed with writing default setup file')
        sleep(sleepTime)
        sys.exit(0)

    # return list of drives
    # =====================
    import win32api

    drives = win32api.GetLogicalDriveStrings()
    drives = drives.split('\000')[:-1]
    if 'S:\\' in drives:
        root_dir_app = 'S:\\tools\\'  # Global Ecosse installed here
        root_dir_user = 'H:\\'  # user files reside here
    else:
        root_dir_app = 'E:\\'
        root_dir_user = 'C:\\AbUniv\\'

    suite_path = root_dir_app + 'GlobalEcosseSuite\\'
    data_path = root_dir_app + 'GlobalEcosseData\\'
    outputs_path = root_dir_app + 'GlobalEcosseOutputs\\'
    root_dir_user += 'GlobalEcosseSuite\\'

    _default_setup = {
        'setup': {
            'config_dir': root_dir_user + 'config',
            'fname_png': join(suite_path + 'Images', 'Tree_of_life.PNG'),
            'hwsd_dir': data_path + 'HWSD_NEW',
            'images_dir': outputs_path + 'images',
            'log_dir': root_dir_user + 'logs',
            'shp_dir': data_path + 'CountryShapefiles',
            'sims_dir': outputs_path + 'EcosseSims',
            'weather_dir': data_path
        },
        'run_settings': {
            'completed_max': 5000000000,
            'start_at_band': 0,
            'space_remaining_limit': 1270,
            'kml_flag': True,
            'soil_test_flag': False,
            'zeros_file': False
        }
    }
    # create setup file
    # =================
    with open(setup_file, 'w') as fsetup:
        json_dump(_default_setup, fsetup, indent=2, sort_keys=True)
        fsetup.close()
        return _default_setup

def write_default_config_file(config_file):
    """
    #        ll_lon,    ll_lat  ur_lon,ur_lat
    # stanza if config_file needs to be created
    """
    _default_config = {
        'minGUI': {
            'aveWthrFlag': False,
            'bbox': BBOX_DEFAULT,
            'cordexFlag': 0,
            'hwsdCsvFname': '',
            'luPiJsonFname': '',
            'snglPntFlag': True,
            'usePolyFlag': False
        },
        'cmnGUI': {
            'climScnr': 'rcp26',
            'eqilMode': '9.5',
            'futStrtYr': '2006',
            'futEndYr': '2015',
            'gridResol': 0,
            'histStrtYr': '1980',
            'histEndYr': '2005',
            'study': ''
        }
    }
    # if config file does not exist then create it...
    with open(config_file, 'w') as fconfig:
        json_dump(_default_config, fconfig, indent=2, sort_keys=True)
        fconfig.close()
        return _default_config


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
    except (OSError, IOError) as e:
        print(e)
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
                transition_lu = land_use  # TODO: not sure what this does
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
            if exists(yield_map_fname):
                print("Yields map " + yield_map_fname + " exists")
            else:
                print("Yields map " + yield_map_fname + " does not exist")
                yield_map_fname = None
    else:
        print('Key {} not present in {} - will assume no Yield map and will set plant C input to zero for each year'
              .format(landuse_pi_key, lu_pi_json_fname))
    form.yield_map_fname = yield_map_fname

    if fileOkFlag:
        mess = land_use + 'landuse and plant input file is valid'
        form.w_create_files.setEnabled(True)

    return mess


def write_runsites_config_file(form):
    """
    read the runsites config file and edit one line
    """
    runsites_config_file = form.runsites_config_file
    try:
        with open(runsites_config_file, 'r') as fconfig:
            config = json_load(fconfig)
            print('Read config file ' + runsites_config_file)
    except (OSError, IOError) as err:
        print(err)
        return False

    # overwrite config file
    # =====================
    if hasattr(form, 'w_study'):
        sims_dir = normpath(join(form.sims_dir, form.w_study.text()))
    else:
        sims_dir = normpath(join(form.sims_dir, form.study))

    config['Simulations']['sims_dir'] = sims_dir
    with open(runsites_config_file, 'w') as fconfig:
        json_dump(config, fconfig, indent=2, sort_keys=True)
        print('Edited ' + runsites_config_file + '\n\twith simulation location: ' + sims_dir)
        QApplication.processEvents()

    return True

def build_and_display_studies(form):
    """
    is called at start up and when user creates a new project
    """

    glbl_ecsse_str = form.glbl_ecsse_str
    config_files = glob(form.config_dir + '/' + glbl_ecsse_str + '*.txt')
    studies = []
    for fname in config_files:
        dummy, remainder = fname.split(glbl_ecsse_str)
        study, dummy = splitext(remainder)
        if study != '':
            studies.append(study)
    form.studies = studies

    # display studies list
    # ====================
    if hasattr(form, 'combo00s'):
        form.combo00s.clear()
        for study in studies:
            form.combo00s.addItem(study)

    return config_files


def check_sims_dir(form):
    """
    called from _initiation func.
    """
    retFlag = False
    sims_dir = form.sims_dir
    form.sims_dir = ''  # default in case of failure

    # make sure directory has write permissions and it exists
    if not lexists(sims_dir):
        mkdir(sims_dir)
        form.sims_dir = sims_dir
        retFlag = True
    else:
        if isdir(sims_dir):
            form.lgr.info('Directory {0} already exists'.format(sims_dir))
            form.sims_dir = sims_dir
            retFlag = True
        elif isfile(sims_dir):
            print('{0} already exists but as a file' + sims_dir)
        else:
            print('{0} is not a directory or a file' + sims_dir)

    return retFlag
