rem 
@set setup_dir=E:\AbUniv\GlobalEcosseSuite\setup\
@set source_dir=G:\AbUnivGit\HoliSoilsSpGlEc\
@set envmdlng_dir=G:\AbUnivGit\EnvMdllngModuls\
@set PYTHONPATH=%envmdlng_dir%EnvModelModules;%source_dir%GlblEcosseModulesLtd

@set initial_working_dir=%cd%
@chdir /D %setup_dir%ltd_data
start cmd.exe /k "E:\Python38\python.exe -W ignore %source_dir%GlblEcssLtdSpVc\GlblEcsseHwsdGUI.py"
@chdir /D %initial_working_dir%
