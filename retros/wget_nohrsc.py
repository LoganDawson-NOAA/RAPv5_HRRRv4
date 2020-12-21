# Author: L Dawson
#
# Script to pull fcst files from HPSS, rename, and save in desired data directory

import numpy as np
import datetime, time, os, sys, subprocess
import dawsonpy


# Determine case name
try:
    case = str(sys.argv[1])
except IndexError:
    case = None

if case is None:
    case = raw_input('Enter a case name (do not included spaces): ')


# Get machine and head directory
machine, hostname = dawsonpy.get_machine()

if machine == 'WCOSS':
    DIR = '/gpfs/'+hostname[0]+'p2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'WCOSS_C':
    DIR = '/gpfs/hps3/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'WCOSS_DELL_P3':
    DIR = '/gpfs/dell2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'HERA':
    raise NameError, 'Need to modify wget_nohrsc.py script to define data/script head directory'


# Make sure data directory doesn't fail
try:
    OBS_DIR = os.environ['OBS_DIR']
except KeyError:
    OBS_DIR = DIR+'/obs'
    print('OBS_DIR not defined in environment. Setting OBS_DIR to '+OBS_DIR)

DATA_DIR = OBS_DIR+'/nohrsc'
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)
os.chdir(DATA_DIR)


# Make sure script directory doesn't fail
try:
    SCRIPT_DIR = os.environ['SCRIPT_DIR']
except KeyError:
    SCRIPT_DIR = DIR+'/'+case+'/runscripts'
    print('SCRIPT_DIR not defined in environment. Setting SCRIPT_DIR to '+SCRIPT_DIR)


# Create list of cycles
try:
    day1_00z_cycle = str(sys.argv[2])
except IndexError:
    day1_00z_cycle = None

if day1_00z_cycle is None:
    day1_00z_cycle = raw_input('Enter model initialization time (YYYYMMDDHH): ')

yyyy = int(day1_00z_cycle[0:4])
mm   = int(day1_00z_cycle[4:6])
dd   = int(day1_00z_cycle[6:8])
hh   = int(day1_00z_cycle[8:10])
date_str = datetime.datetime(yyyy,mm,dd,hh,0,0)

nhrs = np.arange(-12,61,6)
vtime_list = [date_str + datetime.timedelta(hours=x) for x in nhrs]


# Settings for NOHRSC
try:
    anl_str = str(sys.argv[3])
except IndexError:
    anl_str = None

if anl_str is None:
    anl_str = 'NOHRSC'


# Loop through each cycle
for this_vtime in vtime_list:

    vtime = this_vtime.strftime('%Y%m%d%H')

    # Only download missing analyses
    if not os.path.exists(DATA_DIR+'/sfav2_CONUS_6h_'+vtime+'_grid184.grb2'):
        print "Working on "+vtime+" NOHRSC analysis"
        os.system('wget --tries=2 https://www.nohrsc.noaa.gov/snowfall_v2/data/'+vtime[0:6]+'/sfav2_CONUS_6h_'+vtime+'_grid184.grb2')



os.system('touch '+SCRIPT_DIR+'/'+str.lower(anl_str)+'_done')
print("Done")
