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
    raise NameError, 'Need to modify htar_st4.py script to define data/script head directory'


# Make sure data directory doesn't fail
try:
    OBS_DIR = os.environ['OBS_DIR']
except KeyError:
    OBS_DIR = DIR+'/obs'
    print('OBS_DIR not defined in environment. Setting OBS_DIR to '+OBS_DIR)

DATA_DIR = OBS_DIR+'/st4'
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

nhrs = np.arange(-12,61,1)
vtime_list = [date_str + datetime.timedelta(hours=x) for x in nhrs]


# Settings for HRRRv3
anl_str = 'ST4'


# Loop through each cycle
for this_vtime in vtime_list:

    vtime = this_vtime.strftime('%Y%m%d%H')
    print "Working on "+vtime+" Stage IV analysis"

    RH_PREFIX = '/NCEPPROD/hpssprod/runhistory/rh'+vtime[0:4]+'/'+vtime[0:6]+'/'+vtime[0:8]
    TAR_PREFIX = 'com2_pcpanl_prod_pcpanl.'
    TAR_SUFFIX  = '.tar'
    FILE_PREFIX  = 'ST4.'
    FILE_SUFFIX  = '.01h.gz'


    htar_ball = RH_PREFIX+'/'+TAR_PREFIX+vtime[0:8]+TAR_SUFFIX
    htar_fname = './'+FILE_PREFIX+vtime+FILE_SUFFIX
    mv_from = './'+FILE_PREFIX+vtime[8:10]+FILE_SUFFIX
    mv_to = './'+FILE_PREFIX+vtime[8:10]+FILE_SUFFIX



    # Get rid of any existing file
    if int(vtime[8:])%24 == 0 or vtime == vtime_list[0].strftime('%Y%m%d%H'):
        req = {'htar_ball':'','htar_fname':[],'mv_from':[],'mv_to':[]}
        # Make temp text file to use for HPSS request
        os.system('rm -f ./'+str.lower(anl_str)+'_temp.txt')

    # Update dict
    req['htar_fname'].append(htar_fname)
    req['mv_from'].append(mv_from)
    req['mv_to'].append(mv_to)

    # Delete analysis file if path exists but file size is less than expected file size (trying to catch empty files)
    fname = DATA_DIR+'/'+FILE_PREFIX+vtime[8:10]+FILE_SUFFIX
    if os.path.exists(fname) and os.stat(fname).st_size < 700000:
        os.system('rm -f '+fname)

    # Only include missing analyses in file download listing
    if not os.path.exists(fname):
        o = open("./"+str.lower(anl_str)+"_temp.txt","w")
        o.write('\n'.join(req['htar_fname']))
        o.close()


    if vtime[8:] == '23' or vtime == vtime_list[-1].strftime('%Y%m%d%H'):
        req['htar_ball'] = htar_ball

        # Submit HPSS request
        os.system('htar -xvf '+req['htar_ball']+' -L '+str.lower(anl_str)+'_temp.txt')
      # print('htar -xvf '+req['htar_ball']+' -L '+str.lower(anl_str)+'_temp.txt')

        # Iterate through every item that was requested
      # for idx, (mv_from, mv_to) in enumerate(zip(req['mv_from'],req['mv_to'])):
      #     os.system("mv "+mv_from+" "+mv_to)
         #  print("mv "+mv_from+" "+mv_to)


os.system('gunzip ./*.gz')
os.system('rm -f ./'+str.lower(anl_str)+'_temp.txt')
os.system('touch '+SCRIPT_DIR+'/'+str.lower(anl_str)+'_done')
print("Done")
