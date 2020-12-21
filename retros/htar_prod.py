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


# Create list of cycles
try:
    cycle = str(sys.argv[2])
except IndexError:
    cycle = None

if cycle is None:
    cycle = raw_input('Enter model initialization time (YYYYMMDDHH): ')

yyyy = int(cycle[0:4])
mm   = int(cycle[4:6])
dd   = int(cycle[6:8])
hh   = int(cycle[8:10])

date_str = datetime.datetime(yyyy,mm,dd,hh,0,0)


# Determine prod model
try:
    model_str = str(sys.argv[3])
except IndexError:
    model_str = None

if model_str is None:
    model_str = raw_input('Enter parallel model (HRRR or RAP): ')


# Get machine and head directory
machine, hostname = dawsonpy.get_machine()

if machine == 'WCOSS':
    DIR = '/gpfs/'+hostname[0]+'p2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'WCOSS_C':
    DIR = '/gpfs/hps3/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'WCOSS_DELL_P3':
    DIR = '/gpfs/dell2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'HERA':
    raise NameError, 'Need to modify htar_prod.py script to define data/script head directory'


# Make sure data directory doesn't fail
try:
    PROD_DIR = os.environ['PROD_DIR']
except KeyError:
    PROD_DIR = DIR+'/prod'
    print('PROD_DIR not defined in environment. Setting PROD_DIR to '+PROD_DIR)

if not os.path.exists(PROD_DIR):
    os.makedirs(PROD_DIR)

# Download data into temporary directory to ensure cycles are kept separate
DATA_DIR = os.path.join(PROD_DIR,model_str+'.'+cycle)
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)
os.chdir(DATA_DIR)


# Make sure script directory doesn't fail
try:
    SCRIPT_DIR = os.environ['SCRIPT_DIR']
except KeyError:
    SCRIPT_DIR = DIR+'/'+case+'/runscripts'
    print('SCRIPT_DIR not defined in environment. Setting SCRIPT_DIR to '+SCRIPT_DIR)



# Settings for RAPv4 and HRRRv3
if str.upper(model_str) == 'HRRR':
    RH_PREFIX = '/NCEPPROD/hpssprod/runhistory/rh'+cycle[0:4]+'/'+cycle[0:6]+'/'+cycle[0:8]
    FILE_PREFIX  = 'hrrr.t'+cycle[8:10]+'z.wrfprsf'
    FILE_PREFIX2 = 'hrrr.'+cycle[0:8]+'.t'+cycle[8:10]+'z.wrfprsf'
    TAR_BALL = 'wrf'
    prod_file_size = 300000000 
elif str.upper(model_str) == 'RAP':
    RH_PREFIX = '/NCEPPROD/hpssprod/runhistory/2year/rh'+cycle[0:4]+'/'+cycle[0:6]+'/'+cycle[0:8]
    FILE_PREFIX  = 'rap.t'+cycle[8:10]+'z.awp130pgrbf'
    FILE_PREFIX2 = 'rap.'+cycle[0:8]+'.t'+cycle[8:10]+'z.awp130pgrbf'
    TAR_BALL = 'awp130'
    prod_file_size = 10000000 

TAR_PREFIX = 'gpfs_hps_nco_ops_com_'+str.lower(model_str)+'_prod_'+str.lower(model_str)+'.'
FILE_SUFFIX  = '.grib2'



# Set up forecast hours to grab
hh = int(cycle[8:10])
if hh%6 == 0:
    runlength = 39
else:
    runlength = 21
if str.upper(model_str) == 'HRRR':
    runlength = runlength - 3

fhrs = np.arange(0,runlength+1,1)

date_list = [date_str + datetime.timedelta(hours=x) for x in fhrs]



# Loop through to create list of forecast files
req = {'htar_ball':'','htar_fname':[],'mv_from':[],'mv_to':[]}

print "Working on "+cycle+" "+str.upper(model_str)+" cycle"

for nf in range(len(date_list)):

    print "Adding "+str(fhrs[nf])+"-h forecast from "+cycle+" "+str.upper(model_str)+" cycle to extraction list"

    if int(cycle[8:]) <= 5:
        TAR_SUFFIX  = '00-05.'+TAR_BALL+'.tar'
    elif int(cycle[8:]) >= 6 and int(cycle[8:]) <= 11:
        TAR_SUFFIX  = '06-11.'+TAR_BALL+'.tar'
    elif int(cycle[8:]) >= 12 and int(cycle[8:]) <= 17:
        TAR_SUFFIX  = '12-17.'+TAR_BALL+'.tar'
    elif int(cycle[8:]) >= 18 and int(cycle[8:]) <= 23:
        TAR_SUFFIX  = '18-23.'+TAR_BALL+'.tar'

    if str.upper(model_str) == 'HRRR':
        TAR_SUFFIX = '_conus'+TAR_SUFFIX

    htar_ball = RH_PREFIX+'/'+TAR_PREFIX+cycle[0:8]+TAR_SUFFIX
    htar_fname = './'+FILE_PREFIX+str(fhrs[nf]).zfill(2)+FILE_SUFFIX
    mv_from = './'+FILE_PREFIX+str(fhrs[nf]).zfill(2)+FILE_SUFFIX
    mv_to = './'+FILE_PREFIX2+str(fhrs[nf]).zfill(2)+FILE_SUFFIX


    # Update dict
    req['htar_ball'] = htar_ball
    req['htar_fname'].append(htar_fname)
    req['mv_from'].append(mv_from)
    req['mv_to'].append(mv_to)


    # Delete prod file if path exists but file size is less than expected file size (trying to catch empty files)
    fname = FILE_PREFIX2+str(fhrs[nf]).zfill(2)+FILE_SUFFIX
    if os.path.exists(DATA_DIR+'/'+fname) and os.stat(DATA_DIR+'/'+fname).st_size < prod_file_size: 
        os.system('rm -f '+DATA_DIR+'/'+fname)
    elif os.path.exists(PROD_DIR+'/'+fname) and os.stat(PROD_DIR+'/'+fname).st_size < prod_file_size: 
        os.system('rm -f '+PROD_DIR+'/'+fname)

    # Make temp text file to use for HPSS request
    if not os.path.exists(DATA_DIR+'/'+fname) and not os.path.exists(PROD_DIR+'/'+fname):
        o = open("./"+str.lower(model_str)+"_"+cycle+"_temp.txt","w")
        o.write('\n'.join(req['htar_fname']))
        o.close()


# Submit HPSS request
os.system('htar -xvf '+req['htar_ball']+' -L '+str.lower(model_str)+'_'+cycle+'_temp.txt')

# Iterate through every item that was requested
for idx, (mv_from, mv_to) in enumerate(zip(req['mv_from'],req['mv_to'])):
    os.system("mv "+mv_from+" "+mv_to)

# Move prod files from temporary directory to PROD_DIR and delete temporary directory
os.chdir(PROD_DIR)
os.system('mv '+DATA_DIR+'/'+str.lower(model_str)+'*.grib2 .')
os.system('rm -fR '+DATA_DIR)

# Create done file that driver can find
os.system('touch '+SCRIPT_DIR+'/'+str.lower(model_str)+'_'+cycle+'_done')
print("Done")

