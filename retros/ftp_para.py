# Author: L Dawson
#
# Script to pull fcst files from HPSS, rename, and save in desired data directory

import numpy as np
import datetime, time, os, sys, subprocess
from ftplib import FTP
import dawsonpy


# Determine case name
try:
    case = str(sys.argv[1])
except IndexError:
    case = None

if case is None:
    case = raw_input('Enter a case name (do not included spaces): ')


# Determine cycle
try:
    cycle = str(sys.argv[2])
except IndexError:
    cycle = None

if cycle is None:
    cycle = raw_input('Enter model initialization time (YYYYMMDDHH): ')


# Determine para model
try:
    model_str = str(sys.argv[3])
except IndexError:
    model_str = None

if model_str is None:
    model_str = raw_input('Enter parallel model (HRRRX or RAPX): ')


# Get machine and head directory
machine, hostname = dawsonpy.get_machine()

if machine == 'WCOSS':
    DIR = '/gpfs/'+hostname[0]+'p2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'WCOSS_C':
    DIR = '/gpfs/hps3/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'WCOSS_DELL_P3':
    DIR = '/gpfs/dell2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
elif machine == 'HERA':
    raise NameError, 'Need to modify ftp_para.py script to define data/script head directory'


# Make sure data directory doesn't fail
try:
    PARA_DIR = os.environ['PARA_DIR']
except KeyError:
    PARA_DIR = DIR+'/para'
    print('PARA_DIR not defined in environment. Setting PARA_DIR to '+PARA_DIR)

if not os.path.exists(PARA_DIR):
    os.makedirs(PARA_DIR)

# Download data into temporary directory to ensure cycles are kept separate
DATA_DIR = os.path.join(PARA_DIR,model_str+'.'+cycle)
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)
os.chdir(DATA_DIR)


# Make sure script directory doesn't fail
try:
    SCRIPT_DIR = os.environ['SCRIPT_DIR']
except KeyError:
    SCRIPT_DIR = DIR+'/'+case+'/runscripts'
    print('SCRIPT_DIR not defined in environment. Setting SCRIPT_DIR to '+SCRIPT_DIR)


# Make sure retro period doesn't fail
try:
    RETRO_PERIOD = os.environ['RETRO_PERIOD']
except KeyError:
    if cycle < '20180901':
        if str.upper(model_str) == 'HRRRX':
            RETRO_PERIOD = 'julaug2018'    
        elif str.upper(model_str) == 'RAPX':
            RETRO_PERIOD = 'jul2018'    
    elif cycle > '20190101' and cycle < '20190401':
        RETRO_PERIOD = 'febmar2019'
    elif cycle > '20190401':
        if str.upper(model_str) == 'HRRRX':
            RETRO_PERIOD = 'aprmay2019'    
        elif str.upper(model_str) == 'RAPX':
            RETRO_PERIOD = 'may2019'    
    print('RETRO_PERIOD not defined in environment. Setting RETRO_PERIOD to '+RETRO_PERIOD)



print("Working on "+cycle+" "+str.upper(model_str)+" cycle")

ftp = FTP('gsdftp.fsl.noaa.gov')
ftp.login(user='anonymous',passwd=str.lower(os.environ['USER'])+'@noaa.gov')

# Change to correct directory on FTP server
if str.upper(model_str) == 'HRRRX':
   ftp_dir = 'retro/hrrr/HRRRv4_'+str.lower(RETRO_PERIOD)+'_retro1/'+cycle+'/wrftwo/hrconus'
   para_file_size = 100000000
elif str.upper(model_str) == 'RAPX':
   ftp_dir = 'retro/rap/RAPv5_'+str.lower(RETRO_PERIOD)+'_retro1/'+cycle+'/130/wrftwo'
   para_file_size = 12000000
ftp.cwd(ftp_dir)


# Get list of filenames on FTP server
listing = []
ftp.retrlines("LIST",listing.append)

# Loop through files, retrieve, and rename     
for line in listing:
    items = line.split()
    ftp_file = items[-1].lstrip()

    # Define local filename
    if str.upper(model_str) == 'HRRRX':
        local_fname = 'hrrr.'+cycle[0:8]+'.t'+cycle[8:10]+'z.wrfprsf'+ftp_file[-4:-2]+'.grib2'
    elif str.upper(model_str) == 'RAPX':
        local_fname = 'rap.'+cycle[0:8]+'.t'+cycle[8:10]+'z.awp130pgrbf'+ftp_file[-4:-2]+'.grib2'


    # Delete para file if path exists but file size is less than expected file size (trying to catch empty files)
    if os.path.exists(DATA_DIR+'/'+local_fname) and os.stat(DATA_DIR+'/'+local_fname).st_size < para_file_size:
        os.system('rm -f '+DATA_DIR+'/'+local_fname)
    elif os.path.exists(PARA_DIR+'/'+local_fname) and os.stat(PARA_DIR+'/'+local_fname).st_size < para_file_size:
        os.system('rm -f '+PARA_DIR+'/'+local_fname)


    # Download para file if it doesn't exist in DATA_DIR or PARA_DIR
    if not os.path.exists(DATA_DIR+'/'+local_fname) and not os.path.exists(PARA_DIR+'/'+local_fname):
        print("Retrieving f"+ftp_file[-4:-2]+" file")
        local_file = open(local_fname,'wb')

        try:
            ftp.retrbinary('RETR '+ftp_file,local_file.write)
            local_file.close()
        except:
            print("Likely timed out. Going to try logging in and trying again.")
            ftp = FTP('gsdftp.fsl.noaa.gov')
            ftp.login(user='anonymous',passwd=str.lower(os.environ['USER'])+'@noaa.gov')
            ftp.cwd(ftp_dir)
            ftp.retrbinary('RETR '+ftp_file,local_file.write)
            local_file.close()


ftp.quit()

# Move para files from temporary directory to PARA_DIR and delete temporary directory
os.chdir(PARA_DIR)
os.system('mv '+DATA_DIR+'/*.grib2 .')
os.system('rm -fR '+DATA_DIR)

# Create done file that driver can find
os.system('touch '+SCRIPT_DIR+'/'+str.lower(model_str)+'_'+cycle+'_done')
print("Done")

