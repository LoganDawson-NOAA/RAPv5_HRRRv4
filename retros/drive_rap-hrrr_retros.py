import sys
import datetime
import shutil
import os, re
import subprocess
import numpy as np
import time
import dawsonpy


# Function to write HPSS job scripts
def write_download_job(data):

    if data[0:4] == 'hrrr':
        data = data+'_'+cycle
    elif data == 'rapv4' or data == 'rapx':
        data = data+'_'+cycle

    with open(SCRIPT_DIR+'/download_'+data+'.sh','w') as f:
        f.write("#!/bin/ksh --login\n")
        f.write("#\n")
        f.write("#BSUB -J download_"+data+"\n")
        f.write("#BSUB -o "+SCRIPT_DIR+"/"+data+".out\n")
        f.write("#BSUB -e "+SCRIPT_DIR+"/"+data+".err\n")
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -W 02:00\n")
        f.write("#BSUB -P HRRR-T2O\n")

        if machine == 'WCOSS':
            f.write("#BSUB -q transfer \n")
        elif machine == 'WCOSS_DELL_P3':
            f.write("#BSUB -q dev_transfer \n")

        f.write("#BSUB -R \"rusage[mem=1000]\"\n")
        f.write("#BSUB -R \"affinity[core]\"\n")
        f.write("#\n")
        f.write("set +x\n")
        f.write("source ~/.bashrc\n")
        f.write("set -x\n")
        f.write("\n")
        if str.lower(data[0:5]) == 'hrrrx':
            f.write("python "+DRIVE_DIR+"/ftp_para.py "+case+" "+cycle+" "+str.upper(data[0:5])+"\n")
        elif str.lower(data[0:4]) == 'rapx':
            f.write("python "+DRIVE_DIR+"/ftp_para.py "+case+" "+cycle+" "+str.upper(data[0:4])+"\n")
        elif str.lower(data[0:4]) == 'hrrr':
            f.write("python "+DRIVE_DIR+"/htar_prod.py "+case+" "+cycle+" "+str.upper(data[0:4])+"\n")
        elif str.lower(data[0:5]) == 'rapv4':
            f.write("python "+DRIVE_DIR+"/htar_prod.py "+case+" "+cycle+" "+str.upper(data[0:3])+"\n")
        elif str.lower(data) == 'nohrsc':
            f.write("python "+DRIVE_DIR+"/wget_nohrsc.py "+case+" "+day1_00z_cycle+" "+str.upper(data)+"\n")
        else:
            f.write("python "+DRIVE_DIR+"/htar_"+data+".py "+case+" "+day1_00z_cycle+"\n")
        f.write("\n")
        f.write("exit\n")
        f.write("\n")

        f.close()



# Function to write plotting scripts
def write_plotting_job(jobname,model,cycle,fhr):

    # Set up graphx directory for each cycle
    GRAPHX_DIR = os.path.join(GRAPHX_HEAD,model,cycle)

    with open(SCRIPT_DIR+'/'+jobname+'.sh','w') as f:
        f.write("#!/bin/ksh --login\n")
        f.write("#\n")
        f.write("#BSUB -J "+jobname+"\n")
        f.write("#BSUB -o "+SCRIPT_DIR+"/"+jobname+".out\n")
        f.write("#BSUB -e "+SCRIPT_DIR+"/"+jobname+".out\n")

        if jobname[0:5] == 'plot0' or jobname[0:5] == 'plot1' or jobname[0:5] == 'plot2':
            f.write("#BSUB -n 2\n")
            f.write("#BSUB -R span[ptile=1]\n")
        else:
            f.write("#BSUB -n 16\n")
            f.write("#BSUB -R span[ptile=4]\n")

        f.write("#BSUB -P "+str.upper(model)+"-T2O\n")

        if str.upper(machine) == 'WCOSS':
            f.write("#BSUB -W 02:00\n")
            f.write("#BSUB -q \"dev2\" \n")
            f.write("#BSUB -R \"affinity[core]\"\n")
            f.write("#BSUB -x\n")
            f.write("#BSUB -a poe\n")

        elif str.upper(machine) == 'WCOSS_DELL_P3':
            # wallclock for QPF and UH plotting
            if jobname[0:5] == 'plot0' or jobname[0:5] == 'plot2':
                f.write("#BSUB -W 02:00\n")

            # wallclock for snow plotting
            elif jobname[0:5] == 'plot1':
                f.write("#BSUB -W 01:10\n")

            # wallclock for plotting by forecast hour
            else:
                f.write("#BSUB -W 02:30\n")

            f.write("#BSUB -q dev_shared \n")
            f.write("#BSUB -R \"affinity[core]\"\n")
            f.write("#BSUB -R \"rusage[mem=6000]\"\n")

        f.write("#\n")
        f.write("set +x\n")
        f.write("\n")
        f.write("module use -a /u/Benjamin.Blake/modulefiles\n")
        f.write("module load anaconda2/latest\n")
        f.write("export GRIB_DEFINITION_PATH="+GRIB_DEFINITION_PATH+"\n")
        f.write("export PYTHONPATH=${PYTHONPATH}:"+NCEPY_PATH+":"+DAWSONPY_PATH+"\n")
        f.write("\n")
        f.write("set -x\n")
        f.write("\n")
        f.write("mkdir -p "+GRAPHX_DIR+"\n")
        f.write("cd "+GRAPHX_DIR+"\n")
        f.write("\n")

        # Do run-total QPF plots
        if jobname[0:5] == 'plot0':
            transfer_script = "transfer0"
            os.environ['RUN_LENGTH'] = str(runlength)
            f.write("python "+DRIVE_DIR+"/plot_retro_alldomains2_qpf.py "+case+" "+cycle+" "+model+"\n")

        # Do run-total snow plots
        elif jobname[0:5] == 'plot1':
            transfer_script = "transfer1"
            os.environ['RUN_LENGTH'] = str(runlength)
            f.write("python "+DRIVE_DIR+"/plot_retro_alldomains2_snow.py "+case+" "+cycle+" "+model+"\n")

        # Do run-total UH plots
        elif jobname[0:5] == 'plot2':
            transfer_script = "transfer2"
            os.environ['RUN_LENGTH'] = str(runlength)
            f.write("python "+DRIVE_DIR+"/plot_retro_alldomains2_uh.py "+case+" "+cycle+" "+model+"\n")

        # Plot F00 thorugh F09
        elif jobname[0:5] == 'plot3':
            transfer_script = "transfer3"
            for fhr in np.arange(0,10,1):
                f.write("python "+DRIVE_DIR+"/plot_retro_alldomains1.py "+case+" "+cycle+" "+str(fhr).zfill(2)+" "+model+"\n")

        # Plot F10 thorugh F19
        elif jobname[0:5] == 'plot4':
            transfer_script = "transfer4"
            for fhr in np.arange(10,20,1):
                f.write("python "+DRIVE_DIR+"/plot_retro_alldomains1.py "+case+" "+cycle+" "+str(fhr).zfill(2)+" "+model+"\n")

        # Plot F20 thorugh F29
        elif jobname[0:5] == 'plot5':
            transfer_script = "transfer5"
            for fhr in np.arange(20,30,1):
                f.write("python "+DRIVE_DIR+"/plot_retro_alldomains1.py "+case+" "+cycle+" "+str(fhr).zfill(2)+" "+model+"\n")

        # Plot F30 thorugh F39
        elif jobname[0:5] == 'plot6':
            transfer_script = "transfer6"
            for fhr in np.arange(30,40,1):
                f.write("python "+DRIVE_DIR+"/plot_retro_alldomains1.py "+case+" "+cycle+" "+str(fhr).zfill(2)+" "+model+"\n")

            

        f.write("\n")
        f.write("bsub < "+SCRIPT_DIR+"/"+transfer_script+"_"+str.lower(model)+"_"+cycle+".sh\n")
        f.write("\n")
        f.write("exit\n")
        f.write("\n")

        f.close()



# Function to write rzdm tranfer job scripts
def write_transfer_job(jobname,model,cycle):

    if date_str.strftime('%Y%m%d') < '20190401':
        retro_dir = 'feb2019'
    elif date_str.strftime('%Y%m%d') > '20190401' and date_str.strftime('%Y%m%d') < '20190615':
        retro_dir = 'may2019'
    elif date_str.strftime('%Y%m%d') > '20190615':
        retro_dir = 'jul2019'

    with open(SCRIPT_DIR+'/transfer'+jobname[4]+'_'+str.lower(model)+'_'+cycle+'.sh','w') as f:
        f.write("#!/bin/ksh --login\n")
        f.write("#\n")
        f.write("#BSUB -J scp"+jobname[4]+"_"+str.lower(model)+"_"+cycle+"\n")
        f.write("#BSUB -o "+SCRIPT_DIR+"/scp"+jobname[4]+"_"+str.lower(model)+"_"+cycle+".out\n")
        f.write("#BSUB -e "+SCRIPT_DIR+"/scp"+jobname[4]+"_"+str.lower(model)+"_"+cycle+".out\n")
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -W 00:10\n")
        f.write("#BSUB -P "+str.upper(model)+"-T2O\n")

        if str.upper(machine) == 'WCOSS':
            f.write("#BSUB -q transfer\n")
        elif str.upper(machine) == 'WCOSS_DELL_P3':
            f.write("#BSUB -q dev_transfer\n")

        f.write("#BSUB -R \"rusage[mem=300]\"\n")
        f.write("#BSUB -R \"affinity[core]\"\n")
        f.write("#\n")
        f.write("set -x\n")
        f.write("\n")
        f.write("cd "+GRAPHX_HEAD+"/"+model+"/"+cycle+"\n")
        f.write("\n")
        f.write("ssh "+RZDM_USER+"@emcrzdm.ncep.noaa.gov \"mkdir -p "+RZDM_DIR+"/"+retro_dir+"/"+case+"/"+str.lower(model)+"/"+cycle+"\"\n")
        
        if jobname[0:5] == 'plot0':
            files = 'compareqpf*.png'
        elif jobname[0:5] == 'plot1':
            files = 'comparesnow*.png'
        elif jobname[0:5] == 'plot2':
            files = 'compareuh*.png'
        elif jobname[0:5] == 'plot3':
            files = '*f0*.png'
        elif jobname[0:5] == 'plot4':
            files = '*f1*.png'
        elif jobname[0:5] == 'plot5':
            files = '*f2*.png'
        elif jobname[0:5] == 'plot6':
            files = '*f3*.png'

        f.write("scp "+files+" "+RZDM_USER+"@emcrzdm.ncep.noaa.gov:"+RZDM_DIR+"/"+retro_dir+"/"+case+"/"+str.lower(model)+"/"+cycle+"\n")

        f.write("\n")
        f.write("exit\n")
        f.write("\n")

        f.close()




# Update template php for rzdm
def create_case_php(template_php, case_php):

    replacements = {
        "%RETRO_PERIOD%"   : str.lower(retro_dir),
        "%CASE%"           : str.lower(case),
        "%DAY1_00Z_CYCLE%" : day1_00z_cycle,
	"%DATE_M1%"        : cycle_list[0].strftime('%m/%d/%y'),
	"%ymd_m1%"         : cycle_list[0].strftime('%Y%m%d'),
        "%DATE%"           : cycle_list[2].strftime('%m/%d/%y'),
        "%ymd%"            : cycle_list[2].strftime('%Y%m%d'),
        "%DATE_P1%"        : cycle_list[-1].strftime('%m/%d/%y'),
        "%ymd_p1%"         : cycle_list[-1].strftime('%Y%m%d'),
    }

    with open(template_php) as f:
        new_text = dawsonpy.multiple_replace(replacements, f.read())

    with open(case_php, "w") as result:
        result.write(new_text)



#########################################################################################
########                        BEGIN WORKING                                   #########
#########################################################################################


# Determine case name
try:
    case = str(sys.argv[1])
except IndexError:
    case = None

if case is None:
    case = raw_input('Enter a case name (do not included spaces): ')


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

nhrs = [-12,-6]
day1_cycles = np.arange(0,25,3)
nhrs.extend(day1_cycles)
cycle_list = [date_str + datetime.timedelta(hours=x) for x in nhrs]


# Create list of valid times
final_vtime = cycle_list[-1] + datetime.timedelta(hours=36)

vtime = cycle_list[0]
x=0
vtime_list=[]
vtime_list2=[]
while vtime < final_vtime:
    vtime = cycle_list[0] + datetime.timedelta(hours=x)
    vtime_str = vtime.strftime('%Y%m%d%H')
    vtime_list.append(vtime_str)
    if int(vtime_str[8:10])%6 == 0:
        vtime_list2.append(vtime_str)
    x += 1


# Define retro period for emcrzdm use
if date_str.strftime('%Y%m%d') < '20180901':
    retro_dir = 'julaug2018'
elif date_str.strftime('%Y%m%d') > '20190101' and date_str.strftime('%Y%m%d') < '20190401':
    retro_dir = 'febmar2019'
elif date_str.strftime('%Y%m%d') > '20190401':
    retro_dir = 'may2019'


# Get machine
machine, hostname = dawsonpy.get_machine()


# Set up working directory
if machine == 'WCOSS':
    DIR = '/gpfs/'+hostname[0]+'p2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
    GRIB_DEFINITION_PATH = '/gpfs/'+hostname[0]+'d1/emc/meso/save/Benjamin.Blake/EXT/grib_api.1.14.4/share/grib_api/definitions'
    NCEPY_PATH = '/gpfs/'+hostname[0]+'d1/emc/meso/save/Jacob.Carley/python/lib'
elif machine == 'WCOSS_C':
    DIR = '/gpfs/hps3/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
    GRIB_DEFINITION_PATH = '/gpfs/hps3/emc/meso/save/Benjamin.Blake/EXT/grib_api-1.15.0/share/grib_api/definitions'
    NCEPY_PATH = '/gpfs/'+hostname[0]+'d1/emc/meso/save/Jacob.Carley/python/lib'
elif machine == 'WCOSS_DELL_P3':
    DIR = '/gpfs/dell2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
    GRIB_DEFINITION_PATH = '/gpfs/dell2/emc/modeling/noscrub/Benjamin.Blake/EXT/grib_api.1.14.4/share/grib_api/definitions'
    NCEPY_PATH = '/gpfs/dell2/emc/modeling/noscrub/Jacob.Carley/python/lib'
    DAWSONPY_PATH = '/gpfs/dell2/emc/verification/noscrub/Logan.Dawson/python'
elif machine == 'HERA':
    DIR = raw_input('Enter scratch directory where you want to save data and plots (do not include spaces): ')
    GRIB_DEFINITION_PATH = '/need/to/define'
    NCEPY_PATH = '/scratch2/NCEPDEV/fv3-cam/Jacob.Carley/python/lib64/python'

if not os.path.exists(DIR):
    os.makedirs(DIR)

# Set up directory for production data
PROD_DIR = os.path.join(DIR,'prod')
if not os.path.exists(PROD_DIR):
    os.makedirs(PROD_DIR)
print(" ")
print("RAPv4 and HRRRv3 data will be downloaded into "+PROD_DIR)

# Set up directory for production data
PARA_DIR = os.path.join(DIR,'para')
if not os.path.exists(PARA_DIR):
    os.makedirs(PARA_DIR)
print("RAPv5 and HRRRv4 data will be downloaded into "+PARA_DIR)

# Set up directory for obs data
OBS_DIR = os.path.join(DIR,'obs')
if not os.path.exists(OBS_DIR):
    os.makedirs(OBS_DIR)
print("Observation data will be downloaded into "+OBS_DIR)

# Set up script directory
SCRIPT_DIR = os.path.join(DIR,case,'runscripts')
if os.path.exists(SCRIPT_DIR):
    shutil.rmtree(SCRIPT_DIR)
os.makedirs(SCRIPT_DIR)
print("Batch job scripts and output files will be written to "+SCRIPT_DIR)

# Set up graphx directory
GRAPHX_HEAD = os.path.join(DIR,case,'graphx')
if not os.path.exists(GRAPHX_HEAD):
    os.makedirs(GRAPHX_HEAD)
print("Plots will be saved in "+GRAPHX_HEAD)
print(" ")
print(" ")

# Get path to directory where driver scripts are 
DRIVE_DIR = os.getcwd()

# Set environment variables for downstream scripts to use
os.environ['DRIVE_DIR'] = DRIVE_DIR
os.environ['PROD_DIR'] = PROD_DIR
os.environ['PARA_DIR'] = PARA_DIR
os.environ['OBS_DIR'] = OBS_DIR
os.environ['SCRIPT_DIR'] = SCRIPT_DIR
os.environ['GRAPHX_HEAD'] = GRAPHX_HEAD


# Settings for emcrzdm
RZDM_DIR = '/home/people/emc/www/htdocs/users/meg/rapv5_hrrrv4/retros'
if os.environ['USER'] == 'Logan.Dawson':
     RZDM_USER = 'ldawson'
elif os.environ['USER'] == 'Geoffrey.Manikin':
     RZDM_USER = 'wd20mg'
elif os.environ['USER'] == 'Alicia.Bentley':
     RZDM_USER = 'abentley'
elif os.environ['USER'] == 'Shannon.Shields':
     RZDM_USER = 'shannon.shields'

try:
    RZDM_USER
except NameError:
    RZDM_USER = None

if RZDM_USER is None:
    RZDM_USER = raw_input('Enter emcrzdm user name for transfer jobs: ')


#####################################################################

# Determine what analysis data need to be downloaded (if any)
anls = ['rap','urma','refl','st4','nohrsc']
products = anls

# Check to see if all RAP analysis files exist
rap_exists    = [True if os.path.exists(OBS_DIR+'/rap/rap.'+vtime[0:8]+'.t'+vtime[8:10]+'z.awp130pgrbf00.grib2') and \
                 os.stat(OBS_DIR+'/rap/rap.'+vtime[0:8]+'.t'+vtime[8:10]+'z.awp130pgrbf00.grib2').st_size > 10000000 else False for vtime in vtime_list]

# Check to see if all URMA analysis files exist
urma_exists   = [True if os.path.exists(OBS_DIR+'/urma/urma2p5.'+vtime[0:8]+'.t'+vtime[8:10]+'z.2dvaranl_ndfd.grb2') and \
                 os.stat(OBS_DIR+'/urma/urma2p5.'+vtime[0:8]+'.t'+vtime[8:10]+'z.2dvaranl_ndfd.grb2').st_size > 80000000 else False for vtime in vtime_list]

# Check to see if all reflecitvity obs files exist
refl_exists   = [True if os.path.exists(OBS_DIR+'/refl/refd3d.'+vtime[0:8]+'.t'+vtime[8:10]+'z.grb2f00') and \
                 os.stat(OBS_DIR+'/refl/refd3d.'+vtime[0:8]+'.t'+vtime[8:10]+'z.grb2f00').st_size > 2500000 else False for vtime in vtime_list]

# Check to see if all Stage IV analysis files exist
st4_exists    = [True if os.path.exists(OBS_DIR+'/st4/ST4.'+vtime+'.01h') and \
                 os.stat(OBS_DIR+'/st4/ST4.'+vtime+'.01h').st_size > 900000 else False for vtime in vtime_list]

# Check to see if all URMA analysis files exist
nohrsc_exists = [True if os.path.exists(OBS_DIR+'/nohrsc/sfav2_CONUS_6h_'+vtime+'_grid184.grb2') and \
                 os.stat(OBS_DIR+'/nohrsc/sfav2_CONUS_6h_'+vtime+'_grid184.grb2').st_size > 300000 else False for vtime in vtime_list2]

download_rap  = True
download_urma = True
download_refl = True
download_st4  = True
download_nohrsc  = True

rap_done = os.path.join(SCRIPT_DIR,'rap_done')
urma_done = os.path.join(SCRIPT_DIR,'urma_done')
refl_done = os.path.join(SCRIPT_DIR,'refl_done')
st4_done = os.path.join(SCRIPT_DIR,'st4_done')
nohrsc_done = os.path.join(SCRIPT_DIR,'nohrsc_done')

if all(rap_exists):
    download_rap = False
    os.system('touch '+rap_done)
    print("All RAP analysis data previously downloaded.")
if all(urma_exists):
    download_urma = False
    os.system('touch '+urma_done)
    print("All URMA analysis data previously downloaded.")
if all(refl_exists):
    download_refl = False
    os.system('touch '+refl_done)
    print("All reflectivity data previously downloaded.")
if all(st4_exists):
    download_st4 = False
    os.system('touch '+st4_done)
    print("All Stage IV analysis data previously downloaded.")
if all(nohrsc_exists):
    download_nohrsc = False
    os.system('touch '+nohrsc_done)
    print("All NOHRSC analysis data previously downloaded.")


# Downlod analysis data from HPSS if needed
download = [download_rap, download_urma, download_refl, download_st4, download_nohrsc]
if any(download) is True:
    download_done = False
    j = 0
    for data in anls:
        job_script = SCRIPT_DIR+'/download_'+data+'.sh'

        if download[j] is True:
            print(" ")
            print("Downloading "+str.upper(data)+".")
            write_download_job(data)
            dawsonpy.submit_job(job_script)
            time.sleep(5)  # wait for 5 seconds between each job submission
        j += 1
else:
    download_done = True
    print("All analysis data previously downloaded.")


# Determine what forecast data need to be downloaded (if any)
download_hrrr = []
download_hrrrx = []
download_rapv4 = []
download_rapx = []
for this_cycle in cycle_list:

    cycle = this_cycle.strftime('%Y%m%d%H')

    # Determine forecast hours for HRRR/X cycle
    hh = int(cycle[8:10])
    if hh%12 == 0:
        runlength = 36 
    elif hh%12 != 0:
        runlength = 18
    fhrs = np.arange(0,runlength+1,1)

    # Check to see if all HRRR/X forecast files exist for current cycle
    hrrr_cycle_exists  = [True if os.path.exists(PROD_DIR+'/hrrr.'+cycle[0:8]+'.t'+cycle[8:10]+'z.wrfprsf'+str(fhr).zfill(2)+'.grib2') and \
                         os.stat(PROD_DIR+'/hrrr.'+cycle[0:8]+'.t'+cycle[8:10]+'z.wrfprsf'+str(fhr).zfill(2)+'.grib2').st_size > 300000000 else False for fhr in fhrs]
    hrrrx_cycle_exists = [True if os.path.exists(PARA_DIR+'/hrrr.'+cycle[0:8]+'.t'+cycle[8:10]+'z.wrfprsf'+str(fhr).zfill(2)+'.grib2') and \
                         os.stat(PARA_DIR+'/hrrr.'+cycle[0:8]+'.t'+cycle[8:10]+'z.wrfprsf'+str(fhr).zfill(2)+'.grib2').st_size > 100000000 else False for fhr in fhrs]

    # Determine forecast hours for RAP/X cycle
    if hh%12 == 9:
        runlength = 39
    elif hh%12 != 9:
        runlength = 21 
    fhrs = np.arange(0,runlength+1,1)

    # Check to see if all RAP/X forecast files exist for current cycle
    rapv4_cycle_exists = [True if os.path.exists(PROD_DIR+'/rap.'+cycle[0:8]+'.t'+cycle[8:10]+'z.awp130pgrbf'+str(fhr).zfill(2)+'.grib2') and \
                         os.stat(PROD_DIR+'/rap.'+cycle[0:8]+'.t'+cycle[8:10]+'z.awp130pgrbf'+str(fhr).zfill(2)+'.grib2').st_size > 10000000 else False for fhr in fhrs]
    rapx_cycle_exists  = [True if os.path.exists(PARA_DIR+'/rap.'+cycle[0:8]+'.t'+cycle[8:10]+'z.awp130pgrbf'+str(fhr).zfill(2)+'.grib2') and \
                         os.stat(PARA_DIR+'/rap.'+cycle[0:8]+'.t'+cycle[8:10]+'z.awp130pgrbf'+str(fhr).zfill(2)+'.grib2').st_size > 10000000 else False for fhr in fhrs]

    download_hrrr_cycle  = True
    download_hrrrx_cycle = True
    download_rapv4_cycle  = True
    download_rapx_cycle = True

    # Define paths for writing 'done' files
    hrrr_cycle_done = os.path.join(SCRIPT_DIR,'hrrr_'+cycle+'_done')
    hrrrx_cycle_done = os.path.join(SCRIPT_DIR,'hrrrx_'+cycle+'_done')
    rap_cycle_done = os.path.join(SCRIPT_DIR,'rapv4_'+cycle+'_done')
    rapx_cycle_done = os.path.join(SCRIPT_DIR,'rapx_'+cycle+'_done')

    # Create 'done' file for HRRR/X and RAP/X cycles if forecast files for each cycle exist
    if all(hrrr_cycle_exists):
        download_hrrr_cycle = False
        os.system('touch '+hrrr_cycle_done)
    if all(hrrrx_cycle_exists):
        download_hrrrx_cycle = False
        os.system('touch '+hrrrx_cycle_done)
    if all(rapv4_cycle_exists):
        download_rapv4_cycle = False
        os.system('touch '+rap_cycle_done)
    if all(rapx_cycle_exists):
        download_rapx_cycle = False
        os.system('touch '+rapx_cycle_done)
    
    # Create list (with size = number of cycles) definining which cycles need to be downloaded
    download_hrrr.append(download_hrrr_cycle)
    download_hrrrx.append(download_hrrrx_cycle)
    download_rapv4.append(download_rapv4_cycle)
    download_rapx.append(download_rapx_cycle)


# Download prod and para forecast data from HPSS/FTP if needed

# These two lines indicate we will only download HRRR and HRRRX data.
download = [download_hrrr, download_hrrrx]
fcsts = ['hrrr','hrrrx']

# Uncomment these next two lines if you want to download RAP and RAPX data in addition to HRRR and HRRRX data.
#download = [download_hrrr, download_hrrrx, download_rapv4, download_rapx]
#fcsts = ['hrrr','hrrrx','rapv4','rapx']

products.extend(fcsts)
hrrr_done = os.path.join(SCRIPT_DIR,'hrrr_done')
hrrrx_done = os.path.join(SCRIPT_DIR,'hrrrx_done')
rapv4_done = os.path.join(SCRIPT_DIR,'rapv4_done')
rapx_done = os.path.join(SCRIPT_DIR,'rapx_done')
i = 0
for fcst in fcsts:
    if any(download[i]) is True:
        download_done = False
    
        j = 0
        for this_cycle in cycle_list:
            cycle = this_cycle.strftime('%Y%m%d%H')
            job_script = SCRIPT_DIR+'/download_'+fcst+'_'+cycle+'.sh'

            if download[i][j] is True:
                print(" ")
                print("Downloading "+cycle+" "+str.upper(fcst)+" cycle.")
                write_download_job(fcst)
                dawsonpy.submit_job(job_script)
                time.sleep(5)  # wait for 5 seconds between each job submission

            j += 1
    else:
        download_done = False
        print("All "+str.upper(fcst)+" data previously downloaded.")
        if fcst == 'hrrr':
            os.system('touch '+rapx_done)
        elif fcst == 'hrrrx':
            os.system('touch '+hrrrx_done)
        elif fcst == 'rapv4':
            os.system('touch '+rapv4_done)
        elif fcst == 'rapx':
            os.system('touch '+rapx_done)

    i += 1


# Wait until all data are downloaded to continue
while download_done is False:

    hrrr_done_files_exist = [True if os.path.exists(SCRIPT_DIR+'/hrrr_'+cycle.strftime('%Y%m%d%H')+'_done') else False for cycle in cycle_list]    
    if all(hrrr_done_files_exist):
        os.system('touch '+hrrr_done)

    hrrrx_done_files_exist = [True if os.path.exists(SCRIPT_DIR+'/hrrrx_'+cycle.strftime('%Y%m%d%H')+'_done') else False for cycle in cycle_list]    
    if all(hrrrx_done_files_exist):
        os.system('touch '+hrrrx_done)

    rap_done_files_exist = [True if os.path.exists(SCRIPT_DIR+'/rapv4_'+cycle.strftime('%Y%m%d%H')+'_done') else False for cycle in cycle_list]    
    if all(rap_done_files_exist):
        os.system('touch '+rapv4_done)

    rapx_done_files_exist = [True if os.path.exists(SCRIPT_DIR+'/rapx_'+cycle.strftime('%Y%m%d%H')+'_done') else False for cycle in cycle_list]    
    if all(rapx_done_files_exist):
        os.system('touch '+rapx_done)

    done_files_exist = [True if os.path.exists(SCRIPT_DIR+'/'+product+'_done') else False for product in products] 
    if all(done_files_exist) is True:
        print(" ")
        print(" ")
        print("All download jobs finished. Prepare to plot.")
        print(" ")
        download_done = True
    else:
        print(" ")
        print(" ")
        print("Waiting for download jobs to finish. You may need to re-run driver to download all data from HPSS/FTP server.")
        print("Once all data exists in prod/para/obs directories, re-run driver script to begin plotting.")
        print(" ")
        download_done = False
        exit()



# Plot retro graphics 
if 'rapv4' in fcsts and 'rapx' in fcsts:
    models = ['HRRR','RAP']
else:
    models = ['HRRR']

# If you need to plot a fewer number of cycles as a time, you can do it using somethng like the line below
#cycle_list = ['2019042912','2019042918','2019043000','2019043003']

# Loop through models and cycles to create/submit plotting jobs
for model in models:
    for this_cycle in cycle_list:

        cycle = this_cycle.strftime('%Y%m%d%H')
        print(" ")
        print "Working on "+cycle+" "+str.upper(model)+" cycle"


        # Set up forecast hours to plot
        hh = int(cycle[8:10])
        if str.upper(model) == 'HRRR' and hh%12 == 0:
            runlength = 36
            njobs = 4 
        elif str.upper(model) == 'HRRR' and hh%12 != 0:
            runlength = 18
            njobs = 2 
        elif str.upper(model) == 'RAP' and hh%12 == 9:
            runlength = 39
            njobs = 4 
        elif str.upper(model) == 'RAP' and hh%12 != 9:
            runlength = 21
            njobs = 3 
        fhrs = np.arange(0,runlength+1,1)

        # Add three more jobs for QPF, snow, and UH accum plotting scripts
        njobs = njobs + 3

        # Write batch jobs that plot graphics and transfer to rzdm
        # Submit plotting jobs. Each plotting job submits its transfer job once done.
        j = 0
        while j < njobs:

            print(" ")
            if j == 0:
                print("Submitting batch job to generate run-total QPF plots for "+cycle+" "+str.upper(model)+" cycle.")
            elif j == 1:
                print("Submitting batch job to generate run-total snow plots for "+cycle+" "+str.upper(model)+" cycle.")
            elif j == 2:
                print("Submitting batch job to generate run-total UH plots for "+cycle+" "+str.upper(model)+" cycle.")
            elif j == 3:
                print("Submitting batch job to generate F00 to F09 plots for "+cycle+" "+str.upper(model)+" cycle.")
            elif j == 4:
                print("Submitting batch job to generate F10 to F19 plots for "+cycle+" "+str.upper(model)+" cycle.")
            elif j == 5:
                print("Submitting batch job to generate F20 to F29 plots for "+cycle+" "+str.upper(model)+" cycle.")
            elif j == 6:
                print("Submitting batch job to generate F30 to F39 plots for "+cycle+" "+str.upper(model)+" cycle.")

            jobname = 'plot'+str(j)+'_'+str.lower(model)+'_'+cycle
            write_plotting_job(jobname,model,cycle,fhrs)
            write_transfer_job(jobname,model,cycle)
            job_script = SCRIPT_DIR+'/plot'+str(j)+'_'+str.lower(model)+'_'+cycle+'.sh'
            dawsonpy.submit_job(job_script)
            time.sleep(5)  # wait for 5 seconds between each job submission

            j += 1

        time.sleep(30)  # wait an extra 30 seconds between each cycle




# Update template php and copy to emcrzdm 
print(" ")
print(" ")
print("Generating new php looper and copying to emcrzdm. Look in "+RZDM_DIR+"/"+str.lower(retro_dir)+"/"+str.lower(case)+" for looper and images.")
print("Domains for "+case+" case are:")
print(dawsonpy.get_domains(case,model))
print(" ")
create_case_php('template.php',str.lower(case)+'.php')
os.system("ssh "+RZDM_USER+"@emcrzdm.ncep.noaa.gov \"mkdir -p "+RZDM_DIR+"/"+str.lower(retro_dir)+"/"+str.lower(case)+"\"")
os.system("scp "+str.lower(case)+".php "+RZDM_USER+"@emcrzdm.ncep.noaa.gov:"+RZDM_DIR+"/"+str.lower(retro_dir)+"/"+str.lower(case))
os.system("ssh "+RZDM_USER+"@emcrzdm.ncep.noaa.gov \"chmod -R 775 "+RZDM_DIR+"/"+str.lower(retro_dir)+"/"+str.lower(case)+"\"")


exit()

