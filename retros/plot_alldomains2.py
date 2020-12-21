import pygrib
import matplotlib
matplotlib.use('Agg')
import cStringIO
import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.image as image
from matplotlib.gridspec import GridSpec
from mpl_toolkits.basemap import Basemap
import numpy as np
import time,os,sys,multiprocessing,itertools
import ncepy
import csv
from scipy import ndimage
#from netCDF4 import Dataset

#--------------Define some functions ------------------#

def clear_plotables(ax,keep_ax_lst,fig):
  #### - step to clear off old plottables but leave the map info - ####
  if len(keep_ax_lst) == 0 :
    print "clear_plotables WARNING keep_ax_lst has length 0. Clearing ALL plottables including map info!"
  cur_ax_children = ax.get_children()[:]
  if len(cur_ax_children) > 0:
    for a in cur_ax_children:
      if a not in keep_ax_lst:
       # if the artist isn't part of the initial set up, remove it
        a.remove()

def compress_and_save(filename):
  #### - compress and save the image - ####
# plt.savefig(filename, format='png', bbox_inches='tight', dpi=150)
  ram = cStringIO.StringIO()
  plt.savefig(ram, format='png', bbox_inches='tight', dpi=150)
  ram.seek(0)
  im = Image.open(ram)
  im2 = im.convert('RGB').convert('P', palette=Image.ADAPTIVE)
  im2.save(filename, format='PNG')




#-------------------------------------------------------#

# Necessary to generate figs when not running an Xserver (e.g. via PBS)
# plt.switch_backend('agg')

# Read date/time and forecast hour from command line
case = str(sys.argv[1])
cycle = str(sys.argv[2])
ymd = cycle[0:8]
year = int(cycle[0:4])
month = int(cycle[4:6])
day = int(cycle[6:8])
hour = int(cycle[8:10])
cyc = str(hour).zfill(2)
print year, month, day, hour


# Get machine and head directory
machine, hostname = dawsonpy.get_machine()

if machine == 'WCOSS':
    DIR = '/gpfs/'+hostname[0]+'p2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
    REP_DIR = '/gpfs/'+hostname[0]+'d1/emc/meso/noscrub/Logan.Dawson/RAPv5_HRRRv4_retros/reports'
elif machine == 'WCOSS_C':
    DIR = '/gpfs/hps3/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
elif machine == 'WCOSS_DELL_P3':
    DIR = '/gpfs/dell2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
    REP_DIR = '/gpfs/dell2/emc/verification/noscrub/Logan.Dawson/RAPv5_HRRRv4_retros/reports'
elif machine == 'HERA':
    raise NameError, 'Need to modify plot_alldomains2.py script to define data/graphx head directory'


# Make sure prod directory doesn't fail
try:
    PROD_DIR = os.environ['PROD_DIR']
except KeyError:
    print('Prod directory not defined in environment')
    PROD_DIR = DIR+'/'+case+'/prod'

# Make sure para directory doesn't fail
try:
    PARA_DIR = os.environ['PARA_DIR']
except KeyError:
    print('Para directory not defined in environment')
    PARA_DIR = DIR+'/'+case+'/para'

# Make sure obs directory doesn't fail
try:
    OBS_DIR = os.environ['OBS_DIR']
except KeyError:
    print('Obs directory not defined in environment')
    OBS_DIR = DIR+'/obs'



# Forecast valid date/time
itime = ymdh
if hour%12 == 0:
    fhrs = np.arange(0,37,1)
else:
    fhrs = np.arange(0,13,1)
vtime_list = [ncepy.ndate(itime,x) for x in fhrs]


# Define prod and para strings
model_str = str(sys.argv[3])
if str.upper(model_str) == 'HRRR':
    prod_str = str.upper(model_str)+'v3'
    para_str = str.upper(model_str)+'v4'
elif str.upper(model_str) == 'RAP':
    prod_str = str.upper(model_str)+'v4'
    para_str = str.upper(model_str)+'v5'


###################################################
# Read in all variables and calculate differences #
###################################################
t1a = time.clock()

qpf_list_1 = []
qpf_list_2 = []
qpf_dif_list = []
qpf_list_3 = []
snow_list_1 = []
snow_list_2 = []
snow_dif_list = []
snow_list_3 = []

ref_rain_list_1 = []
ref_zr_list_1 = []
ref_ip_list_1 = []
ref_snow_list_1 = []
ref_mix_list_1 = []
ref_rain_list_2 = []
ref_zr_list_2 = []
ref_ip_list_2 = []
ref_snow_list_2 = []
ref_mix_list_2 = []

crain_list_1 = []
czr_list_1 = []
cip_list_1 = []
csnow_list_1 = []
cmix_list_1 = []
crain_list_2 = []
czr_list_2 = []
cip_list_2 = []
csnow_list_2 = []
cmix_list_2 = []

uh25_list_1 = []
uh25_list_2 = []
uh25_dif_list = []

torn_lats = []
torn_lons = []
torn_times = []
hail_lats = []
hail_lons = []
hail_times = []
wind_lats = []
wind_lons = []
wind_times = []

#for j in range(len(vtime_list)):
for j in range(len(vtime_list[0:2])):

    fhr = fhrs[j]
    fhour = str(fhr).zfill(2)
    print 'fhour '+fhour
    vtime = vtime_list[j]


    # Define the output files
    # Prod, para
    if str.upper(model_str) == 'HRRR':
        data1 = pygrib.open(PROD_DIR+'/hrrr.'+ymd+'.t'+cyc+'z.wrfprsf'+fhour+'.grib2')
     #  data2 = pygrib.open(PARA_DIR+'/hrrrx.'+ymd+'.t'+cyc+'z.wrfprsf'+fhour+'.grib2')
    elif str.upper(model_str) == 'RAP':
        data1 = pygrib.open(PROD_DIR+'/rap.'+ymd+'.t'+cyc+'z.awp130pgrbf'+fhour+'.grib2')
      # data2 = pygrib.open(PARA_DIR+'/rapx.'+ymd+'.t'+cyc+'z.awp130pgrbf'+fhour+'.grib2')
    data2 = data1

    # Stage IV
    data3 = pygrib.open(OBS_DIR+'/st4/ST4.'+vtime+'.01h')

    # Total precipitation
    qpf_1 = data1.select(name='Total Precipitation',lengthOfTimeRange=fhr)[0].values * 0.0393701
    qpf_2 = data2.select(name='Total Precipitation',lengthOfTimeRange=fhr)[0].values * 0.0393701
    qpf_dif = qpf_2 - qpf_1
    qpf_3 = data3.select(name='Total Precipitation')[0].values * 0.0393701 
    if j == 0:
        qpf_3 = qpf_3*0.

    qpf_list_1.append(qpf_1)
    qpf_list_2.append(qpf_2)
    qpf_dif_list.append(qpf_dif)
    qpf_list_3.append(qpf_3)

    # Snow depth
    snow_1 = data1.select(name='Snow depth')[0].values * 39.3701
    snow_2 = data2.select(name='Snow depth')[0].values * 39.3701
    snow_dif = snow_2 - snow_1

    snow_list_1.append(snow_1)
    snow_list_2.append(snow_2)
    snow_dif_list.append(snow_dif)

    # NOHRSC
    if int(vtime[8:10])%6 == 0:
        data4 = pygrib.open(OBS_DIR+'/nohrsc/sfav2_CONUS_6h_'+vtime+'_grid184.grb2')

        snow_3 = data4.select(name='Total snowfall')[0].values * 39.3701
        if j == 0:
            snow_3 = snow_3*0.
        snow_list_3.append(snow_3)


    # Updraft helicity
    if (fhr > 0):
        uh25_1 = data1.select(stepType='max',parameterName="199",topLevel=5000,bottomLevel=2000)[0].values
        uh25_2 = data2.select(stepType='max',parameterName="199",topLevel=5000,bottomLevel=2000)[0].values
        uh25_1[uh25_1 < 10] = 0
        uh25_2[uh25_2 < 10] = 0
    else:
        uh25_1 = data1.select(name='Total Precipitation',lengthOfTimeRange=fhr)[0].values * 0.
        uh25_2 = data2.select(name='Total Precipitation',lengthOfTimeRange=fhr)[0].values * 0.

    uh25_dif = uh25_2 - uh25_1
    uh25_list_1.append(uh25_1)
    uh25_list_2.append(uh25_2)
    uh25_dif_list.append(uh25_dif)


    # Precip type
    ref1 = np.asarray(data1.select(name='Maximum/Composite radar reflectivity',level=0)[0].values)
    crain1 = np.asarray(data1.select(name='Categorical rain',level=0)[0].values)
    czr1 = np.asarray(data1.select(name='Categorical freezing rain',level=0)[0].values)
    cip1 = np.asarray(data1.select(name='Categorical ice pellets',level=0)[0].values)
    csnow1 = np.asarray(data1.select(name='Categorical snow',level=0)[0].values)
    ptype1 = crain1 + czr1*3 + cip1*5 + csnow1*7
    cmix1 = np.copy(ptype1)
    crain1[ptype1!=1] = -1
    czr1[ptype1!=3] = -1
    cip1[ptype1!=5] = -1
    csnow1[ptype1!=7] = -1
    cmix1[ptype1==0] = -1
    cmix1[ptype1==1] = -1
    cmix1[ptype1==3] = -1
    cmix1[ptype1==5] = -1
    cmix1[ptype1==7] = -1

    ref_rain1 = np.copy(ref1)
    ref_zr1 = np.copy(ref1)
    ref_ip1 = np.copy(ref1)
    ref_snow1 = np.copy(ref1)
    ref_mix1 = np.copy(ref1)
    ref_rain1[ptype1!=1] = -1
    ref_zr1[ptype1!=3] = -1
    ref_ip1[ptype1!=5] = -1
    ref_snow1[ptype1!=7] = -1
    ref_mix1[ptype1==0] = -1
    ref_mix1[ptype1==1] = -1
    ref_mix1[ptype1==3] = -1
    ref_mix1[ptype1==5] = -1
    ref_mix1[ptype1==7] = -1
 #  cmix1 = np.copy(ref_mix1)
 #  cmix1[ref_mix1>0] = 1
    cmix1[cmix1>0] = 1

    ref2 = np.asarray(data2.select(name='Maximum/Composite radar reflectivity',level=0)[0].values)
    crain2 = np.asarray(data2.select(name='Categorical rain',level=0)[0].values)
    czr2 = np.asarray(data2.select(name='Categorical freezing rain',level=0)[0].values)
    cip2 = np.asarray(data2.select(name='Categorical ice pellets',level=0)[0].values)
    csnow2 = np.asarray(data2.select(name='Categorical snow',level=0)[0].values)
    ptype2 = crain2 + czr2*3 + cip2*5 + csnow2*7
    cmix2 = np.copy(ptype2)
    crain2[ptype2!=1] = -1
    czr2[ptype2!=3] = -1
    cip2[ptype2!=5] = -1
    csnow2[ptype2!=7] = -1
    cmix2[ptype2==0] = -1
    cmix2[ptype2==1] = -1
    cmix2[ptype2==3] = -1
    cmix2[ptype2==5] = -1
    cmix2[ptype2==7] = -1

    ref_rain2 = np.copy(ref2)
    ref_zr2 = np.copy(ref2)
    ref_ip2 = np.copy(ref2)
    ref_snow2 = np.copy(ref2)
    ref_mix2 = np.copy(ref2)
    ref_rain2[ptype2!=1] = -1
    ref_zr2[ptype2!=3] = -1
    ref_ip2[ptype2!=5] = -1
    ref_snow2[ptype2!=7] = -1
    ref_mix2[ptype2==0] = -1
    ref_mix2[ptype2==1] = -1
    ref_mix2[ptype2==3] = -1
    ref_mix2[ptype2==5] = -1
    ref_mix2[ptype2==7] = -1
 #  cmix2 = np.copy(ref_mix2)
    cmix2[ref_mix2>0] = 1

    crain_list_1.append(crain1)
    crain_list_2.append(crain2)
    czr_list_1.append(czr1)
    czr_list_2.append(czr2)
    cip_list_1.append(cip1)
    cip_list_2.append(cip2)
    csnow_list_1.append(csnow1)
    csnow_list_2.append(csnow2)
    cmix_list_1.append(cmix1)
    cmix_list_2.append(cmix2)

    ref_rain_list_1.append(ref_rain1)
    ref_rain_list_2.append(ref_rain2)
    ref_zr_list_1.append(ref_zr1)
    ref_zr_list_2.append(ref_zr2)
    ref_ip_list_1.append(ref_ip1)
    ref_ip_list_2.append(ref_ip2)
    ref_snow_list_1.append(ref_snow1)
    ref_snow_list_2.append(ref_snow2)
    ref_mix_list_1.append(ref_mix1)
    ref_mix_list_2.append(ref_mix2)



# Storm reports
if str.lower(case[0:3]) == 'jul' or str.lower(case[0:3]) == 'aug':
   retro_period = 'jul2018'
elif str.lower(case[0:3]) == 'feb' or str.lower(case[0:3]) == 'mar':
   retro_period = 'febmar2019'
elif str.lower(case[0:3]) == 'may':
   retro_period = 'may2019'

reports = ['torn','hail','wind']
h = 0
for hazard in reports:
    report_file = REP_DIR+'/'+retro_period+'_'+hazard+'.csv'

    with open(report_file,'r') as f:
        reader = csv.reader(f)
        for row in reader:
            rowtime = row[0]+row[1]+row[2]+row[3]   # in YYYYMMDDHHMM
            begin_val = vtime+'00'
            if rowtime >= vtime_list[0]+'00' and rowtime <= vtime_list[-1]+'00':
                if h == 0:
                    torn_times.append(rowtime)
                    torn_lats.append(float(row[4]))
                    torn_lons.append(float(row[5]))
                elif h == 1:
                    hail_times.append(rowtime)
                    hail_lats.append(float(row[4]))
                    hail_lons.append(float(row[5]))
                elif h == 2:
                    wind_times.append(rowtime)
                    wind_lats.append(float(row[4]))
                    wind_lons.append(float(row[5]))

    h += 1


# Get the lats and lons
lat,lon = data1.select(name='2 metre temperature')[0].latlons()
lat2,lon2 = data2.select(name='2 metre temperature')[0].latlons()
lat3,lon3 = data3.select(name='Total Precipitation')[0].latlons()
lat4,lon4 = data4.select(name='Total snowfall')[0].latlons()
Lon0 = data1[1]['LoVInDegrees']
Lat0 = data1[1]['LaDInDegrees']

t2a = time.clock()
t3a = round(t2a-t1a, 3)
print("%.3f seconds to read all messages") % t3a


# Specify plotting domains
splains_cases   = ['may01','may02','may05','may07','may17','may18','may20','may23','may29']
cplains_cases   = ['may04','may05','may06','may17','may21','may22','may23','may26','may27','may28']
nplains_cases   = ['may15']
midwest_cases   = ['may19','may22','may26','may27','may28','may29']
northeast_cases = ['may19','may23','may26','may28','may29']
southeast_cases = ['may04','may05','may13']

if str.upper(model_str) != 'HRRR-AK' and str.upper(model_str) != 'RAP-AK':
    domains = ['conus']

    if case in splains_cases:
        domains.extend(['splains'])
    elif case in nplain_cases:
        domains.extend(['nplains'])
    elif case in midwest_cases:
        domains.extend(['midwest'])
    elif case in northeast_cases:
        domains.extend(['northeast'])
    elif case in southeast_cases:
        domains.extend(['southeast'])

else:
    domains = ['alaska'] 


plots = [n for n in itertools.product(domains,fhrs[0:2])]

# colors for difference plots, only need to define once
difcolors = ['blue','#1874CD','dodgerblue','deepskyblue','turquoise','white','white','#EEEE00','#EEC900','darkorange','orangered','red']

# colors for ptype plots
snowhex=["#64B3E8", "#5197D7", "#3E7CC6", "#2B60B5", "#1945A4"]
rainhex=["#5EE240", "#4DC534", "#3DA828", "#2D8B1C", "#1D6F11"]
iphex=["#947EEC", "#7F62CB", "#6B47AB", "#562B8A", "#42106A"]
zrhex=["#E65956", "#DF4A48", "#D93B3A", "#D22C2C", "#CC1E1E"]
mixhex=["#E75FD5", "#D54DBB", "#C33BA2", "#B12989", "#A01870"]

########################################
#    START PLOTTING FOR EACH DOMAIN    #
########################################

def main():

  # Number of processes must coincide with the number of domains to plot
# pool = multiprocessing.Pool(len(fhrs[0:2]))
# pool.map(plot_all,fhrs[0:2])
  pool = multiprocessing.Pool(len(plots))
  print plots
  pool.map(plot_all,plots)

def plot_all(plot):

  thing = np.asarray(plot)
  dom = thing[0]
  fhr = int(thing[1])

  fhour = str(fhr).zfill(2)
  vtime = vtime_list[fhr]

  t1dom = time.clock()
  print('Working on '+dom+' for fhr '+fhour)

  # create figure and axes instances
  fig = plt.figure()
  gs = GridSpec(9,9,wspace=0.0,hspace=0.0)
  ax1 = fig.add_subplot(gs[0:4,0:4])
  ax2 = fig.add_subplot(gs[0:4,5:])
  ax3 = fig.add_subplot(gs[5:,0:4])
  ax4 = fig.add_subplot(gs[5:,5:])
  axes = [ax1, ax2, ax3, ax4]
  im = image.imread('/gpfs/hps3/emc/meso/save/Benjamin.Blake/python.raphrrr/noaa.png')
  par = 1

  # Map corners for each domain
  if dom == 'conus':
    llcrnrlon = -120.5
    llcrnrlat = 21.0 
    urcrnrlon = -64.5
    urcrnrlat = 49.0
    lat_0 = 35.4
    lon_0 = -97.6
    xscale=0.15
    yscale=0.2
  elif dom == 'splains':
    llcrnrlon = -105.0
    llcrnrlat = 28.5 
    urcrnrlon = -88.0
    urcrnrlat = 40.
    lat_0 = 33.5
    lon_0 = -96.5
    xscale=0.17
    yscale=0.18
  elif dom == 'nplains':
    llcrnrlon = -105.0
    llcrnrlat = 38.0 
    urcrnrlon = -88.0
    urcrnrlat = 49.0
    lat_0 = 33.5
    lon_0 = -96.5
    xscale=0.17
    yscale=0.18
  elif dom == 'midwest':
    llcrnrlon = -96.5
    llcrnrlat = 36.0 
    urcrnrlon = -79.5
    urcrnrlat = 47.5
    lat_0 = 33.5
    lon_0 = -88.0
    xscale=0.17
    yscale=0.18
  elif dom == 'northeast':
    llcrnrlon = -84.5
    llcrnrlat = 36.0 
    urcrnrlon = -64.5
    urcrnrlat = 47.5
    lat_0 = 33.5
    lon_0 = -74.5
    xscale=0.17
    yscale=0.18
  elif dom == 'southeast':
    llcrnrlon = -94.0
    llcrnrlat = 26.5 
    urcrnrlon = -75.0
    urcrnrlat = 38.0
    lat_0 = 33.5
    lon_0 = -84.5
    xscale=0.17
    yscale=0.18

  elif dom == 'BN':
    llcrnrlon = -75.75
    llcrnrlat = 40.0
    urcrnrlon = -69.5
    urcrnrlat = 43.0
    lat_0 = 41.0
    lon_0 = -74.6
    xscale=0.14
    yscale=0.19
  elif dom == 'CE':
    llcrnrlon = -103.0
    llcrnrlat = 32.5
    urcrnrlon = -88.5
    urcrnrlat = 41.5
    lat_0 = 35.0
    lon_0 = -97.0
    xscale=0.15
    yscale=0.18
  elif dom == 'CO':
    llcrnrlon = -110.5
    llcrnrlat = 35.0
    urcrnrlon = -100.5
    urcrnrlat = 42.0
    lat_0 = 38.0
    lon_0 = -105.0
    xscale=0.17
    yscale=0.18
  elif dom == 'LA':
    llcrnrlon = -121.0
    llcrnrlat = 32.0
    urcrnrlon = -114.0
    urcrnrlat = 37.0
    lat_0 = 34.0
    lon_0 = -114.0
    xscale=0.16
    yscale=0.18
  elif dom == 'MA':
    llcrnrlon = -82.0
    llcrnrlat = 36.5
    urcrnrlon = -73.5
    urcrnrlat = 42.0
    lat_0 = 37.5
    lon_0 = -80.0
    xscale=0.18
    yscale=0.18
  elif dom == 'NC':
    llcrnrlon = -111.0
    llcrnrlat = 39.0
    urcrnrlon = -93.5
    urcrnrlat = 49.0
    lat_0 = 44.5
    lon_0 = -102.0
    xscale=0.16
    yscale=0.18
  elif dom == 'NE':
    llcrnrlon = -80.0     
    llcrnrlat = 40.5
    urcrnrlon = -66.0
    urcrnrlat = 47.5
    lat_0 = 42.0
    lon_0 = -80.0
    xscale=0.16
    yscale=0.18
  elif dom == 'NW':
    llcrnrlon = -125.5     
    llcrnrlat = 40.5
    urcrnrlon = -109.0
    urcrnrlat = 49.5
    lat_0 = 44.0
    lon_0 = -116.0
    xscale=0.15
    yscale=0.18
  elif dom == 'OV':
    llcrnrlon = -91.5 
    llcrnrlat = 34.75
    urcrnrlon = -80.0
    urcrnrlat = 43.0
    lat_0 = 38.0
    lon_0 = -87.0          
    xscale=0.18
    yscale=0.17
  elif dom == 'SC':
    llcrnrlon = -108.0 
    llcrnrlat = 25.0
    urcrnrlon = -88.0
    urcrnrlat = 37.0
    lat_0 = 32.0
    lon_0 = -98.0      
    xscale=0.14
    yscale=0.18
  elif dom == 'SE':
    llcrnrlon = -91.5 
    llcrnrlat = 24.0
    urcrnrlon = -74.0
    urcrnrlat = 36.5
    lat_0 = 34.0
    lon_0 = -85.0
    xscale=0.17
    yscale=0.18
  elif dom == 'SF':
    llcrnrlon = -123.25 
    llcrnrlat = 37.25
    urcrnrlon = -121.25
    urcrnrlat = 38.5
    lat_0 = 37.5
    lon_0 = -121.0
    xscale=0.16
    yscale=0.19
  elif dom == 'SP':
    llcrnrlon = -125.0
    llcrnrlat = 45.0
    urcrnrlon = -119.5
    urcrnrlat = 49.2
    lat_0 = 46.0
    lon_0 = -120.0
    xscale=0.19
    yscale=0.18
  elif dom == 'SW':
    llcrnrlon = -125.0 
    llcrnrlat = 30.0
    urcrnrlon = -108.0
    urcrnrlat = 42.5
    lat_0 = 37.0
    lon_0 = -113.0
    xscale=0.17
    yscale=0.18
  elif dom == 'UM':
    llcrnrlon = -96.75 
    llcrnrlat = 39.75
    urcrnrlon = -81.0
    urcrnrlat = 49.0
    lat_0 = 44.0
    lon_0 = -91.5
    xscale=0.18
    yscale=0.18

  # Create basemap instance and set the dimensions
  for ax in axes:
    if dom == 'BN' or dom == 'LA' or dom == 'SF' or dom == 'SP':
      m = Basemap(ax=ax,projection='gnom',lat_0=lat_0,lon_0=lon_0,\
                  llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,\
                  llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,\
                  resolution='h')
    else:
      m = Basemap(ax=ax,projection='gnom',lat_0=lat_0,lon_0=lon_0,\
                  llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,\
                  llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,\
                  resolution='l')
    m.fillcontinents(color='LightGrey',zorder=0)
    m.drawcoastlines(linewidth=0.75)
    m.drawstates(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
##  parallels = np.arange(0.,90.,10.)
##  map.drawparallels(parallels,labels=[1,0,0,0],fontsize=6)
##  meridians = np.arange(180.,360.,10.)
##  map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=6)
    x,y = m(lon,lat)
    x2,y2 = m(lon2,lat2)
    x3,y3 = m(lon3,lat3)
    x4,y4 = m(lon4,lat4)
  
  # Map/figure has been set up here, save axes instances for use again later
    if par == 1:
      keep_ax_lst_1 = ax.get_children()[:]
    elif par == 2:
      keep_ax_lst_2 = ax.get_children()[:]
    elif par == 3:
      keep_ax_lst_3 = ax.get_children()[:]
    elif par == 4:
      keep_ax_lst_4 = ax.get_children()[:]

    par += 1
  par = 1


#################################
  # Plot Total QPF
#################################
  if (fhr >= 0):
    t1 = time.clock()
    print('Working on total qpf for '+dom)

    qpf_1 = qpf_list_1[fhr]
    qpf_2 = qpf_list_2[fhr]
    qpf_dif = qpf_dif_list[fhr]
    qpf_3 = np.sum(np.array(qpf_list_3[0:fhr+1]),axis=0)

    units = 'in'
    clevs = [0.01,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,4,5,7,10,15,20]
    clevsdif = [-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3]
    colorlist = ['chartreuse','limegreen','green','blue','dodgerblue','deepskyblue','cyan','mediumpurple','mediumorchid','darkmagenta','darkred','crimson','orangered','darkorange','goldenrod','gold','yellow']  

    for ax in axes:
      xmin, xmax = ax.get_xlim()
      ymin, ymax = ax.get_ylim()
      xmax = int(round(xmax))
      ymax = int(round(ymax))

      if par == 1:
        cs_1 = m.contourf(x,y,qpf_1,clevs,colors=colorlist,extend='max',ax=ax)
        cs_1.cmap.set_over('pink')
        cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05)
        cbar1.set_label(units,fontsize=6)
        cbar1.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,prod_str+' '+fhour+'-hr Accumulated Precipitation ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 2:
        cs_2 = m.contourf(x2,y2,qpf_2,clevs,colors=colorlist,extend='max',ax=ax)
        cs_2.cmap.set_over('pink')
        cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05)
        cbar2.set_label(units,fontsize=6)
        cbar2.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' '+fhour+'-hr Accumulated Precipitation ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 3:
        cs = m.contourf(x2,y2,qpf_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
        cs.cmap.set_under('darkblue')
        cs.cmap.set_over('darkred')
        cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
        cbar3.set_label(units,fontsize=6)
        cbar3.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' - '+prod_str+' '+fhour+'-hr Accumulated Precipitation ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))         
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 4:
        qpf_3[qpf_3 > 1000.0] = 0    # avoid plotting missing values with pink
        cs_3 = m.contourf(x3,y3,qpf_3,clevs,colors=colorlist,extend='max',ax=ax)
        cs_3.cmap.set_over('pink')
        cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05)
        cbar4.set_label(units,fontsize=6)
        cbar4.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,'Stage IV '+fhour+'-hr Accumulated Precipitation ('+units+') \n valid: '+itime+' to '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      par += 1
    par = 1

    compress_and_save('compareqpf_'+dom+'_f'+fhour+'.png')
#    plt.savefig('./compareqpf_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
    t2 = time.clock()
    t3 = round(t2-t1, 3)
    print('%.3f seconds to plot total qpf for: '+dom) % t3


#################################
  # Plot QPF3
#################################
  if (fhr >= 0):
    t1 = time.clock()
    print('Working on qpf3 for '+dom)

    # Clear off old plottables but keep all the map info
    clear_plotables(ax1,keep_ax_lst_1,fig)
    clear_plotables(ax2,keep_ax_lst_2,fig)
    clear_plotables(ax3,keep_ax_lst_3,fig)
    clear_plotables(ax4,keep_ax_lst_4,fig)

    if fhr > 3:
        qpf3_1 = qpf_1 - qpf_list_1[fhr-3]
        qpf3_2 = qpf_2 - qpf_list_2[fhr-3]
        qpf3_dif = qpf_dif - qpf_dif_list[fhr-3]
        qpf3_3 = qpf_3 - np.sum(np.array(qpf_list_3[0:fhr-2]),axis=0)
        accum_str = '03'
        vtime_m3 = vtime_list[fhr-3]
    else:
        qpf3_1 = qpf_1
        qpf3_2 = qpf_2
        qpf3_dif = qpf_dif
        qpf3_3 = qpf_3
        accum_str = fhour.zfill(1)
        vtime_m3 = itime

    units = 'in'
    clevs = [0.01,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,4,5,7,10,15,20]
    clevsdif = [-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5]
    colorlist = ['chartreuse','limegreen','green','blue','dodgerblue','deepskyblue','cyan','mediumpurple','mediumorchid','darkmagenta','darkred','crimson','orangered','darkorange','goldenrod','gold','yellow']  
   
    for ax in axes:
      xmin, xmax = ax.get_xlim()
      ymin, ymax = ax.get_ylim()
      xmax = int(round(xmax))
      ymax = int(round(ymax))

      if par == 1:
        cs_1 = m.contourf(x,y,qpf3_1,clevs,colors=colorlist,extend='max',ax=ax)
        cs_1.cmap.set_over('pink')
        cbar1.remove()
        cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05)
        cbar1.set_label(units,fontsize=6)
        cbar1.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,prod_str+' '+accum_str+'-hr Accumulated Precipitation ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 2:
        cs_2 = m.contourf(x2,y2,qpf3_2,clevs,colors=colorlist,extend='max',ax=ax)
        cs_2.cmap.set_over('pink')
        cbar2.remove()
        cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05)
        cbar2.set_label(units,fontsize=6)
        cbar2.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' '+accum_str+'-hr Accumulated Precipitation ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 3:
        cs = m.contourf(x2,y2,qpf3_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
        cs.cmap.set_under('darkblue')
        cs.cmap.set_over('darkred')
        cbar3.remove()
        cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
        cbar3.set_label(units,fontsize=6)
        cbar3.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' - '+prod_str+' '+accum_str+'-hr Accumulated Precipitation ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))         
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 4:
        qpf_3[qpf_3 > 1000.0] = 0    # avoid plotting missing values with pink
        cs_3 = m.contourf(x3,y3,qpf3_3,clevs,colors=colorlist,extend='max',ax=ax)
        cs_3.cmap.set_over('pink')
        cbar4.remove()
        cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05)
        cbar4.set_label(units,fontsize=6)
        cbar4.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,'Stage IV '+accum_str+'-hr Accumulated Precipitation ('+units+') \n valid: '+vtime_m3+' to '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      par += 1
    par = 1

    compress_and_save('compareqpf3_'+dom+'_f'+fhour+'.png')
#    plt.savefig('./compareqpf3_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
    t2 = time.clock()
    t3 = round(t2-t1, 3)
    print('%.3f seconds to plot qpf3 for: '+dom) % t3


#####################################
  # Plot total accumulated snow depth
#####################################
  if (fhr%6 == 0):
    t1 = time.clock()
    print('Working on total accumulated snow depth for '+dom)

    # Clear off old plottables but keep all the map info
    clear_plotables(ax1,keep_ax_lst_1,fig)
    clear_plotables(ax2,keep_ax_lst_2,fig)
    clear_plotables(ax3,keep_ax_lst_3,fig)
    clear_plotables(ax4,keep_ax_lst_4,fig)

    snow_1 = snow_list_1[fhr] - snow_list_1[0]
    snow_2 = snow_list_2[fhr] - snow_list_2[0]
    snow_dif = snow_dif_list[fhr] - snow_dif_list[0]
    snow_3 = np.sum(np.array(snow_list_3[0:(fhr/6)+1]),axis=0)

    units = 'in'
    clevs = [0.1,1,2,3,6,9,12,18,24,36,48]
    clevsdif = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]
    cm = ncepy.ncl_perc_11Lev()
    norm = matplotlib.colors.BoundaryNorm(clevs, cm.N)

    for ax in axes:
      xmin, xmax = ax.get_xlim()
      ymin, ymax = ax.get_ylim()
      xmax = int(round(xmax))
      ymax = int(round(ymax))

      if par == 1:
        cs_1 = m.contourf(x,y,snow_1,clevs,cmap=cm,norm=norm,extend='both',ax=ax)
        cs_1.cmap.set_under('white')
        cbar1.remove()
        cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
        cbar1.set_label(units,fontsize=6)
        cbar1.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,prod_str+' Accumulated Snow Depth ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 2:
        cs_2 = m.contourf(x2,y2,snow_2,clevs,cmap=cm,norm=norm,extend='both',ax=ax)
        cs_2.cmap.set_under('white')
        cbar2.remove()
        cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
        cbar2.set_label(units,fontsize=6)
        cbar2.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' Accumulated Snow Depth ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 3:
        cs = m.contourf(x2,y2,snow_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
        cs.cmap.set_under('darkblue')
        cs.cmap.set_over('darkred')
        cbar3.remove()
        cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
        cbar3.set_label(units,fontsize=6)
        cbar3.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' - '+prod_str+' Accumulated Snow Depth ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 4:
        snow_3[snow_3 > 1000.0] = 0    # avoid plotting missing values with pink
        cs_3 = m.contourf(x4,y4,snow_3,clevs,cmap=cm,norm=norm,extend='both',ax=ax)
        cs_3.cmap.set_under('white')
        cbar4.remove()
        cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=clevs)
        cbar4.set_label(units,fontsize=6)
        cbar4.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,'NOHRSC '+fhour+'-hr Accumulated Snow Depth ('+units+') \n valid: '+itime+' to '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      par += 1
    par = 1

    compress_and_save('comparesnow_'+dom+'_f'+fhour+'.png')
#   plt.savefig('./comparesnow_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
    t2 = time.clock()
    t3 = round(t2-t1, 3)
    print('%.3f seconds to plot accumulated snow depth for: '+dom) % t3


#################################
  # Plot 6-hr change in snow depth
#################################
#  if (fhr % 3 == 0) and (fhr >= 6):
    t1 = time.clock()
    print('Working on 6-hr change in snow depth for '+dom)

    # Clear off old plottables but keep all the map info
    clear_plotables(ax1,keep_ax_lst_1,fig)
    clear_plotables(ax2,keep_ax_lst_2,fig)
    clear_plotables(ax3,keep_ax_lst_3,fig)
    clear_plotables(ax4,keep_ax_lst_4,fig)

    if fhr == 0:
        vtime_m6 = vtime_list[fhr]
        snow6_1 = snow_list_1[fhr]*0.
        snow6_2 = snow_list_2[fhr]*0.
        snow6_dif = snow_dif_list[fhr]*0.
        snow6_3 = snow_list_3[fhr]*0.
    else:
        vtime_m6 = vtime_list[fhr-6]
        snow6_1 = snow_list_1[fhr] - snow_list_1[fhr-6]
        snow6_2 = snow_list_2[fhr] - snow_list_2[fhr-6]
        snow6_dif = snow_dif_list[fhr] - snow_dif_list[fhr-6]
        print(fhr, fhr/6, (fhr-6)/6)
        snow6_3 = snow_list_3[fhr/6] - snow_list_3[(fhr-6)/6]

    units = 'in'
    clevs = [-6,-4,-3,-2,-1,-0.5,0,0.5,1,2,3,4,6]
    clevsdif = [-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3]
    colorlist = ['blue','#1874CD','dodgerblue','deepskyblue','turquoise','white','white','#EEEE00','#EEC900','darkorange','orangered','red']

    for ax in axes:
      xmin, xmax = ax.get_xlim()
      ymin, ymax = ax.get_ylim()
      xmax = int(round(xmax))
      ymax = int(round(ymax))

      if par == 1:
        cs_1 = m.contourf(x,y,snow6_1,clevs,colors=colorlist,extend='both',ax=ax)
        cs_1.cmap.set_under('darkblue')
        cs_1.cmap.set_over('darkred')
        cbar1.remove()
        cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
        cbar1.set_label(units,fontsize=6)
        cbar1.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,prod_str+' 6-hr Change in Snow Depth ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 2:
        cs_2 = m.contourf(x2,y2,snow6_2,clevs,colors=colorlist,extend='both',ax=ax)
        cs_2.cmap.set_under('darkblue')
        cs_2.cmap.set_over('darkred')
        cbar2.remove()
        cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
        cbar2.set_label(units,fontsize=6)
        cbar2.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' 6-hr Change in Snow Depth ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 3:
        cs = m.contourf(x2,y2,snow6_dif,clevsdif,colors=colorlist,extend='both',ax=ax)
        cs.cmap.set_under('darkblue')
        cs.cmap.set_over('darkred')
        cbar3.remove()
        cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
        cbar3.set_label(units,fontsize=6)
        cbar3.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' - '+prod_str+' 6-hr Change in Snow Depth ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 4:
        cs_3 = m.contourf(x4,y4,snow6_3,clevs,colors=colorlist,extend='both',ax=ax)
        cs_3.cmap.set_under('darkblue')
        cs_3.cmap.set_over('darkred')
        cbar4.remove()
        cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=clevs)
        cbar4.set_label(units,fontsize=6)
        cbar4.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,'NOHRSC 6-hr Change in Snow Depth ('+units+') \n valid: '+vtime_m6+' to '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      par += 1
    par = 1

    compress_and_save('comparesnow6_'+dom+'_f'+fhour+'.png')
#   plt.savefig('./comparesnow6_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
    t2 = time.clock()
    t3 = round(t2-t1, 3)
    print('%.3f seconds to plot 6-hr change in snow depth for: '+dom) % t3





#################################
  # Plot Max Hourly 2-5 km UH
#################################
  if (fhr >= 0):
    t1 = time.clock()
    print('Working on Max Hourly 2-5 km UH for '+dom)

    # Clear off old plottables but keep all the map info
    clear_plotables(ax1,keep_ax_lst_1,fig)
    clear_plotables(ax2,keep_ax_lst_2,fig)
    clear_plotables(ax3,keep_ax_lst_3,fig)
    clear_plotables(ax4,keep_ax_lst_4,fig)

    uh25_1 = uh25_list_1[fhr]
    uh25_2 = uh25_list_2[fhr]
    uh25_dif = uh25_dif_list[fhr]

    units = 'm${^2}$ s$^{-2}$'
    units = '$\mathregular{m^{2}}$ $\mathregular{s^{-2}}$'
    clevs = [25,50,75,100,150,200,250,300]
    clevsdif = [-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60]
    colorlist = ['turquoise','dodgerblue','lime','limegreen','yellow','darkorange','red','firebrick','fuchsia']

    if dom == 'conus':
        markersize = 5
    else:
        markersize = 7

    for ax in axes:
      xmin, xmax = ax.get_xlim()
      ymin, ymax = ax.get_ylim()
      xmax = int(round(xmax))
      ymax = int(round(ymax))

      if par == 1:
        cs_1 = m.contourf(x,y,uh25_1,clevs,colors=colorlist,extend='max',ax=ax)
        cs_1.cmap.set_over('fuchsia')
        cbar1.remove()
        cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05)
        cbar1.set_label(units,fontsize=6)
        cbar1.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,prod_str+' 1-h Max 2-5 km Updraft Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 2:
        cs_2 = m.contourf(x2,y2,uh25_2,clevs,colors=colorlist,extend='max',ax=ax)
        cs_2.cmap.set_over('fuchsia')
        cbar2.remove()
        cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05)
        cbar2.set_label(units,fontsize=6)
        cbar2.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' 1-h Max 2-5 km Updraft Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 3:
        cs = m.contourf(x2,y2,uh25_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
        cs.cmap.set_under('darkblue')
        cs.cmap.set_over('darkred')
        cbar3.remove()
        cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
        cbar3.set_label(units,fontsize=6)
        cbar3.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' - '+prod_str+' 1-h Max 2-5 km Updraft Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 4:
        if fhr == 0:
            vtime_m1 = vtime
        else:
            vtime_m1 = vtime_list[fhr-1]

        # Determine tornado reports to plot
        plot_torn_lats = []
        plot_torn_lons = []
        for i in xrange(len(torn_times)):
            if torn_times[i] >= vtime_m1 and torn_times[i] <= vtime:
                plot_torn_lats.append(torn_lats[i])
                plot_torn_lons.append(torn_lons[i])

        # Determine hail reports to plot
        plot_hail_lats = []
        plot_hail_lons = []
        for i in xrange(len(hail_times)):
            if hail_times[i] >= vtime_m1 and hail_times[i] <= vtime:
                plot_hail_lats.append(hail_lats[i])
                plot_hail_lons.append(hail_lons[i])

        # Determine wind reports to plot
        plot_wind_lats = []
        plot_wind_lons = []
        for i in xrange(len(wind_times)):
            if wind_times[i] >= vtime_m1 and wind_times[i] <= vtime:
                plot_wind_lats.append(wind_lats[i])
                plot_wind_lons.append(wind_lons[i])


        xt, yt = m(plot_torn_lons, plot_torn_lats)
        xh, yh = m(plot_hail_lons, plot_hail_lats)
        xw, yw = m(plot_wind_lons, plot_wind_lats)

        wind_reps = m.scatter(xw, yw, marker='s', color='b', s=markersize, ax=ax) 
        hail_reps = m.scatter(xh, yh, marker='o', color='g', s=markersize, ax=ax) 
        torn_reps = m.scatter(xt, yt, marker='v', color='r', s=markersize, ax=ax) 

        cbar4.remove()
        ax.text(.5,1.03,'Local Storm Reports valid: '+vtime_m1+' to '+vtime+'\n Tornado (red), Hail (green), Wind (blue)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      par += 1
    par = 1

    compress_and_save('compareuh25_'+dom+'_f'+fhour+'.png')
    t2 = time.clock()
    t3 = round(t2-t1, 3)
    print('%.3f seconds to plot Max Hourly 2-5 km UH for: '+dom) % t3


#################################
  # Plot Run-Total 2-5 km UH
#################################
  if (fhr >= 0):
    t1 = time.clock()
    print('Working on Run-Total 2-5 km UH for '+dom)

    # Clear off old plottables but keep all the map info
    clear_plotables(ax1,keep_ax_lst_1,fig)
    clear_plotables(ax2,keep_ax_lst_2,fig)
    clear_plotables(ax3,keep_ax_lst_3,fig)
    clear_plotables(ax4,keep_ax_lst_4,fig)

    uh25_1 = np.amax(np.array(uh25_list_1[0:fhr+1]),axis=0)
    uh25_2 = np.amax(np.array(uh25_list_2[0:fhr+1]),axis=0)
    uh25_dif = np.amax(np.array(uh25_dif_list[0:fhr+1]),axis=0)

    units = '$\mathregular{m^{2}}$ $\mathregular{s^{-2}}$'
    clevs = [25,50,75,100,150,200,250,300]
    clevsdif = [-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60]
    colorlist = ['turquoise','dodgerblue','lime','limegreen','yellow','darkorange','red','firebrick','fuchsia']
  # colorlist = ['white','skyblue','mediumblue','green','orchid','firebrick','#EEC900','DarkViolet']

    for ax in axes:
      xmin, xmax = ax.get_xlim()
      ymin, ymax = ax.get_ylim()
      xmax = int(round(xmax))
      ymax = int(round(ymax))

      if par == 1:
        cs_1 = m.contourf(x,y,uh25_1,clevs,colors=colorlist,extend='max',ax=ax)
        cs_1.cmap.set_over('fuchsia')
        cbar1.remove()
        cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05)
        cbar1.set_label(units,fontsize=6)
        cbar1.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,prod_str+' '+fhour+'-h Max 2-5 km Updraft Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 2:
        cs_2 = m.contourf(x2,y2,uh25_2,clevs,colors=colorlist,extend='max',ax=ax)
        cs_2.cmap.set_over('fuchsia')
        cbar2.remove()
        cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05)
        cbar2.set_label(units,fontsize=6)
        cbar2.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' '+fhour+'-h Max 2-5 km Updraft Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 3:
        cs = m.contourf(x2,y2,uh25_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
        cs.cmap.set_under('darkblue')
        cs.cmap.set_over('darkred')
        cbar3.remove()
        cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
        cbar3.set_label(units,fontsize=6)
        cbar3.ax.tick_params(labelsize=6)
        ax.text(.5,1.03,para_str+' - '+prod_str+' '+fhour+'-h Max 2-5 km Updraft Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 4:

        # Determine tornado reports to plot
        plot_torn_lats = []
        plot_torn_lons = []
        for i in xrange(len(torn_times)):
            if torn_times[i] >= itime and torn_times[i] <= vtime:
                plot_torn_lats.append(torn_lats[i])
                plot_torn_lons.append(torn_lons[i])

        # Determine hail reports to plot
        plot_hail_lats = []
        plot_hail_lons = []
        for i in xrange(len(hail_times)):
            if hail_times[i] >= itime and hail_times[i] <= vtime:
                plot_hail_lats.append(hail_lats[i])
                plot_hail_lons.append(hail_lons[i])

        # Determine wind reports to plot
        plot_wind_lats = []
        plot_wind_lons = []
        for i in xrange(len(wind_times)):
            if wind_times[i] >= itime and wind_times[i] <= vtime:
                plot_wind_lats.append(wind_lats[i])
                plot_wind_lons.append(wind_lons[i])


        xt, yt = m(plot_torn_lons, plot_torn_lats)
        xh, yh = m(plot_hail_lons, plot_hail_lats)
        xw, yw = m(plot_wind_lons, plot_wind_lats)

        wind_reps = m.scatter(xw, yw, marker='s', color='b', s=markersize, ax=ax) 
        hail_reps = m.scatter(xh, yh, marker='o', color='g', s=markersize, ax=ax) 
        torn_reps = m.scatter(xt, yt, marker='v', color='r', s=markersize, ax=ax) 

        ax.text(.5,1.03,'Local Storm Reports valid: '+itime+' to '+vtime+'\n Tornado (red), Hail (green), Wind (blue)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      par += 1
    par = 1

    compress_and_save('compareuh25_accum_'+dom+'_f'+fhour+'.png')
    t2 = time.clock()
    t3 = round(t2-t1, 3)
    print('%.3f seconds to plot Run-Total 2-5 km UH for: '+dom) % t3


#################################
  # Plot precip type
#################################
  t1 = time.clock()
  print('Working on ptype for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  crain1 = crain_list_1[fhr]
  czr1 = czr_list_1[fhr]
  cip1 = cip_list_1[fhr]
  csnow1 = csnow_list_1[fhr]
  cmix1 = cmix_list_1[fhr]

  crain2 = crain_list_2[fhr]
  czr2 = czr_list_2[fhr]
  cip2 = cip_list_2[fhr]
  csnow2 = csnow_list_2[fhr]
  cmix2 = cmix_list_2[fhr]

  ref_rain1 = ref_rain_list_1[fhr]
  ref_zr1 = ref_zr_list_1[fhr]
  ref_ip1 = ref_ip_list_1[fhr]
  ref_snow1 = ref_snow_list_1[fhr]
  ref_mix1 = ref_mix_list_1[fhr]

  ref_rain2 = ref_rain_list_2[fhr]
  ref_zr2 = ref_zr_list_2[fhr]
  ref_ip2 = ref_ip_list_2[fhr]
  ref_snow2 = ref_snow_list_2[fhr]
  ref_mix2 = ref_mix_list_2[fhr]

  cbar1.remove()
  cbar2.remove()
  cbar3.remove()

  clevs = [0,10,20,30,40]
  clevs2 = [0.5,1.5]

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      csrain=m.contourf(x,y,ref_rain1,clevs,colors=rainhex,extend='max',ax=ax)
      csmix=m.contourf(x,y,ref_mix1,clevs,colors=mixhex,extend='max',ax=ax)
      cssnow=m.contourf(x,y,ref_snow1,clevs,colors=snowhex,extend='max',ax=ax)
      cssleet=m.contourf(x,y,ref_ip1,clevs,colors=iphex,extend='max',ax=ax)
      csfrzra=m.contourf(x,y,ref_zr1,clevs,colors=zrhex,extend='max',ax=ax)
      ax.text(.5,1.03,prod_str+' Composite Reflectivity by Precip Type \n initialized: '+itime +' valid: '+ vtime + ' (F'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      csrain=m.contourf(x,y,ref_rain2,clevs,colors=rainhex,extend='max',ax=ax)
      csmix=m.contourf(x,y,ref_mix2,clevs,colors=mixhex,extend='max',ax=ax)
      cssnow=m.contourf(x,y,ref_snow2,clevs,colors=snowhex,extend='max',ax=ax)
      cssleet=m.contourf(x,y,ref_ip2,clevs,colors=iphex,extend='max',ax=ax)
      csfrzra=m.contourf(x,y,ref_zr2,clevs,colors=zrhex,extend='max',ax=ax)
      ax.text(.5,1.03,para_str+' Composite Reflectivity by Precip Ttype \n initialized: '+itime +' valid: '+ vtime + ' (F'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      caxrain=fig.add_axes([.09,.53,.1,.03])
      cbrain=fig.colorbar(csrain,cax=caxrain,ticks=clevs,orientation='horizontal')
      cbrain.ax.tick_params(labelsize=5)
      cbrain.ax.set_xticklabels(['Light Rain','','','','Heavy Rain'])

      caxsnow=fig.add_axes([.27,.53,.1,.03])
      cbsnow=fig.colorbar(cssnow,cax=caxsnow,ticks=clevs,orientation='horizontal')
      cbsnow.ax.tick_params(labelsize=5)
      cbsnow.ax.set_xticklabels(['Light Snow','','','','Heavy Snow'])

      caxsleet=fig.add_axes([.45,.53,.1,.03])
      cbsleet=fig.colorbar(cssleet,cax=caxsleet,ticks=clevs,orientation='horizontal')
      cbsleet.ax.tick_params(labelsize=5)
      cbsleet.ax.set_xticklabels(['Light Sleet','','','','Heavy Sleet'])

      caxfrzra=fig.add_axes([.63,.53,.1,.03])
      cbfrzra=fig.colorbar(csfrzra,cax=caxfrzra,ticks=clevs,orientation='horizontal')
      cbfrzra.ax.tick_params(labelsize=5)
      cbfrzra.ax.set_xticklabels(['Light Freezing Rain','','','','Heavy Freezing Rain'])

      caxmix=fig.add_axes([.81,.53,.1,.03])
      cbmix=fig.colorbar(csmix,cax=caxmix,ticks=clevs,orientation='horizontal')
      cbmix.ax.tick_params(labelsize=5)
      cbmix.ax.set_xticklabels(['Light Mix','','','','Heavy Mix'])

    elif par == 3:
      csrain=m.contourf(x,y,crain1,clevs2,colors=rainhex[-3],ax=ax)
      csmix=m.contourf(x,y,cmix1,clevs2,colors=mixhex[-3],ax=ax)
      cssnow=m.contourf(x,y,csnow1,clevs2,colors=snowhex[-3],ax=ax)
      cssleet=m.contourf(x,y,cip1,clevs2,colors=iphex[-3],ax=ax)
      csfrzra=m.contourf(x,y,czr1,clevs2,colors=zrhex[-3],ax=ax)
      ax.text(.5,1.03,prod_str+' Categorical Precipitation Type \n initialized: '+itime +' valid: '+ vtime + ' (F'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      csrain=m.contourf(x,y,crain2,clevs2,colors=rainhex[-3],ax=ax)
      csmix=m.contourf(x,y,cmix2,clevs2,colors=mixhex[-3],ax=ax)
      cssnow=m.contourf(x,y,csnow2,clevs2,colors=snowhex[-3],ax=ax)
      cssleet=m.contourf(x,y,cip2,clevs2,colors=iphex[-3],ax=ax)
      csfrzra=m.contourf(x,y,czr2,clevs2,colors=zrhex[-3],ax=ax)
      ax.text(.5,1.03,para_str+' Categorical Precipitation Type \n initialized: '+itime +' valid: '+ vtime + ' (F'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      caxrain2=fig.add_axes([.09,.08,.1,.03])
      cbrain2=fig.colorbar(csrain,cax=caxrain2,ticks=[1],orientation='horizontal')
      cbrain2.ax.tick_params(labelsize=6)
      cbrain2.ax.set_xticklabels(['Rain'])

      caxsnow2=fig.add_axes([.27,.08,.1,.03])
      cbsnow2=fig.colorbar(cssnow,cax=caxsnow2,ticks=[1],orientation='horizontal')
      cbsnow2.ax.tick_params(labelsize=6)
      cbsnow2.ax.set_xticklabels(['Snow'])

      caxsleet2=fig.add_axes([.45,.08,.1,.03])
      cbsleet2=fig.colorbar(cssleet,cax=caxsleet2,ticks=[1],orientation='horizontal')
      cbsleet2.ax.tick_params(labelsize=6)
      cbsleet2.ax.set_xticklabels(['Sleet'])

      caxfrzra2=fig.add_axes([.63,.08,.1,.03])
      cbfrzra2=fig.colorbar(csfrzra,cax=caxfrzra2,ticks=[1],orientation='horizontal')
      cbfrzra2.ax.tick_params(labelsize=6)
      cbfrzra2.ax.set_xticklabels(['Freezing Rain'])

      caxmix2=fig.add_axes([.81,.08,.1,.03])
      cbmix2=fig.colorbar(csmix,cax=caxmix2,ticks=[1],orientation='horizontal')
      cbmix2.ax.tick_params(labelsize=6)
      cbmix2.ax.set_xticklabels(['Mix'])

    par += 1
  par = 1

  compress_and_save('compareptype_'+dom+'_f'+fhour+'.png')
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot ptype for: '+dom) % t3








######################################################

  t3dom = round(t2-t1dom, 3)
  print("%.3f seconds to plot all variables for: "+dom) % t3dom
  plt.clf()

######################################################

main()
#plot_all('conus')
#plot_all(6)




