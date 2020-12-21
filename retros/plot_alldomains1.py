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
import time,os,sys,multiprocessing
import ncepy
import dawsonpy
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

def cmap_t2m():
 # Create colormap for 2-m temperature
 # Modified version of the ncl_t2m colormap from Jacob's ncepy code
    r=np.array([255,128,0,  70, 51, 0,  255,0, 0,  51, 255,255,255,255,255,171,128,128,36,162,255])
    g=np.array([0,  0,  0,  70, 102,162,255,92,128,185,255,214,153,102,0,  0,  0,  68, 36,162,255])
    b=np.array([255,128,128,255,255,255,255,0, 0,  102,0,  112,0,  0,  0,  56, 0,  68, 36,162,255])
    xsize=np.arange(np.size(r))
    r = r/255.
    g = g/255.
    b = b/255.
    red = []
    green = []
    blue = []
    for i in range(len(xsize)):
        xNorm=np.float(i)/(np.float(np.size(r))-1.0)
        red.append([xNorm,r[i],r[i]])
        green.append([xNorm,g[i],g[i]])
        blue.append([xNorm,b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    cmap_t2m_coltbl = matplotlib.colors.LinearSegmentedColormap('CMAP_T2M_COLTBL',colorDict)
    return cmap_t2m_coltbl


def cmap_t850():
 # Create colormap for 850-mb equivalent potential temperature
    r=np.array([255,128,0,  70, 51, 0,  0,  0, 51, 255,255,255,255,255,171,128,128,96,201])
    g=np.array([0,  0,  0,  70, 102,162,225,92,153,255,214,153,102,0,  0,  0,  68, 96,201])
    b=np.array([255,128,128,255,255,255,162,0, 102,0,  112,0,  0,  0,  56, 0,  68, 96,201])
    xsize=np.arange(np.size(r))
    r = r/255.
    g = g/255.
    b = b/255.
    red = []
    green = []
    blue = []
    for i in range(len(xsize)):
        xNorm=np.float(i)/(np.float(np.size(r))-1.0)
        red.append([xNorm,r[i],r[i]])
        green.append([xNorm,g[i],g[i]])
        blue.append([xNorm,b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    cmap_t850_coltbl = matplotlib.colors.LinearSegmentedColormap('CMAP_T850_COLTBL',colorDict)
    return cmap_t850_coltbl


def cmap_terra():
 # Create colormap for terrain height
 # Emerald green to light green to tan to gold to dark red to brown to light brown to white
    r=np.array([0,  152,212,188,127,119,186])
    g=np.array([128,201,208,148,34, 83, 186])
    b=np.array([64, 152,140,0,  34, 64, 186])
    xsize=np.arange(np.size(r))
    r = r/255.
    g = g/255.
    b = b/255.
    red = []
    green = []
    blue = []
    for i in range(len(xsize)):
        xNorm=np.float(i)/(np.float(np.size(r))-1.0)
        red.append([xNorm,r[i],r[i]])
        green.append([xNorm,g[i],g[i]])
        blue.append([xNorm,b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    cmap_terra_coltbl = matplotlib.colors.LinearSegmentedColormap('CMAP_TERRA_COLTBL',colorDict)
    cmap_terra_coltbl.set_over(color='#E0EEE0')
    return cmap_terra_coltbl


def extrema(mat,mode='wrap',window=100):
    # find the indices of local extrema (max only) in the input array.
    mx = ndimage.filters.maximum_filter(mat,size=window,mode=mode)
    # (mat == mx) true if pixel is equal to the local max
    return np.nonzero(mat == mx)

#-------------------------------------------------------#

# Necessary to generate figs when not running an Xserver (e.g. via PBS)
plt.switch_backend('agg')

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

fhr = int(sys.argv[3])
fhour = str(fhr).zfill(2)
print 'fhour '+fhour

# Forecast valid date/time
itime = ymdh
vtime = ncepy.ndate(itime,int(fhr))


# Get machine and head directory
machine, hostname = dawsonpy.get_machine()

if machine == 'WCOSS':
    DIR = '/gpfs/'+hostname[0]+'p2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
elif machine == 'WCOSS_C':
    DIR = '/gpfs/hps3/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
elif machine == 'WCOSS_DELL_P3':
    DIR = '/gpfs/dell2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros/'
elif machine == 'HERA':
    raise NameError, 'Need to modify plot_alldomains2.py script to define data/graphx head directory'


# Make sure prod directory doesn't fail
try:
    PROD_DIR = os.environ['PROD_DIR']
except KeyError:
    print('Prod directory not defined in environment')
    PROD_DIR = DIR+'/'+case+'prod'

# Make sure para directory doesn't fail
try:
    PARA_DIR = os.environ['PARA_DIR']
except KeyError:
    print('Para directory not defined in environment')
    PARA_DIR = DIR+'/'+case+'para'

# Make sure obs directory doesn't fail
try:
    OBS_DIR = os.environ['OBS_DIR']
except KeyError:
    print('Obs directory not defined in environment')
    OBS_DIR = DIR+'/obs'



# Define the output files
# Prod, para
model_str = str(sys.argv[4])
if str.upper(model_str) == 'HRRR':
    prod_str = str.upper(model_str)+'v3'
    para_str = str.upper(model_str)+'v4'
    data1 = pygrib.open(PROD_DIR+'/hrrr.'+ymd+'.t'+cyc+'z.wrfprsf'+fhour+'.grib2')
   #data2 = pygrib.open(PARA_DIR+'/hrrrx.'+ymd+'.t'+cyc+'z.wrfprsf'+fhour+'.grib2')
elif str.upper(model_str) == 'RAP':
    prod_str = str.upper(model_str)+'v4'
    para_str = str.upper(model_str)+'v5'
    data1 = pygrib.open(PROD_DIR+'/rap.'+ymd+'.t'+cyc+'z.awp130pgrbf'+fhour+'.grib2')
   #data2 = pygrib.open(PARA_DIR+'/rapx.'+ymd+'.t'+cyc+'z.awp130pgrbf'+fhour+'.grib2')
data2 = data1

# RAP, URMA, refc obs
data3 = pygrib.open(OBS_DIR+'/rap/rap.'+vtime[0:8]+'.t'+vtime[8:10]+'z.awp130pgrbf00.grib2')
data4 = pygrib.open(OBS_DIR+'/urma/urma2p5.'+vtime[0:8]+'.t'+vtime[8:10]+'z.2dvaranl_ndfd.grb2')
data5 = pygrib.open(OBS_DIR+'/refl/refd3d.'+vtime[0:8]+'.t'+vtime[8:10]+'z.grb2f00')

# Get the lats and lons
lat,lon = data1.select(name='2 metre temperature')[0].latlons()
lat2,lon2 = data2.select(name='2 metre temperature')[0].latlons()
lat3,lon3 = data3.select(name='2 metre temperature')[0].latlons()
lat4,lon4 = data4.select(name='2 metre temperature')[0].latlons()
lat5,lon5 = data5.select(name='Maximum/Composite radar reflectivity')[0].latlons()
Lon0 = data1[1]['LoVInDegrees']
Lat0 = data1[1]['LaDInDegrees']
Lon0_3 = data3[1]['LoVInDegrees']
Lat0_3 = data3[1]['LaDInDegrees']
Lon0_4 = data4[1]['LoVInDegrees']
Lat0_4 = data4[1]['LaDInDegrees']


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


###################################################
# Read in all variables and calculate differences #
###################################################
t1a = time.clock()


# Sea level pressure
slp_1 = data1.select(name='MSLP (MAPS System Reduction)')[0].values * 0.01
slpsmooth1 = ndimage.filters.gaussian_filter(slp_1, 13.78)
slp_2 = data2.select(name='MSLP (MAPS System Reduction)')[0].values * 0.01
slpsmooth2 = ndimage.filters.gaussian_filter(slp_2, 13.78)
slp_dif = slp_2 - slp_1
slp_3 = data3.select(name='MSLP (MAPS System Reduction)')[0].values * 0.01
slpsmooth3 = ndimage.filters.gaussian_filter(slp_3, 13.78)

# 2-m temperature
tmp2m_1 = data1.select(name='2 metre temperature')[0].values
tmp2m_1 = (tmp2m_1 - 273.15)*1.8 + 32.0
tmp2m_2 = data2.select(name='2 metre temperature')[0].values
tmp2m_2 = (tmp2m_2 - 273.15)*1.8 + 32.0
tmp2m_dif = tmp2m_2 - tmp2m_1
tmp2m_4 = data4.select(name='2 metre temperature')[0].values
tmp2m_4 = (tmp2m_4 - 273.15)*1.8 + 32.0

# Surface temperature
tmpsfc_1 = data1.select(name='Temperature',typeOfLevel='surface')[0].values
tmpsfc_1 = (tmpsfc_1 - 273.15)*1.8 + 32.0
tmpsfc_2 = data2.select(name='Temperature',typeOfLevel='surface')[0].values
tmpsfc_2 = (tmpsfc_2 - 273.15)*1.8 + 32.0
tmpsfc_dif = tmpsfc_2 - tmpsfc_1
tmpsfc_3 = data3.select(name='Temperature',typeOfLevel='surface')[0].values
tmpsfc_3 = (tmpsfc_3 - 273.15)*1.8 + 32.0

# 2-m dew point temperature
dew2m_1 = data1.select(name='2 metre dewpoint temperature')[0].values
dew2m_1 = (dew2m_1 - 273.15)*1.8 + 32.0
dew2m_2 = data2.select(name='2 metre dewpoint temperature')[0].values
dew2m_2 = (dew2m_2 - 273.15)*1.8 + 32.0
dew2m_dif = dew2m_2 - dew2m_1
dew2m_4 = data4.select(name='2 metre dewpoint temperature')[0].values
dew2m_4 = (dew2m_4 - 273.15)*1.8 + 32.0

# 10-m wind speed
uwind_1 = data1.select(name='10 metre U wind component')[0].values * 1.94384
uwind_2 = data2.select(name='10 metre U wind component')[0].values * 1.94384
uwind_4 = data4.select(name='10 metre U wind component')[0].values * 1.94384
vwind_1 = data1.select(name='10 metre V wind component')[0].values * 1.94384
vwind_2 = data2.select(name='10 metre V wind component')[0].values * 1.94384
vwind_4 = data4.select(name='10 metre V wind component')[0].values * 1.94384
# Rotate winds from grid relative to Earth relative
uwind_1, vwind_1 = ncepy.rotate_wind(Lat0,Lon0,lon,uwind_1,vwind_1,'lcc',inverse=False)
uwind_2, vwind_2 = ncepy.rotate_wind(Lat0,Lon0,lon2,uwind_2,vwind_2,'lcc',inverse=False)
uwind_4, vwind_4 = ncepy.rotate_wind(Lat0_4,Lon0_4,lon4,uwind_4,vwind_4,'lcc',inverse=False)
wspd10m_1 = np.sqrt(uwind_1**2 + vwind_1**2)
wspd10m_2 = np.sqrt(uwind_2**2 + vwind_2**2)
wspd10m_dif = wspd10m_2 - wspd10m_1
wspd10m_4 = np.sqrt(uwind_4**2 + vwind_4**2)


# Most unstable CAPE
mucape_1 = data1.select(name='Convective available potential energy',topLevel=18000)[0].values
mucape_2 = data2.select(name='Convective available potential energy',topLevel=18000)[0].values
mucape_dif = mucape_2 - mucape_1
mucape_3 = data3.select(name='Convective available potential energy',topLevel=18000)[0].values

# Most Unstable CIN
mucin_1 = data1.select(name='Convective inhibition',topLevel=18000)[0].values
mucin_2 = data2.select(name='Convective inhibition',topLevel=18000)[0].values
mucin_dif = mucin_2 - mucin_1
mucin_3 = data3.select(name='Convective inhibition',topLevel=18000)[0].values

# Surface-based CAPE
cape_1 = data1.select(name='Convective available potential energy',typeOfLevel='surface')[0].values
cape_2 = data2.select(name='Convective available potential energy',typeOfLevel='surface')[0].values
cape_dif = cape_2 - cape_1
cape_3 = data3.select(name='Convective available potential energy',typeOfLevel='surface')[0].values

# Surface-based CIN
sfcin_1 = data1.select(name='Convective inhibition',typeOfLevel='surface')[0].values
sfcin_2 = data2.select(name='Convective inhibition',typeOfLevel='surface')[0].values
sfcin_dif = sfcin_2 - sfcin_1
sfcin_3 = data3.select(name='Convective inhibition',typeOfLevel='surface')[0].values

# Mixed Layer CAPE
mlcape_1 = data1.select(name='Convective available potential energy',topLevel=9000)[0].values
mlcape_2 = data2.select(name='Convective available potential energy',topLevel=9000)[0].values
mlcape_dif = mlcape_2 - mlcape_1
mlcape_3 = data3.select(name='Convective available potential energy',topLevel=9000)[0].values

# Mixed Layer CIN
mlcin_1 = data1.select(name='Convective inhibition',topLevel=9000)[0].values
mlcin_2 = data2.select(name='Convective inhibition',topLevel=9000)[0].values
mlcin_dif = mlcin_2 - mlcin_1
mlcin_3 = data3.select(name='Convective inhibition',topLevel=9000)[0].values

# 850-mb temperature and equivalent potential temperature
t850_1 = data1.select(name='Temperature',level=850)[0].values
rh850_1 = data1.select(name='Relative humidity',level=850)[0].values
es_1 = 6.11*np.exp((53.49*np.ones(t850_1.shape))-((6808.*np.ones(t850_1.shape))/t850_1)-5.09*np.log(t850_1))
ws_1 = 0.622*(es_1/850.)
thetae_1 = (t850_1+((2450000/1005.7)*((rh850_1/100)*ws_1))) * (1000/850)**(0.2854)

t850_2 = data2.select(name='Temperature',level=850)[0].values
rh850_2 = data2.select(name='Relative humidity',level=850)[0].values
es_2 = 6.11*np.exp((53.49*np.ones(t850_2.shape))-((6808.*np.ones(t850_2.shape))/t850_2)-5.09*np.log(t850_2))
ws_2 = 0.622*(es_2/850.)
thetae_2 = (t850_2+((2450000/1005.7)*((rh850_2/100)*ws_2))) * (1000/850)**(0.2854)
thetae_dif = thetae_2 - thetae_1
t850_dif = t850_2 - t850_1

t850_3 = data3.select(name='Temperature',level=850)[0].values
rh850_3 = data3.select(name='Relative humidity',level=850)[0].values
es_3 = 6.11*np.exp((53.49*np.ones(t850_3.shape))-((6808.*np.ones(t850_3.shape))/t850_3)-5.09*np.log(t850_3))
ws_3 = 0.622*(es_3/850.)
thetae_3 = (t850_3+((2450000/1005.7)*((rh850_3/100)*ws_3))) * (1000/850)**(0.2854)

t850_1 = t850_1 - 273.15
t850_2 = t850_2 - 273.15
t850_3 = t850_3 - 273.15

# 850 mb heights
z850_1 = data1.select(name='Geopotential Height',level=850)[0].values * 0.1
z850_1 = ndimage.filters.gaussian_filter(z850_1, 6.89)
z850_2 = data2.select(name='Geopotential Height',level=850)[0].values * 0.1
z850_2 = ndimage.filters.gaussian_filter(z850_2, 6.89)
z850_dif = z850_2 - z850_1
z850_3 = data3.select(name='Geopotential Height',level=850)[0].values * 0.1
z850_3 = ndimage.filters.gaussian_filter(z850_3, 6.89)

# 850-mb winds
u850_1 = data1.select(name='U component of wind',level=850)[0].values * 1.94384
u850_2 = data2.select(name='U component of wind',level=850)[0].values * 1.94384
u850_3 = data3.select(name='U component of wind',level=850)[0].values * 1.94384
v850_1 = data1.select(name='V component of wind',level=850)[0].values * 1.94384
v850_2 = data2.select(name='V component of wind',level=850)[0].values * 1.94384
v850_3 = data3.select(name='V component of wind',level=850)[0].values * 1.94384
# Rotate winds from grid relative to Earth relative
u850_1, v850_1 = ncepy.rotate_wind(Lat0,Lon0,lon,u850_1,v850_1,'lcc',inverse=False)
u850_2, v850_2 = ncepy.rotate_wind(Lat0,Lon0,lon2,u850_2,v850_2,'lcc',inverse=False)
u850_3, v850_3 = ncepy.rotate_wind(Lat0_3,Lon0_3,lon3,u850_3,v850_3,'lcc',inverse=False)
wspd850_1 = np.sqrt(u850_1**2 + v850_1**2)
wspd850_2 = np.sqrt(u850_2**2 + v850_2**2)
wspd850_3 = np.sqrt(u850_3**2 + v850_3**2)
wspd850_dif = wspd850_2 - wspd850_1

# 700-mb omega and relative humidity
omg700_1 = data1.select(name='Vertical velocity',level=700)[0].values
omg700_2 = data2.select(name='Vertical velocity',level=700)[0].values
omg700_3 = data3.select(name='Vertical velocity',level=700)[0].values
rh700_1 = data1.select(name='Relative humidity',level=700)[0].values
rh700_2 = data2.select(name='Relative humidity',level=700)[0].values
rh700_3 = data3.select(name='Relative humidity',level=700)[0].values
rh700_dif = rh700_2 - rh700_1

# 500 mb height, wind, vorticity
z500_1 = data1.select(name='Geopotential Height',level=500)[0].values * 0.1
z500_1 = ndimage.filters.gaussian_filter(z500_1, 6.89)
z500_2 = data2.select(name='Geopotential Height',level=500)[0].values * 0.1
z500_2 = ndimage.filters.gaussian_filter(z500_2, 6.89)
z500_3 = data3.select(name='Geopotential Height',level=500)[0].values * 0.1
z500_3 = ndimage.filters.gaussian_filter(z500_3, 6.89)
z500_dif = z500_2 - z500_1
vort500_1 = data1.select(name='Absolute vorticity',level=500)[0].values * 100000
vort500_1 = ndimage.filters.gaussian_filter(vort500_1,1.7225)
vort500_2 = data2.select(name='Absolute vorticity',level=500)[0].values * 100000
vort500_2 = ndimage.filters.gaussian_filter(vort500_2,1.7225)
vort500_3 = data3.select(name='Absolute vorticity',level=500)[0].values * 100000
vort500_3 = ndimage.filters.gaussian_filter(vort500_3,1.7225)
u500_1 = data1.select(name='U component of wind',level=500)[0].values * 1.94384
u500_2 = data2.select(name='U component of wind',level=500)[0].values * 1.94384
u500_3 = data3.select(name='U component of wind',level=500)[0].values * 1.94384
v500_1 = data1.select(name='V component of wind',level=500)[0].values * 1.94384
v500_2 = data2.select(name='V component of wind',level=500)[0].values * 1.94384
v500_3 = data3.select(name='V component of wind',level=500)[0].values * 1.94384
# Rotate winds from grid relative to Earth relative
u500_1, v500_1 = ncepy.rotate_wind(Lat0,Lon0,lon,u500_1,v500_1,'lcc',inverse=False)
u500_2, v500_2 = ncepy.rotate_wind(Lat0,Lon0,lon2,u500_2,v500_2,'lcc',inverse=False)
u500_3, v500_3 = ncepy.rotate_wind(Lat0_3,Lon0_3,lon3,u500_3,v500_3,'lcc',inverse=False)

# 250 mb winds
u250_1 = data1.select(name='U component of wind',level=250)[0].values * 1.94384
u250_2 = data2.select(name='U component of wind',level=250)[0].values * 1.94384
u250_3 = data3.select(name='U component of wind',level=250)[0].values * 1.94384
v250_1 = data1.select(name='V component of wind',level=250)[0].values * 1.94384
v250_2 = data2.select(name='V component of wind',level=250)[0].values * 1.94384
v250_3 = data3.select(name='V component of wind',level=250)[0].values * 1.94384
# Rotate winds from grid relative to Earth relative
u250_1, v250_1 = ncepy.rotate_wind(Lat0,Lon0,lon,u250_1,v250_1,'lcc',inverse=False)
u250_2, v250_2 = ncepy.rotate_wind(Lat0,Lon0,lon2,u250_2,v250_2,'lcc',inverse=False)
u250_3, v250_3 = ncepy.rotate_wind(Lat0_3,Lon0_3,lon3,u250_3,v250_3,'lcc',inverse=False)
wspd250_1 = np.sqrt(u250_1**2 + v250_1**2)
wspd250_2 = np.sqrt(u250_2**2 + v250_2**2)
wspd250_3 = np.sqrt(u250_3**2 + v250_3**2)
wspd250_dif = wspd250_2 - wspd250_1

# Precipitable water
pw_1 = data1.select(name='Precipitable water',level=0)[0].values * 0.0393701
pw_2 = data2.select(name='Precipitable water',level=0)[0].values * 0.0393701
pw_dif = pw_2 - pw_1
pw_3 = data3.select(name='Precipitable water',level=0)[0].values * 0.0393701

# Visibility
vis_1 = data1.select(name='Visibility',typeOfLevel='surface')[0].values * 0.000621371
vis_2 = data2.select(name='Visibility',typeOfLevel='surface')[0].values * 0.000621371
vis_dif = vis_2 - vis_1
vis_4 = data4.select(name='Visibility',typeOfLevel='surface')[0].values * 0.000621371

# Cloud Ceiling Height
ceil_1 = data1.select(name='Geopotential Height',nameOfFirstFixedSurface='215')[0].values * (3.28084/1000)
ceil_2 = data2.select(name='Geopotential Height',nameOfFirstFixedSurface='215')[0].values * (3.28084/1000)
ceil_dif = ceil_2 - ceil_1
ceil_4 = data4.select(name='Ceiling')[0].values * (3.28084/1000)

# Composite reflectivity
refc_1 = data1.select(name='Maximum/Composite radar reflectivity')[0].values 
refc_2 = data2.select(name='Maximum/Composite radar reflectivity')[0].values 
refc_5 = data5.select(name='Maximum/Composite radar reflectivity')[0].values 

# 1-km AGL reflectivity
ref1km_1 = data1.select(name='Derived radar reflectivity',level=1000)[0].values
ref1km_2 = data2.select(name='Derived radar reflectivity',level=1000)[0].values
ref1km_5 = data5.select(name='Derived radar reflectivity',level=1000)[0].values 

# PBL height
hpbl_1 = data1.select(name='Planetary boundary layer height')[0].values
hpbl_2 = data2.select(name='Planetary boundary layer height')[0].values
hpbl_dif = hpbl_2 - hpbl_1
hpbl_3 = data3.select(name='Planetary boundary layer height')[0].values

# 0-3 km Storm Relative Helicity
hel3km_1 = data1.select(name='Storm relative helicity',topLevel=3000,bottomLevel=0)[0].values
hel3km_2 = data2.select(name='Storm relative helicity',topLevel=3000,bottomLevel=0)[0].values
hel3km_dif = hel3km_2 - hel3km_1
hel3km_3 = data3.select(name='Storm relative helicity',topLevel=3000,bottomLevel=0)[0].values

# 0-1 km Storm Relative Helicity
hel1km_1 = data1.select(name='Storm relative helicity',topLevel=1000,bottomLevel=0)[0].values
hel1km_2 = data2.select(name='Storm relative helicity',topLevel=1000,bottomLevel=0)[0].values
hel1km_dif = hel1km_2 - hel1km_1
hel1km_3 = data3.select(name='Storm relative helicity',topLevel=1000,bottomLevel=0)[0].values



t2a = time.clock()
t3a = round(t2a-t1a, 3)
print("%.3f seconds to read all messages") % t3a

# colors for difference plots, only need to define once
difcolors = ['blue','#1874CD','dodgerblue','deepskyblue','turquoise','white','white','#EEEE00','#EEC900','darkorange','orangered','red']

# colors for cloud cover plots
cdccolors = ['#FFFFFF','#F0F0F0','#E0E0E0','#D8D8D8','#C8C8C8','#B8B8B8','#A8A8A8','#909090','#787878','#696969']

########################################
#    START PLOTTING FOR EACH DOMAIN    #
########################################

def main():

  # Number of processes must coincide with the number of domains to plot
  pool = multiprocessing.Pool(len(domains))
  pool.map(plot_all,domains)

def plot_all(dom):

  t1dom = time.clock()
  print('Working on '+dom)

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
    x5,y5 = m(lon5,lat5)
  
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

################################
  # Plot SLP
################################
  t1 = time.clock()
  print('Working on slp for '+dom)

  units = 'mb'
  clevs = [976,980,984,988,992,996,1000,1004,1008,1012,1016,1020,1024,1028,1032,1036,1040,1044,1048,1052]
  clevsdif = [-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12]

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs1_a = m.contourf(x,y,slpsmooth1,clevs,cmap=plt.cm.Spectral_r,extend='both',ax=ax)  
      cbar1 = m.colorbar(cs1_a,ax=ax,location='bottom',pad=0.05)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      cs1_b = m.contour(x,y,slpsmooth1,np.arange(940,1060,4),colors='black',linewidths=1.25,ax=ax)
      plt.clabel(cs1_b,np.arange(940,1060,4),inline=1,fmt='%d',fontsize=6,zorder=12,ax=ax)
      ax.text(.5,1.03,prod_str+' SLP ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)
  # plot highs and lows - window parameter controls the number of highs and lows detected
      ncepy.plt_highs_and_lows(m,slp_1,lon,lat,mode='reflect',window='500')

    elif par == 2:
      cs2_a = m.contourf(x2,y2,slpsmooth2,clevs,cmap=plt.cm.Spectral_r,extend='both',ax=ax)
      cbar2 = m.colorbar(cs2_a,ax=ax,location='bottom',pad=0.05)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      cs2_b = m.contour(x2,y2,slpsmooth2,np.arange(940,1060,4),colors='black',linewidths=1.25,ax=ax)
      plt.clabel(cs2_b,np.arange(940,1060,4),inline=1,fmt='%d',fontsize=6,zorder=12,ax=ax)
      ax.text(.5,1.03,para_str+' SLP ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)
  # plot highs and lows - window parameter controls the number of highs and lows detected
      ncepy.plt_highs_and_lows(m,slp_2,lon,lat,mode='reflect',window='500')

    elif par == 3:
      cs = m.contourf(x2,y2,slp_dif,clevsdif,colors=difcolors,extend='both',ax=ax)  
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' SLP ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs3_a = m.contourf(x3,y3,slpsmooth3,clevs,cmap=plt.cm.Spectral_r,extend='both',ax=ax)
      cbar4 = m.colorbar(cs3_a,ax=ax,location='bottom',pad=0.05)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      cs3_b = m.contour(x3,y3,slpsmooth3,np.arange(940,1060,4),colors='black',linewidths=1.25,ax=ax)
      plt.clabel(cs3_b,np.arange(940,1060,4),inline=1,fmt='%d',fontsize=6,zorder=12,ax=ax)
      ax.text(.5,1.03,'RAP SLP ('+units+') \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)
  # plot highs and lows - window parameter controls the number of highs and lows detected
      ncepy.plt_highs_and_lows(m,slp_3,lon3,lat3,mode='reflect',window='75')

    par += 1
  par = 1

  compress_and_save('compareslp_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./compareslp_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot slp for: '+dom) % t3


#################################
  # Plot 2-m T
#################################
  t1 = time.clock()
  print('Working on t2m for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = u'\xb0''F'
  clevs = np.linspace(-16,134,51)
  clevsdif = [-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12]
  cm = cmap_t2m()
  norm = matplotlib.colors.BoundaryNorm(clevs, cm.N)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,tmp2m_1,clevs,cmap=cm,extend='both',ax=ax)
      cs_1.cmap.set_under('white')
      cs_1.cmap.set_over('white')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=[-16,-4,8,20,32,44,56,68,80,92,104,116,128])
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' 2-m Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)
 
    elif par == 2:
      cs_2 = m.contourf(x2,y2,tmp2m_2,clevs,cmap=cm,extend='both',ax=ax)
      cs_2.cmap.set_under('white')
      cs_2.cmap.set_over('white')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=[-16,-4,8,20,32,44,56,68,80,92,104,116,128])
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' 2-m Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))       
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,tmp2m_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 2-m Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2')) 
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_4 = m.contourf(x4,y4,tmp2m_4,clevs,cmap=cm,extend='both',ax=ax)
      cs_4.cmap.set_under('white')
      cs_4.cmap.set_over('white')
      cbar4.remove()
      cbar4 = m.colorbar(cs_4,ax=ax,location='bottom',pad=0.05,ticks=[-16,-4,8,20,32,44,56,68,80,92,104,116,128])
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'URMA 2-m Temperature ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))       
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compare2mt_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./compare2mt_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 2mt for: '+dom) % t3


#################################
# Plot SFC T
#################################
  t1 = time.clock()
  print('Working on tsfc for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = u'\xb0''F'
  clevs = np.linspace(-16,134,51)
  clevsdif = [-18,-15,-12,-9,-6,-3,0,3,6,9,12,15,18]
  cm = cmap_t2m()
  norm = matplotlib.colors.BoundaryNorm(clevs, cm.N)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,tmpsfc_1,clevs,cmap=cm,extend='both',ax=ax)
      cs_1.cmap.set_under('white')
      cs_1.cmap.set_over('white')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=[-16,-4,8,20,32,44,56,68,80,92,104,116,128])
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' 2-m Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)
 
    elif par == 2:
      cs_2 = m.contourf(x2,y2,tmpsfc_2,clevs,cmap=cm,extend='both',ax=ax)
      cs_2.cmap.set_under('white')
      cs_2.cmap.set_over('white')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=[-16,-4,8,20,32,44,56,68,80,92,104,116,128])
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' 2-m Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))       
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,tmpsfc_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' Skin Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,tmpsfc_3,clevs,cmap=cm,extend='both',ax=ax)
      cs_3.cmap.set_under('white')
      cs_3.cmap.set_over('white')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=[-16,-4,8,20,32,44,56,68,80,92,104,116,128])
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'RAP Surface Temperature ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))       
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparetsfc_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./compare2mt_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot tsfc for: '+dom) % t3


#################################
  # Plot 2-m Dew Point
#################################
  t1 = time.clock()
  print('Working on 2mdew for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = u'\xb0''F'
  clevs = np.linspace(-5,80,35)
  clevsdif = [-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12]
  cm = ncepy.cmap_q2m()
  norm = matplotlib.colors.BoundaryNorm(clevs, cm.N)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,dew2m_1,clevs,cmap=cm,extend='both',ax=ax)
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' 2-m Dew Point Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,dew2m_2,clevs,cmap=cm,extend='both',ax=ax)
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' 2-m Dew Point Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,dew2m_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 2-m Dew Point Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_4 = m.contourf(x4,y4,dew2m_4,clevs,cmap=cm,extend='both',ax=ax)
      cbar4.remove()
      cbar4 = m.colorbar(cs_4,ax=ax,location='bottom',pad=0.05)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'URMA 2-m Dew Point Temperature ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compare2mdew_'+dom+'_f'+fhour+'.png')
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 2mdew for: '+dom) % t3

#################################
  # Plot 10-m WSPD
#################################
  t1 = time.clock()
  print('Working on 10mwspd for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'kts'
  if dom == 'conus':
    skip = 80
  elif dom == 'SE':
    skip = 35
  elif dom == 'CO' or dom == 'LA' or dom == 'MA':
    skip = 12
  elif dom == 'BN':
    skip = 10
  elif dom == 'SP':
    skip = 9
  elif dom == 'SF':
    skip = 3
  else:
    skip = 20
  rap_skip = 15
  urma_skip = 75 
  barblength = 4
  clevs = [5,10,15,20,25,30,35,40,45,50,55,60]
  clevsdif = [-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12]
  colorlist = ['turquoise','dodgerblue','blue','#FFF68F','#E3CF57','peru','brown','crimson','red','fuchsia','DarkViolet']

  # Rotate winds to gnomonic projection
  urot_1, vrot_1 = m.rotate_vector(uwind_1,vwind_1,lon,lat)
  urot_2, vrot_2 = m.rotate_vector(uwind_2,vwind_2,lon2,lat2)
  urot_4, vrot_4 = m.rotate_vector(uwind_4,vwind_4,lon4,lat4)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,wspd10m_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_1[::skip,::skip],vrot_1[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,prod_str+' 10-m Winds ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)
    
    elif par == 2:
      cs_2 = m.contourf(x2,y2,wspd10m_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_2[::skip,::skip],vrot_2[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,para_str+' 10-m Winds ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,wspd10m_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label('kts',fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 10-m Wind Speed (kts) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))       
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_4 = m.contourf(x4,y4,wspd10m_4,clevs,colors=colorlist,extend='max',ax=ax)
      cs_4.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_4,ax=ax,location='bottom',pad=0.05)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      m.barbs(lon4[::urma_skip,::urma_skip],lat4[::urma_skip,::urma_skip],urot_4[::urma_skip,::urma_skip],vrot_4[::urma_skip,::urma_skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,'URMA 10-m Winds ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compare10mwind_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./compare10mwind_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 10mwspd for: '+dom) % t3

#################################
  # Plot Most Unstable CAPE/CIN
#################################
  t1 = time.clock()
  print('Working on mucapecin for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'J/kg'
  clevs = [100,250,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]
  clevs2 = [-2000,-500,-250,-100,-25]
  clevsdif = [-2000,-1500,-1000,-500,-250,-100,0,100,250,500,1000,1500,2000]
  colorlist = ['blue','dodgerblue','cyan','mediumspringgreen','#FAFAD2','#EEEE00','#EEC900','darkorange','crimson','darkred','darkviolet']

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,mucape_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=4)
      cs_1b = m.contourf(x,y,mucin_1,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,prod_str+' MUCAPE (shaded) and MUCIN (hatched) ('+units+') \n <-500 (*), -500<-250 (+), -250<-100 (/), -100<-25 (.) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,mucape_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=4)
      cs_2b = m.contourf(x2,y2,mucin_2,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,para_str+' MUCAPE (shaded) and MUCIN (hatched) ('+units+') \n <-500 (*), -500<-250 (+), -250<-100 (/), -100<-25 (.) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,mucape_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05,ticks=clevsdif)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=4)
      ax.text(.5,1.03,para_str+' - '+prod_str+' Most Unstable CAPE ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,mucape_3,clevs,colors=colorlist,extend='max',ax=ax)
      cs_3.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=4)
      cs_3b = m.contourf(x3,y3,mucin_3,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,'RAP MUCAPE (shaded) and MUCIN (hatched) ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparemucape_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparemucape_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot mucapecin for: '+dom) % t3

#################################
  # Plot Surface-Based CAPE/CIN
#################################
  t1 = time.clock()
  print('Working on sfcapecin for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,cape_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=4)
      cs_1b = m.contourf(x,y,sfcin_1,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,prod_str+' SFCAPE (shaded) and SFCIN (hatched) ('+units+') \n <-500 (*), -500<-250 (+), -250<-100 (/), -100<-25 (.) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,cape_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=4)
      cs_2b = m.contourf(x2,y2,sfcin_2,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,para_str+' SFCAPE (shaded) and SFCIN (hatched) ('+units+') \n <-500 (*), -500<-250 (+), -250<-100 (/), -100<-25 (.) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,cape_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05,ticks=clevsdif)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=4)
      ax.text(.5,1.03,para_str+' - '+prod_str+' Surface-Based CAPE ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,cape_3,clevs,colors=colorlist,extend='max',ax=ax)
      cs_3.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=4)
      cs_3b = m.contourf(x3,y3,sfcin_3,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,'RAP SBCAPE (shaded) and SBCIN (hatched) ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparesfcape_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparesfcape_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot sfcapecin for: '+dom) % t3


#################################
  # Plot Mixed Layer CAPE/CIN
#################################
  t1 = time.clock()
  print('Working on mlcapecin for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,mlcape_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=4)
      cs_1b = m.contourf(x,y,mlcin_1,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,prod_str+' MLCAPE (shaded) and MLCIN (hatched) ('+units+') \n <-500 (*), -500<-250 (+), -250<-100 (/), -100<-25 (.) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,mlcape_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=4)
      cs_2b = m.contourf(x2,y2,mlcin_2,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,para_str+' MLCAPE (shaded) and MLCIN (hatched) ('+units+') \n <-500 (*), -500<-250 (+), -250<-100 (/), -100<-25 (.) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,mlcape_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05,ticks=clevsdif)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=4)
      ax.text(.5,1.03,para_str+' - '+prod_str+' Mixed Layer CAPE ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,mlcape_3,clevs,colors=colorlist,extend='max',ax=ax)
      cs_3.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=4)
      cs_3b = m.contourf(x3,y3,mlcin_3,clevs2,colors='none',hatches=['**','++','////','..'],ax=ax)
      ax.text(.5,1.05,'RAP MLCAPE (shaded) and MLCIN (hatched) ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparemlcape_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparemlcape_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot mlcapecin for: '+dom) % t3


#################################
  # Plot 850-mb THETAE/WIND
#################################
  t1 = time.clock()
  print('Working on 850 mb Theta-e for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'K'
# Wind barb density settings for 850, 500, and 250 mb plots
  if dom == 'conus':
    skip = 100
  elif dom == 'SE':
    skip = 40
  elif dom == 'CO' or dom == 'LA' or dom =='MA':
    skip = 18
  elif dom == 'BN':
    skip = 15
  elif dom == 'SP':
    skip = 13
  elif dom == 'SF':
    skip = 4
  else:
    skip = 30
  barblength = 4

  clevs = np.linspace(250,340,31)
  clevsdif = [-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12]
  cm = cmap_t850()
  norm = matplotlib.colors.BoundaryNorm(clevs, cm.N)

  # Rotate winds to gnomonic projection
  urot_1, vrot_1 = m.rotate_vector(u850_1,v850_1,lon,lat)
  urot_2, vrot_2 = m.rotate_vector(u850_2,v850_2,lon2,lat2)
  urot_3, vrot_3 = m.rotate_vector(u850_3,v850_3,lon3,lat3)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,thetae_1,clevs,cmap=cm,extend='both',ax=ax)
      cs_1.cmap.set_under('white')
      cs_1.cmap.set_over('white')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=[250,256,262,268,274,280,286,292,298,304,310,316,322,328,334,340])
      cbar1.set_label(units,fontsize=6)   
      cbar1.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_1[::skip,::skip],vrot_1[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,prod_str+' 850 mb $\Theta$e ('+units+') and Winds (kts) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,thetae_2,clevs,cmap=cm,extend='both',ax=ax)
      cs_2.cmap.set_under('white')
      cs_2.cmap.set_over('white')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=[250,256,262,268,274,280,286,292,298,304,310,316,322,328,334,340])
      cbar2.set_label(units,fontsize=6)   
      cbar2.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_2[::skip,::skip],vrot_2[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,para_str+' 850 mb $\Theta$e ('+units+') and Winds (kts) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)
    
    elif par == 3:
      cs = m.contourf(x2,y2,thetae_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)   
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 850 mb $\Theta$e ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,thetae_3,clevs,cmap=cm,extend='both',ax=ax)
      cs_3.cmap.set_under('white')
      cs_3.cmap.set_over('white')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=[250,256,262,268,274,280,286,292,298,304,310,316,322,328,334,340])
      cbar4.set_label(units,fontsize=6)   
      cbar4.ax.tick_params(labelsize=6)
      m.barbs(lon3[::rap_skip,::rap_skip],lat3[::rap_skip,::rap_skip],urot_3[::rap_skip,::rap_skip],vrot_3[::rap_skip,::rap_skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,'RAP 850 mb $\Theta$e ('+units+') and Winds (kts) \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)
    
    par += 1
  par = 1

  compress_and_save('compare850t_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./compare850t_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 850 mb Theta-e for: '+dom) % t3


#################################
  # Plot 850 mb TEMP
#################################
  t1 = time.clock()
  print('Working on 850 mb TEMP for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = u'\xb0''C'
  clevs = np.linspace(-45,45,31)
  clevsdif = [-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12]
  cm = cmap_t850()
  norm = matplotlib.colors.BoundaryNorm(clevs, cm.N)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,t850_1,clevs,cmap=cm,extend='both',ax=ax)
      cs_1.cmap.set_under('white')
      cs_1.cmap.set_over('white')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=[-42,-36,-30,-24,-18,-12,-6,-0,6,12,18,24,30,36,42])
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_1[::skip,::skip],vrot_1[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,prod_str+' 850-mb Temperaature ('+units+') and Winds (kts) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,t850_2,clevs,cmap=cm,extend='both',ax=ax)
      cs_2.cmap.set_under('white')
      cs_2.cmap.set_over('white')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=[-42,-36,-30,-24,-18,-12,-6,-0,6,12,18,24,30,36,42])
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_2[::skip,::skip],vrot_2[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,para_str+' 850-mb Temperaature ('+units+') and Winds (kts) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,t850_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 850-mb Temperature ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,t850_3,clevs,cmap=cm,extend='both',ax=ax)
      cs_3.cmap.set_under('white')
      cs_3.cmap.set_over('white')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=[-42,-36,-30,-24,-18,-12,-6,-0,6,12,18,24,30,36,42])
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      m.barbs(lon3[::rap_skip,::rap_skip],lat3[::rap_skip,::rap_skip],urot_3[::rap_skip,::rap_skip],vrot_3[::rap_skip,::rap_skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,'RAP 850-mb Temperature ('+units+') and Winds (kts) \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compare850temp_'+dom+'_f'+fhour+'.png')
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 850 temp for: '+dom) % t3


#################################
  # Plot 850 mb WIND/HGT
#################################
  t1 = time.clock()
  print('Working on 850 mb WIND/HGT for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'kts'
  clevs = [10,20,30,40,50,60,70,80,90,100,110]
  clevsdif = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30]
  colorlist = ['turquoise','deepskyblue','dodgerblue','#1874CD','blue','beige','khaki','peru','brown','crimson']

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,wspd850_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('red')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_1[::skip,::skip],vrot_1[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
    # x,y = m(lon,lat)  # need to redefine to avoid index error
      cs1_b = m.contour(x,y,z850_1,np.arange(90,420,3),colors='black',linewidths=0.75,ax=ax)
      plt.clabel(cs1_b,np.arange(90,420,3),inline_spacing=1,fmt='%d',fontsize=6,dorder=12,ax=ax)
      ax.text(.5,1.03,prod_str+' 850-mb Heights (dam) and Winds ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,wspd850_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('red')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_2[::skip,::skip],vrot_2[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      cs2_b = m.contour(x2,y2,z850_2,np.arange(90,420,3),colors='black',linewidths=0.75,ax=ax)
      plt.clabel(cs2_b,np.arange(90,420,3),inline_spacing=1,fmt='%d',fontsize=6,dorder=12,ax=ax)
      ax.text(.5,1.03,para_str+' 850-mb Heights (dam) and Winds ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,wspd850_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+'  850-mb Winds ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,wspd850_3,clevs,colors=colorlist,extend='max',ax=ax)
      cs_3.cmap.set_over('red')
      cbar4.remove()
      cbar4 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      m.barbs(lon3[::rap_skip,::rap_skip],lat3[::rap_skip,::rap_skip],urot_3[::rap_skip,::rap_skip],vrot_3[::rap_skip,::rap_skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      cs3_b = m.contour(x3,y3,z850_3,np.arange(90,420,3),colors='black',linewidths=0.75,ax=ax)
      plt.clabel(cs3_b,np.arange(90,420,3),inline_spacing=1,fmt='%d',fontsize=6,dorder=12,ax=ax)
      ax.text(.5,1.03,'RAP 850-mb Heights (dam) and Winds ('+units+') \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compare850wind_'+dom+'_f'+fhour+'.png')
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 850 wind for: '+dom) % t3


#################################
  # Plot 700-mb OMEGA and RH
#################################
  t1 = time.clock()
  print('Working on 700 mb omega and RH for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = '%'
  clevs = [50,60,70,80,90,100]
  clevsw = [-100,-5]
  clevsdif = [-30,-25,-20,-15,-10,-5,-0,5,10,15,20,25,30]

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs1_a = m.contourf(x,y,rh700_1,clevs,cmap=plt.cm.BuGn,ax=ax)
      cbar1.remove()
      cbar1 = m.colorbar(cs1_a,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6) 
      cbar1.ax.tick_params(labelsize=6)
      cs1_b = m.contourf(x,y,omg700_1,clevsw,colors='blue',ax=ax)
      ax.text(.5,1.03,prod_str+' 700 mb $\omega$ (rising motion in blue) and RH ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs2_a = m.contourf(x2,y2,rh700_2,clevs,cmap=plt.cm.BuGn,ax=ax)
      cbar2.remove()
      cbar2 = m.colorbar(cs2_a,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6) 
      cbar2.ax.tick_params(labelsize=6)
      cs2_b = m.contourf(x2,y2,omg700_2,clevsw,colors='blue',ax=ax)
      ax.text(.5,1.03,para_str+' 700 mb $\omega$ (rising motion in blue) and RH ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,rh700_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      cbar3.remove()
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6) 
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 700 mb RH ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs3_a = m.contourf(x3,y3,rh700_3,clevs,cmap=plt.cm.BuGn,ax=ax)
      cbar4.remove()
      cbar4 = m.colorbar(cs3_a,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6) 
      cbar4.ax.tick_params(labelsize=6)
      cs3_b = m.contourf(x3,y3,omg700_3,clevsw,colors='blue',ax=ax)
      ax.text(.5,1.03,'RAP 700 mb $\omega$ (rising motion in blue) and RH ('+units+') \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compare700_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./compare700_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 700 mb $\omega$ and RH for: '+dom) % t3


#################################
  # Plot 500 mb HGT/WIND/VORT
#################################
  t1 = time.clock()
  print('Working on 500 mb Hgt/Wind/Vort for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'x10${^5}$ s${^{-1}}$'
  units = '$\mathregular{x10^{-5}}$ $\mathregular{s^{-1}}$'
  vortlevs = [16,20,24,28,32,36,40]
  clevsdif = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]
  colorlist = ['yellow','gold','goldenrod','orange','orangered','red']

  # Rotate winds to gnomonic projection
  urot_1, vrot_1 = m.rotate_vector(u500_1,v500_1,lon,lat)
  urot_2, vrot_2 = m.rotate_vector(u500_2,v500_2,lon2,lat2)
  urot_3, vrot_3 = m.rotate_vector(u500_3,v500_3,lon3,lat3)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs1_a = m.contourf(x,y,vort500_1,vortlevs,colors=colorlist,extend='both',ax=ax)
      cs1_a.cmap.set_under('white')
      cs1_a.cmap.set_over('darkred')
      cbar1.remove()
      cbar1 = m.colorbar(cs1_a,ax=ax,location='bottom',pad=0.05,ticks=vortlevs)
      cbar1.set_label(units,fontsize=6)   
      cbar1.ax.tick_params(labelsize=6)

      # plot vorticity maxima as black X's
      local_max = extrema(vort500_1,mode='wrap',window=100)
      xhighs = lon[local_max]
      yhighs = lat[local_max]
      highvals = vort500_1[local_max]
      xyplotted = []
      # don't plot if there is already a X within dmin meters
      yoffset = 0.022*(m.ymax - m.ymin)
      dmin = yoffset
      for x,y,p in zip(xhighs, yhighs, highvals):
        if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin and p > 35:
          dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
          if not dist or min(dist) > dmin:
            ax.text(x,y,'x',fontsize=6,fontweight='bold',\
                    ha='center',va='center',color='black',zorder=10)
            xyplotted.append((x,y))

      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_1[::skip,::skip],vrot_1[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='steelblue',ax=ax)
      x,y = m(lon,lat)	# need to redefine to avoid index error
      cs1_b = m.contour(x,y,z500_1,np.arange(486,600,6),colors='black',linewidths=1,ax=ax)
      plt.clabel(cs1_b,np.arange(486,600,6),inline_spacing=1,fmt='%d',fontsize=6,dorder=12,ax=ax)
      ax.text(.5,1.03,prod_str+' 500 mb Heights (dam), Winds (kts), and $\zeta$ ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs2_a = m.contourf(x2,y2,vort500_2,vortlevs,colors=colorlist,extend='both',ax=ax)
      cs2_a.cmap.set_under('white')
      cs2_a.cmap.set_over('darkred')
      cbar2.remove()
      cbar2 = m.colorbar(cs2_a,ax=ax,location='bottom',pad=0.05,ticks=vortlevs)
      cbar2.set_label(units,fontsize=6)   
      cbar2.ax.tick_params(labelsize=6)

      # plot vorticity maxima as black X's
      local_max = extrema(vort500_2,mode='wrap',window=100)
      xhighs = lon[local_max]
      yhighs = lat[local_max]
      highvals = vort500_2[local_max]
      xyplotted = []
      # don't plot if there is already a X within dmin meters
      yoffset = 0.022*(m.ymax - m.ymin)
      dmin = yoffset
      for x,y,p in zip(xhighs, yhighs, highvals):
        if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin and p > 35:
          dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
          if not dist or min(dist) > dmin:
            ax.text(x,y,'x',fontsize=6,fontweight='bold',\
                    ha='center',va='center',color='black',zorder=10)
            xyplotted.append((x,y))

      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_2[::skip,::skip],vrot_2[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='steelblue',ax=ax)
#      x,y = m(lon,lat)	# need to redefine to avoid index error
      cs2_b = m.contour(x2,y2,z500_2,np.arange(486,600,6),colors='black',linewidths=1,ax=ax)
      plt.clabel(cs2_b,np.arange(486,600,6),inline_spacing=1,fmt='%d',fontsize=6,dorder=12,ax=ax)
      ax.text(.5,1.03,para_str+' 500 mb Heights (dam), Winds (kts), and $\zeta$ ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,z500_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      cbar3.remove()
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6) 
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 500 mb Heights (dam) \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs3_a = m.contourf(x3,y3,vort500_3,vortlevs,colors=colorlist,extend='both',ax=ax)
      cs3_a.cmap.set_under('white')
      cs3_a.cmap.set_over('darkred')
      cbar4.remove()
      cbar4 = m.colorbar(cs3_a,ax=ax,location='bottom',pad=0.05,ticks=vortlevs)
      cbar4.set_label(units,fontsize=6)   
      cbar4.ax.tick_params(labelsize=6)

      # plot vorticity maxima as black X's
      local_max = extrema(vort500_3,mode='wrap',window=100)
      xhighs = lon[local_max]
      yhighs = lat[local_max]
      highvals = vort500_3[local_max]
      xyplotted = []
      # don't plot if there is already a X within dmin meters
      yoffset = 0.022*(m.ymax - m.ymin)
      dmin = yoffset
      for x,y,p in zip(xhighs, yhighs, highvals):
        if x < m.xmax and x > m.xmin and y < m.ymax and y > m.ymin and p > 35:
          dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
          if not dist or min(dist) > dmin:
            ax.text(x3,y3,'x',fontsize=6,fontweight='bold',\
                    ha='center',va='center',color='black',zorder=10)
            xyplotted.append((x3,y3))

      m.barbs(lon3[::rap_skip,::rap_skip],lat3[::rap_skip,::rap_skip],urot_3[::rap_skip,::rap_skip],vrot_3[::rap_skip,::rap_skip],latlon=True,length=barblength,linewidth=0.5,color='steelblue',ax=ax)
      x3,y3 = m(lon3,lat3)	# need to redefine to avoid index error
      cs3_b = m.contour(x3,y3,z500_3,np.arange(486,600,6),colors='black',linewidths=1,ax=ax)
      plt.clabel(cs3_b,np.arange(486,600,6),inline_spacing=1,fmt='%d',fontsize=6,dorder=12,ax=ax)
      ax.text(.5,1.03,'RAP 500 mb Heights (dam), Winds (kts), and $\zeta$ ('+units+') \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compare500_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./compare500_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 500 mb Hgt/Wind/Vort for: '+dom) % t3


#################################
  # Plot 250 mb WIND
#################################
  t1 = time.clock()
  print('Working on 250 mb WIND for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'kts'
  clevs = [50,60,70,80,90,100,110,120,130,140,150]
  clevsdif = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30]
  colorlist = ['turquoise','deepskyblue','dodgerblue','#1874CD','blue','beige','khaki','peru','brown','crimson']

  # Rotate winds to gnomonic projection
  urot_1, vrot_1 = m.rotate_vector(u250_1,v250_1,lon,lat)
  urot_2, vrot_2 = m.rotate_vector(u250_2,v250_2,lon2,lat2)
  urot_3, vrot_3 = m.rotate_vector(u250_3,v250_3,lon3,lat3)

  x,y = m(lon,lat)	# need to redefine to avoid index error
  x3,y3 = m(lon3,lat3)	# need to redefine to avoid index error

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,wspd250_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('red')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_1[::skip,::skip],vrot_1[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,prod_str+' 250 mb Winds ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,wspd250_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('red')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      m.barbs(lon[::skip,::skip],lat[::skip,::skip],urot_2[::skip,::skip],vrot_2[::skip,::skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,para_str+' 250 mb Winds ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,wspd250_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      cbar3.remove()
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6) 
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 250 mb Winds ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,wspd250_3,clevs,colors=colorlist,extend='max',ax=ax)
      cs_3.cmap.set_over('red')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      m.barbs(lon3[::rap_skip,::rap_skip],lat3[::rap_skip,::rap_skip],urot_3[::rap_skip,::rap_skip],vrot_3[::rap_skip,::rap_skip],latlon=True,length=barblength,linewidth=0.5,color='black',ax=ax)
      ax.text(.5,1.03,'RAP 250 mb Winds ('+units+') \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1
   
  compress_and_save('compare250wind_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./compare250wind_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 250 mb WIND for: '+dom) % t3


#################################
  # Plot PW
#################################
  t1 = time.clock()
  print('Working on PW for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'in'
  clevs = [0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25]
  clevsdif = [-1.25,-1,-.75,-.5,-.25,-.1,0.,.1,.25,.50,.75,1,1.25]
  colorlist = ['lightsalmon','khaki','palegreen','cyan','turquoise','cornflowerblue','mediumslateblue','darkorchid','deeppink']

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,pw_1,clevs,colors=colorlist,extend='both',ax=ax)
      cs_1.cmap.set_under('white')
      cs_1.cmap.set_over('hotpink')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' Precipitable Water ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,pw_2,clevs,colors=colorlist,extend='both',ax=ax)
      cs_2.cmap.set_under('white')
      cs_2.cmap.set_over('hotpink')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' Precipitable Water ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,pw_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      cbar3.remove()
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05,ticks=clevsdif)
      cbar3.set_label(units,fontsize=6) 
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' Precipitable Water ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,pw_3,clevs,colors=colorlist,extend='both',ax=ax)
      cs_3.cmap.set_under('white')
      cs_3.cmap.set_over('hotpink')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'RAP Precipitable Water ('+units+') \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparepw_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparepw_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot PW for: '+dom) % t3


#################################
  # Plot Surface Visibility
#################################
  t1 = time.clock()
  print('Working on Surface Visibility for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'miles'
  clevs = [0.25,0.5,1,2,3,4,5,10]
  clevsdif = [-15,-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12.5,15]
  colorlist = ['salmon','goldenrod','#EEEE00','palegreen','darkturquoise','blue','mediumpurple']

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,vis_1,clevs,colors=colorlist,extend='min',ax=ax)
      cs_1.cmap.set_under('firebrick')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' Surface Visibility ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,vis_2,clevs,colors=colorlist,extend='min',ax=ax)
      cs_2.cmap.set_under('firebrick')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' Surface Visibility ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,vis_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' Surface Visibility ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_4 = m.contourf(x4,y4,vis_4,clevs[0:7],colors=colorlist[0:6],extend='min',ax=ax)
      cs_4.cmap.set_under('firebrick')
      cbar4.remove()
      cbar4 = m.colorbar(cs_4,ax=ax,location='bottom',pad=0.05,ticks=clevs[0:7])
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'URMA Surface Visibility ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparevis_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparevis_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot Surface Visibility for: '+dom) % t3


#################################
  # Plot Cloud Ceiling Height
#################################
  t1 = time.clock()
  print('Working on Cloud Ceiling Height for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'kft'
  clevs = [0,0.1,0.3,0.5,1,5,10,15,20,25,30,35,40]
  clevsdif = [-12,-10,-8,-6,-4,-2,0.,2,4,6,8,10,12]
  colorlist = ['firebrick','tomato','salmon','lightsalmon','goldenrod','khaki','gold','yellow','palegreen','mediumspringgreen','lime','limegreen']

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,ceil_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('white')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.set_xticklabels(clevs)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' Ceiling ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,ceil_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('white')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.set_xticklabels(clevs)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' Ceiling ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,ceil_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      try:
          cbar3.remove()
      except:
          pass
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05,ticks=clevsdif)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' - '+para_str+' Ceiling ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_4 = m.contourf(x4,y4,ceil_4,clevs,colors=colorlist,ax=ax)
    # cs_4 = m.contourf(x4,y4,ceil_4,clevs,colors=colorlist,extend='max',ax=ax)
    # cs_4.cmap.set_over('white')
      cbar4.remove()
      cbar4 = m.colorbar(cs_4,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'URMA Ceiling ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compareceil_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparezceil_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock() 
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot Ceiling for: '+dom) % t3


#################################
  # Plot PBL height
#################################
  t1 = time.clock()
  print('Working on PBL height for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'm'
  clevs = [50,100,250,500,1000,1500,2000,2500,3000,3500,4000,4500,5000]
  clevsdif = [-1800,-1500,-1200,-900,-600,-300,0,300,600,900,1200,1500,1800]
  colorlist= ['gray','blue','dodgerblue','cyan','mediumspringgreen','#FAFAD2','#EEEE00','#EEC900','darkorange','crimson','darkred','darkviolet']

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,hpbl_1,clevs,colors=colorlist,extend='both',ax=ax)
      cs_1.cmap.set_under('white')
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,ticks=clevs,location='bottom',pad=0.05)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=5)
      ax.text(.5,1.03,prod_str+' PBL Height ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,hpbl_2,clevs,colors=colorlist,extend='both',ax=ax)
      cs_2.cmap.set_under('white')
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,ticks=clevs,location='bottom',pad=0.05)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=5)
      ax.text(.5,1.03,para_str+' PBL Height ('+units+') \n initialized: '+itime+' valid: '+vtime+' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,hpbl_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      cbar3.remove()
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' PBL Height ('+units+') \n initialized: '+itime+' valid: '+vtime+' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,hpbl_3,clevs,colors=colorlist,extend='both',ax=ax)
      cs_3.cmap.set_under('white')
      cs_3.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,ticks=clevs,location='bottom',pad=0.05)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=5)
      ax.text(.5,1.03,'RAP PBL Height ('+units+') \n valid: '+vtime + ' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparehpbl_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparehpbl_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot PBL height for: '+dom) % t3


#################################
  # Plot 0-3 km Storm Relative Helicity
#################################
  t1 = time.clock()
  print('Working on 0-3 km SRH for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'm${^2}$ s$^{-2}$'
  units = '$\mathregular{m^{2}}$ $\mathregular{s^{-2}}$'
  clevs = [50,100,150,200,250,300,400,500,600,700,800]
  clevsdif = [-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120]
  colorlist = ['mediumblue','dodgerblue','chartreuse','limegreen','darkgreen','#EEEE00','orange','orangered','firebrick','darkmagenta']

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,hel3km_1,clevs,colors=colorlist,extend='both',ax=ax)
      cs_1.cmap.set_under('white')
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' 0-3 km Storm Relative Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime+' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,hel3km_2,clevs,colors=colorlist,extend='both',ax=ax)
      cs_2.cmap.set_under('white')
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' 0-3 km Storm Relative Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime+' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,hel3km_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      cbar3.remove()
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 0-3 km Storm Relative Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime+' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,hel3km_3,clevs,colors=colorlist,extend='both',ax=ax)
      cs_3.cmap.set_under('white')
      cs_3.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'RAP 0-3 km Storm Relative Helicity ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparehel3km_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparehel3km_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 0-3 km SRH for: '+dom) % t3


#################################
  # Plot 0-1 km Storm Relative Helicity
#################################
  t1 = time.clock()
  print('Working on 0-1 km SRH for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,hel1km_1,clevs,colors=colorlist,extend='both',ax=ax)
      cs_1.cmap.set_under('white')
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' 0-1 km Storm Relative Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime+' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,hel1km_2,clevs,colors=colorlist,extend='both',ax=ax)
      cs_2.cmap.set_under('white')
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' 0-1 km Storm Relative Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime+' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      cs = m.contourf(x2,y2,hel1km_dif,clevsdif,colors=difcolors,extend='both',ax=ax)
      cs.cmap.set_under('darkblue')
      cs.cmap.set_over('darkred')
      cbar3.remove()
      cbar3 = m.colorbar(cs,ax=ax,location='bottom',pad=0.05)
      cbar3.set_label(units,fontsize=6)
      cbar3.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' - '+prod_str+' 0-1 km Storm Relative Helicity ('+units+') \n initialized: '+itime+' valid: '+vtime+' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_3 = m.contourf(x3,y3,hel1km_3,clevs,colors=colorlist,extend='both',ax=ax)
      cs_3.cmap.set_under('white')
      cs_3.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_3,ax=ax,location='bottom',pad=0.05)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'RAP 0-1 km Storm Relative Helicity ('+units+') \n valid: '+vtime+' (f00)',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparehel1km_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparehel1km_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 0-1 km SRH for: '+dom) % t3


#################################
  # Plot composite reflectivity
#################################
  t1 = time.clock()
  print('Working on composite reflectivity for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'dBZ'
  clevs = np.linspace(5,70,14)
  clevsdif = [20,1000]
  colorlist = ['turquoise','dodgerblue','mediumblue','lime','limegreen','green','#EEEE00','#EEC900','darkorange','red','firebrick','darkred','fuchsia']
  
  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,refc_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' Composite Reflectivity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,refc_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' Composite Reflectivity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      csdif = m.contourf(x,y,refc_1,clevsdif,colors='red')
      csdif2 = m.contourf(x,y,refc_2,clevsdif,colors='dodgerblue')
      cbar3.remove()
      ax.text(.5,1.03,prod_str+' (red) and '+para_str+' (blue) Composite Reflectivity > 20 ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_5 = m.contourf(x5,y5,refc_5,clevs,colors=colorlist,extend='max',ax=ax)
      cs_5.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_5,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'Composite Reflectivity Observations ('+units+') \n valid: '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('comparerefc_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparerefc_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot composite reflectivity for: '+dom) % t3


#################################
  # Plot 1-km AGL reflectivity
#################################
  t1 = time.clock()
  print('Working on 1-km AGL reflectivity for '+dom)

  # Clear off old plottables but keep all the map info
  clear_plotables(ax1,keep_ax_lst_1,fig)
  clear_plotables(ax2,keep_ax_lst_2,fig)
  clear_plotables(ax3,keep_ax_lst_3,fig)
  clear_plotables(ax4,keep_ax_lst_4,fig)

  units = 'dBZ'
  clevs = np.linspace(5,70,14)
  clevsdif = [20,1000]
  colorlist = ['turquoise','dodgerblue','mediumblue','lime','limegreen','green','#EEEE00','#EEC900','darkorange','red','firebrick','darkred','fuchsia']
  
  for ax in axes:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmax = int(round(xmax))
    ymax = int(round(ymax))

    if par == 1:
      cs_1 = m.contourf(x,y,ref1km_1,clevs,colors=colorlist,extend='max',ax=ax)
      cs_1.cmap.set_over('black')
      cbar1.remove()
      cbar1 = m.colorbar(cs_1,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar1.set_label(units,fontsize=6)
      cbar1.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,prod_str+' 1-km AGL Reflectivity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 2:
      cs_2 = m.contourf(x2,y2,ref1km_2,clevs,colors=colorlist,extend='max',ax=ax)
      cs_2.cmap.set_over('black')
      cbar2.remove()
      cbar2 = m.colorbar(cs_2,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar2.set_label(units,fontsize=6)
      cbar2.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,para_str+' 1-km AGL Reflectivity ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 3:
      csdif = m.contourf(x,y,ref1km_1,clevsdif,colors='red')
      csdif2 = m.contourf(x,y,ref1km_2,clevsdif,colors='dodgerblue')
      ax.text(.5,1.03,prod_str+' (red) and '+para_str+' (blue) 1-km AGL Reflectivity > 20 ('+units+') \n initialized: '+itime+' valid: '+vtime + ' (f'+fhour+')',horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    elif par == 4:
      cs_5 = m.contourf(x5,y5,ref1km_5,clevs,colors=colorlist,extend='max',ax=ax)
      cs_5.cmap.set_over('black')
      cbar4.remove()
      cbar4 = m.colorbar(cs_5,ax=ax,location='bottom',pad=0.05,ticks=clevs)
      cbar4.set_label(units,fontsize=6)
      cbar4.ax.tick_params(labelsize=6)
      ax.text(.5,1.03,'1-km AGL Reflectivity Observations ('+units+') \n valid: '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
      ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

    par += 1
  par = 1

  compress_and_save('compareref1km_'+dom+'_f'+fhour+'.png')
#  plt.savefig('./comparerefc_'+dom+'_f'+fhour+'.png', bbox_inches='tight',dpi=150)
  t2 = time.clock()
  t3 = round(t2-t1, 3)
  print('%.3f seconds to plot 1-km AGL reflectivity for: '+dom) % t3




######################################################

  t3dom = round(t2-t1dom, 3)
  print("%.3f seconds to plot all variables for: "+dom) % t3dom
  plt.clf()

######################################################

#main()
plot_all('conus')
