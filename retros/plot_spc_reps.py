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
import dawsonpy
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
itime = cycle
if hour%12 == 0:
    fhrs = np.arange(0,25,1)
else:
    fhrs = np.arange(0,13,1)
vtime_list = [ncepy.ndate(itime,x) for x in fhrs]


# Define prod and para strings
model_str = str(sys.argv[3])


###################################################
# Read in all variables and calculate differences #
###################################################
t1a = time.clock()

torn_lats = []
torn_lons = []
torn_times = []
hail_lats = []
hail_lons = []
hail_times = []
wind_lats = []
wind_lons = []
wind_times = []


# Storm reports
if str.lower(case[0:3]) == 'jul' or str.lower(case[0:3]) == 'aug':
   retro_period = 'jul2018'
elif str.lower(case[0:3]) == 'feb' or str.lower(case[0:3]) == 'mar':
   retro_period = 'febmar2019'
elif str.lower(case[0:3]) == 'may':
   retro_period = 'may2019'

vtime = vtime_list[-1]
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



t2a = time.clock()
t3a = round(t2a-t1a, 3)
print("%.3f seconds to read all messages") % t3a


# Specify plotting domains
domains = dawsonpy.get_domains(case, model_str)
#domains = ['conus','nplains','cplains','splains','midwest','northeast','southeast']

plots = [n for n in itertools.product(domains,fhrs[-1:])]

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
  im = image.imread('/gpfs/dell2/emc/modeling/noscrub/Benjamin.Blake/python.fv3/noaa.png')
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
  elif dom == 'cplains':
    llcrnrlon = -105.0
    llcrnrlat = 32.5 
    urcrnrlon = -88.0
    urcrnrlat = 43.5
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
  elif dom == 'northwest':
    llcrnrlon = -125.0
    llcrnrlat = 37.0 
    urcrnrlon = -102.0
    urcrnrlat = 49.5
    lat_0 = 45.0
    lon_0 = -113.5
    xscale=0.15
    yscale=0.18
  elif dom == 'southwest':
    llcrnrlon = -123.5
    llcrnrlat = 30.0 
    urcrnrlon = -100.5
    urcrnrlat = 42.5
    lat_0 = 37.0
    lon_0 = -112.0
    xscale=0.15
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
  # Plot Run-Total 2-5 km UH
#################################
  if (fhr >= 0):
    t1 = time.clock()
    print('Working on Run-Total 2-5 km UH for '+dom)

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

      elif par == 2:

        torn_reps = m.scatter(xt, yt, marker='v', color='r', s=markersize, ax=ax) 

        ax.text(.5,1.03,'Tornado Reports valid: '+itime+' to '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 3:

        hail_reps = m.scatter(xh, yh, marker='o', color='g', s=markersize, ax=ax) 

        ax.text(.5,1.03,'Hail Reports valid: '+itime+' to '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)

      elif par == 4:

        wind_reps = m.scatter(xw, yw, marker='s', color='b', s=markersize, ax=ax) 

        ax.text(.5,1.03,'Wind Reports valid: '+itime+' to '+vtime,horizontalalignment='center',fontsize=6,transform=ax.transAxes,bbox=dict(facecolor='white',alpha=0.85,boxstyle='square,pad=0.2'))
        ax.imshow(im,aspect='equal',alpha=0.5,origin='upper',extent=(0,int(round(xmax*xscale)),0,int(round(ymax*yscale))),zorder=4)


      par += 1
    par = 1

    compress_and_save('reports_'+itime+'_'+dom+'.png')
    t2 = time.clock()
    t3 = round(t2-t1, 3)
    print('%.3f seconds to plot Run-Total 2-5 km UH for: '+dom) % t3



######################################################

  t3dom = round(t2-t1dom, 3)
  print("%.3f seconds to plot all variables for: "+dom) % t3dom
  plt.clf()

######################################################

main()
#plot_all('conus')
#plot_all(6)




