# Author: L Dawson
#
# Script to pull fcst files from HPSS, rename, and save in desired data directory

import numpy as np
import datetime, time, os, sys, subprocess
import urllib2
import re, csv, glob
import dawsonpy


#Determine initial date/time
try:
   retro_period = str(sys.argv[1])
except IndexError:
   retro_period = None

if retro_period is None:
   retro_period = raw_input('Enter retro period (jul2018, febmar2019, or may2019): ')

if str.lower(retro_period) == 'jul2018':
    first_date = datetime.datetime(2018,07,14,00,00)
    last_date  = datetime.datetime(2018,8,16,00,00)
elif str.lower(retro_period) == 'febmar2019':
    first_date = datetime.datetime(2019,01,31,00,00)
    last_date  = datetime.datetime(2019,03,11,00,00)
elif str.lower(retro_period) == 'may2019':
    first_date = datetime.datetime(2019,04,30,00,00)
    last_date  = datetime.datetime(2019,06,01,00,00)


# Build list of dates for downloading files
date_list = []
new_date = first_date
while new_date <= last_date:
    date_list.append(new_date.strftime('%Y%m%d'))
    new_date = new_date + datetime.timedelta(hours=24)


# Decide where to write report file
machine, hostname = dawsonpy.get_machine()

if machine == 'HERA':
   pass
elif machine == 'WCOSS':
    REP_DIR = '/gpfs/'+hostname[0]+'d1/emc/meso/noscrub/Logan.Dawson/RAPv5_HRRRv4_retros/reports'
elif machine == 'WCOSS_C':
   pass
elif machine == 'WCOSS_DELL_P3':
    REP_DIR = '/gpfs/dell2/emc/verification/noscrub/Logan.Dawson/RAPv5_HRRRv4_retros/reports'

if not os.path.exists(REP_DIR):
    os.makedirs(REP_DIR)


# Begin processing SPC reports
hazards = ['torn','hail','wind']

for hazard in hazards:
    times = []
    lats = []
    lons = []
    mags = []

    for date in date_list:

        print("Reading "+date+" "+hazard+" reports file")
        url = 'https://www.spc.noaa.gov/climo/reports/'+date[2:]+'_rpts_filtered_'+hazard+'.csv'
     
        f = urllib2.urlopen(url)
        reader = csv.reader(f)
      # headers = next(reader)

        for row in reader:
            if row[0] != 'Time':   # skip header row
                try:
                    float(row[5])
                    good_report = True
                except:
                    good_report = False
                    print("Bad report in "+date+" "+hazard+" reports file. Skipping this report.")

                if good_report:
                    if int(row[0][0:2]) < 12:    # correct the date for 00Z to 12Z reports
                        real_date = datetime.datetime.strptime(date,'%Y%m%d')
                        real_date = real_date + datetime.timedelta(days=1)
                        times.append(real_date.strftime('%Y%m%d')+row[0])
                    else:
                        times.append(date+row[0])
                    lats.append(float(row[5]))
                    lons.append(float(row[6]))
                    if hazard == 'hail':
                        mags.append(float(row[1])/100.)
                    elif hazard == 'wind':
                        try:
                            mags.append(int(row[1]))
                        except:
                            mags.append(np.nan)
                    elif hazard == 'torn':
                        mags.append(np.nan)
 

    print("Writing "+str.upper(retro_period)+" "+hazard+" reports file")
    f = open(REP_DIR+'/'+str.lower(retro_period)+'_'+hazard+'.csv','wt')

    try:
        writer = csv.writer(f)
        i = 0
        for i in xrange(len(times)):
            writer.writerow([times[i][0:4],times[i][4:6],times[i][6:8],times[i][8:],lats[i],lons[i],mags[i]])
            i += 1
    finally:
        f.close()


print("Done")
