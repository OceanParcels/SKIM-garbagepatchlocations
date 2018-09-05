# -*- coding: utf-8 -*-
"""
Created on Mon Sep 03 15:16:10 2018

@author: Victor Onink
"""
from netCDF4 import Dataset
import os
import datetime
import numpy as np

def netCDFcreator(filename,avgU,avgV,lon,lat):
    root_grp = Dataset(filename, 'w', format='NETCDF4')
    #Create Dimensions
    latSize=lat.size
    lonSize=lon.size
    root_grp.createDimension('time', 1)
    root_grp.createDimension('latitude', latSize)
    root_grp.createDimension('longitude', lonSize)
    # variables
    time = root_grp.createVariable('time', 'f8', ('time',))
    latitude = root_grp.createVariable('latitude', 'f8', ('latitude',))
    longitude = root_grp.createVariable('longitude', 'f8', ('longitude',))
    Uavg = root_grp.createVariable('Uavg', 'f8', ('time', 'latitude', 'longitude',))
    Vavg = root_grp.createVariable('Vavg','f8',('time','latitude','longitude'))
    #Now assigning values to the variables
    latitude[:]=lat
    longitude[:]=lon
    time[0]=0
    Uavg[0,:,:]=avgU
    Vavg[0,:,:]=avgV
    #Finally, close the new field
    root_grp.close()
#%%

s='/' #this marks the directories, not \ since that fails with numbers in the directory names
standard='-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v03.0-fv01.0.nc'
meanU=np.zeros((720,1440))
meanV=np.zeros((720,1440))
count=0
for i in range(2002,2015):
        direcyear=str(i)
        if i==2004 or i==2008 or i==2012:
            days=366
        else:
            days=365 
        for j in range(1,days+1):
            if j<10:
                day='00'+str(j)
            elif j>=10 and j<100:
                day='0'+str(j)
            else:
                day=str(j)
            Date=datetime.datetime(i, 1, 1) + datetime.timedelta(j - 1)
            [y, m, d]=[Date.year,Date.month,Date.day]
            if m<10:
                if d<10:
                    datestamp=str(y)+'0'+str(m)+'0'+str(d)
                else:
                    datestamp=str(y)+'0'+str(m)+str(d)
            else:
                if d<10:
                    datestamp=str(y)+str(m)+'0'+str(d)
                else:
                    datestamp=str(y)+str(m)+str(d)
            count+=1
            File=datestamp+standard
            dataset=Dataset(File)
            U=dataset.variables['eastward_eulerian_current_velocity'][0,:,:]
            V=dataset.variables['northward_eulerian_current_velocity'][0,:,:]
            U[U==np.nan]=0
            V[V==np.nan]=0
            meanU+=U
            meanV+=V
meanU=meanU/count
meanV=meanV/count             
lat=dataset.variables['lat'][:]
lon=dataset.variables['lon'][:]            
filename='/scratch/Victor/AverageTotalCurrents.nc'
netCDFcreator(filename,meanU,meanV,lon,lat)    