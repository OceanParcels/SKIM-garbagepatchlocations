# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 16:58:23 2018

@author: Victor Onink
computing daily mean wind data from wavewatch
"""

from netCDF4 import Dataset
import os
import datetime
import numpy as np

def netCDFcreator(filename,avgU,avgV,lon,lat,timeStamp,label):
    root_grp = Dataset(filename, 'w', format='NETCDF4')
    root_grp.description = label
    #Create Dimensions
    latSize=lat.size
    lonSize=lon.size
    root_grp.createDimension('time', 1)
    root_grp.createDimension('latitude', latSize)
    root_grp.createDimension('longitude', lonSize)
    # variables
    time = root_grp.createVariable('time', 'f8', ('time',))
    time.units='days since 1990-01-01 00:00:00'
    time.calender='gregorian'
    #time.conventions-'relative julian days with decimal part (as parts of the day)'
    latitude = root_grp.createVariable('latitude', 'f8', ('latitude',))
    longitude = root_grp.createVariable('longitude', 'f8', ('longitude',))
    uwnd = root_grp.createVariable('uwnd', 'f8', ('time', 'latitude', 'longitude',))
    vwnd = root_grp.createVariable('vwnd','f8',('time','latitude','longitude'))
    #Now assigning values to the variables
    latitude[:]=lat
    longitude[:]=lon
    time[0]=timeStamp
    uwnd[0,:,:]=avgU
    vwnd[0,:,:]=avgV
    #Finally, close the new field
    root_grp.close()

#Now we need to recompute all the Stokes files so that we can get daily means and
#save those into netCDF files for Parcels to read in.

for i in range(0,37983,8):
    uussM=np.zeros((8,317,720))
    vussM=np.zeros((8,317,720))
    for j in range(8):
        digits=len(str(i))
        number=str(i)
        while len(number)<5:
            number='0'+number
        File='WaveWatchWind0'+number+'.nc'
        time=Dataset(File).variables['time'][:]
        lon=Dataset(File).variables['longitude'][:]
        lat=Dataset(File).variables['latitude'][:]
        uussM[j,:,:]=Dataset(File).variables['uwnd'][0,:,:]
        vussM[j,:,:]=Dataset(File).variables['vwnd'][0,:,:]
    uussMean=np.nanmean(uussM,axis=0)
    vussMean=np.nanmean(vussM,axis=0)
    uuss,vuss=np.zeros((1,317,720)),np.zeros((1,317,720))
    uuss[0,:,:]=uussMean
    uuss[uuss<-50]=np.nan
    vuss[0,:,:]=vussMean
    vuss[vuss<-50]=np.nan
    year=datetime.datetime(1990, 1, 1)
    Date=year+ datetime.timedelta(days=time[0])
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
    filename='/scratch/Victor/WaveWatchWind24h/WaveWatchWindDailyMean'+datestamp+'.nc'
    label='Daily Mean WaveWatch III Wind Fields '+datestamp
    netCDFcreator(filename,uuss,vuss,lon,lat,time,label)
        
