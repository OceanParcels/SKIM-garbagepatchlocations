# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 16:58:23 2018

@author: Victor Onink
Loading in the downloaded Stokes data so that we can compute daily means and save those into seperate files
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
    #time.conventions='relative julian days with decimal part (as parts of the day)'
    latitude = root_grp.createVariable('latitude', 'f8', ('latitude',))
    longitude = root_grp.createVariable('longitude', 'f8', ('longitude',))
    uuss = root_grp.createVariable('uuss', 'f8', ('time', 'latitude', 'longitude',))
    vuss = root_grp.createVariable('vuss','f8',('time','latitude','longitude'))
    #Now assigning values to the variables
    latitude[:]=lat
    longitude[:]=lon
    time[0]=timeStamp
    uuss[0,:,:]=avgU
    vuss[0,:,:]=avgV
    #Finally, close the new field
    root_grp.close()

#Now we need to recompute all the Stokes files so that we can get daily means and
#save those into netCDF files for Parcels to read in.

for i in range(0,37983,8):
    uussMean=np.zeros((8,317,720))
    vussMean=np.zeros((8,317,720))
    for j in range(8):
        digits=len(str(i))
        number=str(i)
        while len(number)<5:
            number='0'+number
        File='Stokes0'+number+'.nc'
        time=Dataset(File).variables['time'][:]
        lon=Dataset(File).variables['longitude'][:]
        lat=Dataset(File).variables['latitude'][:]
        uussMean[j,:,:]=Dataset(File).variables['uuss'][0,:,:]
        vussMean[j,:,:]=Dataset(File).variables['vuss'][0,:,:]
    uussM=np.nanmean(uussMean,axis=0)
    vussM=np.nanmean(vussMean,axis=0)
    uuss,vuss=np.zeros((1,317,720)),np.zeros((1,317,720))
    uuss[0,:,:]=uussM
    uuss[uuss<-50]=np.nan
    vuss[0,:,:]=vussM
    vuss[vuss<-50]=np.nan
    year=datetime.datetime(1990, 1, 1)
    Date=year + datetime.timedelta(days=time[0])
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
    filename='/scratch/Victor/StokeData24h/StokesDriftDailyMean'+datestamp+'.nc'
    label='Daily Mean Stokes Drift '+datestamp
    netCDFcreator(filename,uuss,vuss,lon,lat,time,label)
        
