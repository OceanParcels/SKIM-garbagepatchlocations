# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 15:34:14 2018

@author: Victor Onink
"""

from netCDF4 import Dataset
import numpy as np
def AreaCalc(sizeLat,sizeLon): #Calculate surface area of grid cells
    deg2rd = np.pi/180.
    r=6378.1
    lon_bins = np.linspace(0,360., sizeLon+1)
    lat_bins = np.linspace(-80, 80, sizeLat+1) 
    Area=np.array([[deg2rd*(lon_bins[i+1]-lon_bins[i])*(np.sin(deg2rd*lat_bins[j+1])
                            - np.sin(deg2rd*lat_bins[j])) for i in range(len(lon_bins)-1)] 
                            for j in range(len(lat_bins)-1)])
    Area=r*r*Area*1000*1000 #convert it to m^2 instead of km^2
    return Area
def HistogramFunction(londata,latdata,correction):
    londata,latdata=londata.reshape(np.size(londata)),latdata.reshape(np.size(latdata))
    binsLon=np.arange(-180,180)
    binsLat=np.arange(-80,80)
    density=np.zeros((len(binsLon),len(binsLat)))
    for i in range(np.array(londata).shape[0]):
        density[np.argmin(np.abs(londata[i]-binsLon)),np.argmin(np.abs(latdata[i]-binsLat))]+=1#*correction
    #Now, normalize it by area
    area=AreaCalc(len(binsLat),len(binsLon)).T
    density/=area
    density[density==0]=np.nan
    return density
#%%
File=[
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal24hBeach.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal24h.nc']
open_ocean=np.array([0,0])
for i in range(len(File)):
    print File[i]
    dataset=Dataset(File[i])
    lat=dataset.variables['lat'][:]
    lon=dataset.variables['lon'][:]
    time=dataset.variables['time'][:]
    lon[lon>180]-=360
    #Remove all the beach particles, which requires them to be stuck for 10 steps, so 20 days
    for k in range(lon.shape[0]):
        if lon[k,-1]==lon[k,-10]:
            if lat[k,-1]==lat[k,-10]:
                lat[k,:]=np.nan
                lon[k,:]=np.nan
    open_ocean[i]=np.count_nonzero(~np.isnan(lon[:,0]))
correction=[float(open_ocean[1])/open_ocean[0],1]
#%%
File=[
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal24hBeach.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal24h.nc']
saveFiles=['AtlanticBeach','AtlanticNoBeach']
location='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles/Onink et al/Densities/'
for i in range(len(File)):
    print File[i]
    dataset=Dataset(File[i])
    lat=dataset.variables['lat'][:]
    lon=dataset.variables['lon'][:]
    time=dataset.variables['time'][:]
    lon[lon>180]-=360
    #Remove all the beach particles, which requires them to be stuck for 10 steps, so 20 days
    for k in range(lon.shape[0]):
        if lon[k,-1]==lon[k,-10]:
            if lat[k,-1]==lat[k,-10]:
                lat[k,:]=np.nan
                lon[k,:]=np.nan
#    open_ocean=np.count_nonzero(~np.isnan(lon[:,0]))
    #Now, we want the last 5 years of particle positions, since in this time
    #the garbage patch has largely been formed
    if lon.shape[1]==4748:
        lonLast=lon[:,-365*2:]
        latLast=lat[:,-365*2:]
        Time=time[:,-365*2:]
    else:
        lonLast=lon[:,-185:]
        latLast=lat[:,-185:]
        Time=time[:,-185:]
    density=np.zeros((lonLast.shape[1],360,160))
    for j in range(lonLast.shape[1]):
        density[j,:,:]=HistogramFunction(lonLast[:,j],latLast[:,j],correction[i])
    density.dump(location+saveFiles[i])
#    Time.dump(location+saveFiles[i]+'Time')