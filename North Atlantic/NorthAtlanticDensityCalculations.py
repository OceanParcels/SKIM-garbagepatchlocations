# -*- coding: utf-8 -*-
"""
Created on Tue May 29 16:08:35 2018

@author: Victor Onink
We have lon-lat positions for all the particles, but for making figures and such
makes more sense to consider densities. Also, I want to be able to consider densities
for more than just the last time point, and instead of generating those each time
it makes more sense to just generate and save them once, so i can then call upon them
when I want to make figures or animations
Another new example of my coding getting better so that I can improve my already existing
codes which are a whole lot less efficient
There will be a similar one for the pacific datasets and global datasets of course
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
def HistogramFunction(londata,latdata):
    londata,latdata=londata.reshape(np.size(londata)),latdata.reshape(np.size(latdata))
    binsLon=np.arange(-180,180)
    binsLat=np.arange(-80,80)
    density=np.zeros((len(binsLon),len(binsLat)))
    for i in range(np.array(londata).shape[0]):
        density[np.argmin(np.abs(londata[i]-binsLon)),np.argmin(np.abs(latdata[i]-binsLat))]+=1
    #Now, normalize it by area
    area=AreaCalc(len(binsLat),len(binsLon)).T
    density/=area
    density[density==0]=np.nan
    return density
#%%
File=[
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal3h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal24h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticMeanTotal.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\NorthAtlanticEkman.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticGeostrophic.nc',
#      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\NorthAtlanticStoke.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticStokeTotal3h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticStokeTotal.nc',
#      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\NorthAtlanticStoke.nc'
      ]
  
saveFiles=[
           'NorthAtlanticTotalDensity3h','NorthAtlanticTotalDensity24h','NorthAtlanticTotalMeanDensity',
           'NorthAtlanticEkmanDensity','NorthAtlanticGeostrophicDensity',
#           'NorthAtlanticStokesDensity',
           'NorthAtlanticStokesTotalDensity3h','NorthAtlanticStokesTotal24h',
           ]
location='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\Densities/'
#location='/scratch/Victor/Densities/'

for i in range(len(File)):
    print File[i]
    dataset=Dataset(File[i])
    lat=dataset.variables['lat'][:]
    lon=dataset.variables['lon'][:]
    time=dataset.variables['time'][:]
    lon[lon>180]-=360
    #Now, we want the last 5 years of particle positions, since in this time
    #the garbage patch has largely been formed
    if lon.shape[1]==4748:
        lonLast=lon[:,-365*5:]
        latLast=lat[:,-365*5:]
        Time=time[:,-365*5:]
    else:
        lonLast=lon[:,-183*5:]
        latLast=lat[:,-183*5:]
        Time=time[:,-183*5:]
    density=np.zeros((lonLast.shape[1],360,160))
    for j in range(lonLast.shape[1]):
        density[j,:,:]=HistogramFunction(lonLast[:,j],latLast[:,j])
    density.dump(location+saveFiles[i])
    Time.dump(location+saveFiles[i]+'Time')

#%%
#File=['D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\AtlanticWindage_1.nc',
#      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\AtlanticWindage_3.nc',
#      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\AtlanticWindage_5.nc']
File=['D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles/Onink et al/AtlanticWindage0.01.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles/Onink et al/AtlanticWindage0.03.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles/Onink et al/AtlanticWindage0.05.nc']
saveFiles=['NorthAtlanticWindage1per','NorthAtlanticWindage3per','NorthAtlanticWindage5per']
location='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles/Onink et al/Densities/'
for i in range(len(File)):
    print File[i]
    dataset=Dataset(File[i])
    lat=dataset.variables['lat'][:]
    lon=dataset.variables['lon'][:]
    time=dataset.variables['time'][:]
    lon[lon>180]-=360
    #Now, we want the last 5 years of particle positions, since in this time
    #the garbage patch has largely been formed
    if lon.shape[1]==4748:
        lonLast=lon[:,-365*5:]
        latLast=lat[:,-365*5:]
        Time=time[:,-365*5:]
    else:
        lonLast=lon[:,-183*5:]
        latLast=lat[:,-183*5:]
        Time=time[:,-183*5:]
    density=np.zeros((lonLast.shape[1],360,160))
    for j in range(lonLast.shape[1]):
        density[j,:,:]=HistogramFunction(lonLast[:,j],latLast[:,j])
    density.dump(location+saveFiles[i])
    Time.dump(location+saveFiles[i]+'Time')
#%%
File=['D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles/Onink et al/AtlanticTotal3h_dt15m.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles/Onink et al/AtlanticTotal3h_dt30m.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal24hBeach.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal24h.nc']
saveFiles=['AtlanticIntegration_30m','AtlanticIntegration_15m','AtlanticBeach','AtlanticNoBeach']
location='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles/Onink et al/Densities/'
for i in [2,3]:#range(len(File)):
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
        density[j,:,:]=HistogramFunction(lonLast[:,j],latLast[:,j])
    density.dump(location+saveFiles[i])
#    Time.dump(location+saveFiles[i]+'Time')