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
    binsLon=np.arange(0,360)
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
File=['D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\NorthPacificTotal3h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\NorthPacificTotal24h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\NorthPacificEkman.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\NorthPacificGeostrophic.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\NorthPacificStokes.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\NorthPacificStokeTotal3h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\NorthPacificStokeTotal.nc'
      ]
#File=['/scratch/Victor/NorthPacificTotal.nc',
#      '/scratch/Victor/NorthPacificEkmanNoBeach_05x05.nc',
#      '/scratch/Victor/NorthPacificGeostrophicNoBeach_05x05.nc',
#      '/scratch/Victor/NorthPacificStokeTotal.nc',
#      '/scratch/Victor/PacificStokeRunNoBeach1x1.nc']

saveFiles=['NorthPacificTotalDensity3h','NorthPacificTotalDensity24h',
           'NorthPacificEkmanDensity','NorthPacificGeostrophicDensity',
           'NorthPacificStokesDensity',
           'NorthPacificStokesTotalDensity3h','NorthPacificStokesTotalDensity24h'
           ]
location='D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\Densities/'
for i in range(len(File)):
    print File[i]
    dataset=Dataset(File[i])
    lat=dataset.variables['lat'][:]
    lon=dataset.variables['lon'][:]
    time=dataset.variables['time'][:]
    lon[lon<0]+=360
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
#    Time.dump(location+saveFiles[i]+'Time')