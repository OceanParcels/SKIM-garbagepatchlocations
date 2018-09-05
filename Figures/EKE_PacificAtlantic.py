# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 11:47:47 2018

@author: Victor Onink
"""

from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates


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
def AreaWeighAverage(dataset,coarseness):
    AreaInitial=AreaCalc(dataset.shape[0],dataset.shape[1])
    AreaFinal=AreaCalc(dataset.shape[0]/coarseness,dataset.shape[1]/coarseness)
    dataset[np.isnan(dataset)]=0
    DataSetArea=np.multiply(AreaInitial,dataset)
    temp=DataSetArea.reshape((DataSetArea.shape[0] // coarseness,coarseness,
                              DataSetArea.shape[1] // coarseness,coarseness))
    finalDataset = np.sum(temp,axis=(1,3))
    finalDataset=np.divide(finalDataset,AreaFinal)
    finalDataset[finalDataset==0]=np.nan
    return finalDataset
def plotParticles(lon, lat):
    latmin,latmax=-5,75
    lonmin,lonmax=-180,-100
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
    my_map.drawcoastlines()
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 10),labels=[0,0,0,1])
    my_map.drawparallels(np.arange(-90, 90, 10),labels=[1,0,0,0])
    my_map.plot(lon,lat,'r.',markersize=3)
#%%
File='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticTotal24h.nc'
dataset=Dataset(File)
trajectory=dataset.variables['trajectory'][:]
lat=dataset.variables['lat'][:]
lon=dataset.variables['lon'][:]
EKE=dataset.variables['EKE'][:]
time=dataset.variables['time'][:]
#%% First we need to see which of the particles are within the garbage patch at the end
lastTraj=trajectory[:,-1]
lastLat=lat[:,-1]
lastLon=lon[:,-1]
#[lastLon1,lastLat1]=[np.array([lo for lo, la in zip(lastLon,lastLat) if  (25<la< 35 and -70<lo<-45)  ]),
#                    np.array([la for lo, la in zip(lastLon,lastLat) if (25<la< 35 and -70<lo<-45)])]    
#plt.figure(1)
#plotParticles(lastLon1,lastLat1)
#trashpatch=lastTraj[(25<lastLat<35) and (-70<lastLon<-45)]
trashpatch=[]
trashpatchXL=[]
outside=[]
for i in range(len(lastTraj)):
    if 25<lastLat[i]<35 and -70<lastLon[i]<-40:
        trashpatch.append(lastTraj[i])
    else:
        outside.append(lastTraj[i]) 
for i in range(len(lastTraj)):
    if 20<lastLat[i]<40 and -75<lastLon[i]<-35:
        trashpatchXL.append(lastTraj[i])
        
#Now for the x axis I am going to want to have time, so we need datetime objects for that
time_axis=[]
for i in range(lat.shape[1]):
    time_axis.append(datetime(2002,1,1,0)+timedelta(seconds=time[0,:][i]))

AmeanEKE=np.mean(EKE[trashpatch,1:],axis=0)
AmeanEKExl=np.mean(EKE[trashpatchXL,1:],axis=0)
AmeanEKEtot=np.mean(EKE[:,1:],axis=0)

#%% Now we repeat the whole thing for the pacific
File='D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\NorthPacificTotal24h.nc'
dataset=Dataset(File)
trajectory=dataset.variables['trajectory'][:]
lat=dataset.variables['lat'][:]
lon=dataset.variables['lon'][:]
EKE=dataset.variables['EKE'][:]
time=dataset.variables['time'][:]
#%% First we need to see which of the particles are within the garbage patch at the end
lastTraj=trajectory[:,-1]
lastLat=lat[:,-1]
lastLon=lon[:,-1]
#[lastLon1,lastLat1]=[np.array([lo for lo, la in zip(lastLon,lastLat) if  (30<la< 40 and -150<lo<-130)  ]),
#                    np.array([la for lo, la in zip(lastLon,lastLat) if (30<la< 40 and -150<lo<-130)])]    
#plt.figure(1)
#plotParticles(lastLon1,lastLat1)
#trashpatch=lastTraj[(25<lastLat<35) and (-70<lastLon<-45)]
trashpatch=[]
trashpatchXL=[]
outside=[]
for i in range(len(lastTraj)):
    if 30<lastLat[i]<40 and -150<lastLon[i]<-130:
        trashpatch.append(lastTraj[i])
    else:
        outside.append(lastTraj[i])
for i in range(len(lastTraj)):
    if 25<lastLat[i]<45 and -155<lastLon[i]<-125:
        trashpatchXL.append(lastTraj[i])
#Now for the x axis I am going to want to have time, so we need datetime objects for that
time_axis=[]
for i in range(lat.shape[1]):
    time_axis.append(datetime(2002,1,1,0)+timedelta(seconds=time[0,:][i]))

PmeanEKE=np.mean(EKE[trashpatch,1:],axis=0)
PmeanEKExl=np.mean(EKE[trashpatchXL,1:],axis=0)
PmeanEKEtot=np.mean(EKE[:,1:],axis=0)

#%% Now we need the actual figure
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')
fig=plt.figure(figsize=(10*1.5,8*1.5))
ax1 = fig.add_subplot(211)
ax1.semilogy(time_axis[1:],AmeanEKE,'r',label='End Within Garbage Patch')
ax1.semilogy(time_axis[1:],AmeanEKExl,'forestgreen',label='End Within Extended Garbage Patch')
ax1.semilogy(time_axis[1:],AmeanEKEtot,'b',label='All Particles')
#formatting!
ax1.xaxis.set_major_locator(years)
ax1.xaxis.set_major_formatter(yearsFmt)
ax1.xaxis.set_minor_locator(months)

axeslabelsize=14
textsize=12

datemin = datetime(2002, 1, 1)
datemax = datetime(2015, 1, 1)
ax1.set_ylabel('EKE (m$^2$ s$^{-1}$)',fontsize=axeslabelsize)
ax1.set_ylim([3e-4,6e-2])
ax1.set_xlim(datemin, datemax)
ax1.tick_params(labelsize=textsize)
ax1.set_xlabel('Time (yr)',fontsize=axeslabelsize)
ax1.set_title('North Atlantic',fontsize=axeslabelsize,fontweight='bold')
ax1.tick_params(which='major',length=7)
ax1.tick_params(which='minor',length=3)

ax2 = fig.add_subplot(212)
ax2.semilogy(time_axis[1:],PmeanEKE,'r',label='End Within Garbage Patch')
ax2.semilogy(time_axis[1:],PmeanEKExl,'forestgreen',label='End Within Extended Garbage Patch')
ax2.semilogy(time_axis[1:],PmeanEKEtot,'b',label='All Particles')
#formatting!
ax2.xaxis.set_major_locator(years)
ax2.xaxis.set_major_formatter(yearsFmt)
ax2.xaxis.set_minor_locator(months)

datemin = datetime(2002, 1, 1)
datemax = datetime(2015, 1, 1)
ax2.set_title('North Pacific',fontsize=axeslabelsize,fontweight='bold')
ax2.set_ylabel('EKE (m$^2$ s$^{-1}$)',fontsize=axeslabelsize)
ax2.set_ylim([3e-4,6e-2])
ax2.set_xlim(datemin, datemax)
ax2.tick_params(labelsize=textsize)
ax2.set_xlabel('Time (yr)',fontsize=axeslabelsize)
ax2.tick_params(which='major',length=7)
ax2.tick_params(which='minor',length=3)
ax2.legend(fontsize=16)
fig.autofmt_xdate()
plt.tight_layout()
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\AtlanticPacificEKE.jpg')

#%% Regression between timeseries
from scipy import stats
_, _, r_A, p_A, _ = stats.linregress(AmeanEKE,AmeanEKEtot)
_, _, r_Axl, p_Axl, _ = stats.linregress(AmeanEKExl,AmeanEKEtot)
_, _, r_AxlA, p_Axl, _ = stats.linregress(AmeanEKExl,AmeanEKE)
_, _, r_P, p_P, _ = stats.linregress(PmeanEKE,PmeanEKEtot)
_, _, r_Pxl, p_Pxl, _ = stats.linregress(PmeanEKExl,PmeanEKEtot)
_, _, r_PxlP, p_PxlP, _ = stats.linregress(PmeanEKExl,PmeanEKE)