# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 13:35:23 2018

@author: Victor Onink
Here we create a figure that has the 24h and the 3h densities
for the North Atlantic
"""

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy import io
import pandas as pd

def AreaCalc(sizeLat,sizeLon): #Calculate surface area of grid cells
    deg2rd = np.pi/180.
    r=6378.1
    lon_bins = np.linspace(0,360., sizeLon+1)
    lat_bins = np.linspace(-90, 90, sizeLat+1) 
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
    lonmin,lonmax=-100,20
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
    my_map.drawcoastlines()
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 30))
    my_map.drawparallels(np.arange(-90, 90, 30))
    my_map.plot(lon,lat,'r.',markersize=3)
    
def HistogramFunction(londata,latdata):
    londata,latdata=londata.reshape(np.size(londata)),latdata.reshape(np.size(latdata))
    binsLon=np.arange(-180,180)
    binsLat=np.arange(-90,90)
    density=np.zeros((len(binsLon),len(binsLat)))
    for i in range(np.array(londata).shape[0]):
        density[np.argmin(np.abs(londata[i]-binsLon)),np.argmin(np.abs(latdata[i]-binsLat))]+=1
    #Now, normalize it by area
    area=AreaCalc(len(binsLat),len(binsLon)).T
    density/=area
    density[density==0]=np.nan
    return density

def plotDensity(typ,lon,lat,dens):
    Lat,Lon=np.meshgrid(lat,lon)
#    Lon,Lat=np.meshgrid(lon,lat)
    latmin,latmax=-5,75
    lonmin,lonmax=-100,20
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
#    my_map.drawcoastlines()
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 30),labels=[0,0,0,1],fontsize=12)
    if (typ+1)%2==1:
        my_map.drawparallels(np.arange(-90, 91, 30),labels=[1,0,0,0],fontsize=12)
    else:
        my_map.drawparallels(np.arange(-90, 91, 30),fontsize=12)
    density=my_map.contourf(Lon,Lat,dens/1e-9,np.linspace(1e-1,4e0,40),
                              #norm=colors.LogNorm(1e-10,1e-9),
                              cmap='rainbow',extend='both')
#    my_map.contour(LonSam,LatSam,sample,[1e4,1e5],colors='k',zorder=100)
#    my_map.contour(LonP,LatP,CountPvS,[5e4],colors='k',zorder=100)
    title=['(a) Total Currents ($\Delta$t=24h)','(b) Total Currents + Stokes Drift ($\Delta$t=24h)',
           '(c) Total Currents ($\Delta$t=3h)','(d) Total Currents + Stokes Drift ($\Delta$t=3h)']
    plt.title(title[typ],fontsize=14,fontweight='bold')
    return density
#    cbar=my_map.colorbar(density)
#    cbar.ax.tick_params(labelsize=12)
#    cbar.set_label("Plastic Counts ($10^{-3}$ # km$^{-2}$)", rotation=90,fontsize=12)

#%%                                    
location='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\Densities/'
File=['NorthAtlanticTotalDensity24h','NorthAtlanticStokesTotal24h',
      'NorthAtlanticTotalDensity3h','NorthAtlanticStokesTotalDensity3h']
axeslabelsize=14
textsize=12
fig,axes=plt.subplots(nrows=2, ncols=1,figsize=(10*1.55,8*1.5))
for i in range(len(File)):
    density=np.load(location+File[i])
    density[np.isnan(density)]=0
    meanFinalYear=np.sum(density[-183:,:,:]/density[-183:,:,:].shape[0],axis=0)
    meanFinalYear[meanFinalYear==0]=np.nan
    latD=np.linspace(-80,80,160)
    lonD=np.linspace(-180,180,360)
    plt.subplot(2,2,i+1)
    density=plotDensity(i,lonD,latD,meanFinalYear)
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.93, 0.13, 0.02, 0.73])
cbar=fig.colorbar(density,cax=cbar_ax)
cbar.ax.tick_params(labelsize=textsize)
cbar.set_label("Plastic Counts ($10^{-3}$ # km$^{-2}$)", rotation=90,fontsize=axeslabelsize)
cbar.ax.set_yticklabels(['<0.1','0.5','0.9','1.3','1.7','2.1','2.5','2.9','3.3','3.7<'])
plt.subplots_adjust(wspace=0.1)
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\Figures\NorthAtlanticTimeStepDensities.jpg',
            bbox_inches='tight') 

#%% Subplot Comparing the two integration timesteps
def plotDensity(typ,lon,lat,dens):
    Lat,Lon=np.meshgrid(lat,lon)
#    Lon,Lat=np.meshgrid(lon,lat)
    latmin,latmax=-5,75
    lonmin,lonmax=-100,20
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
#    my_map.drawcoastlines()
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 30),labels=[0,0,0,1],fontsize=22)
    my_map.drawparallels(np.arange(-90, 90, 30),labels=[1,0,0,0],fontsize=22)
    density=my_map.contourf(Lon,Lat,dens/1e-9,np.linspace(1e-1,4e0,40),
                              #norm=colors.LogNorm(1e-10,1e-9),
                              cmap='rainbow',extend='both')
#    my_map.contour(LonSam,LatSam,sample,[1e4,1e5],colors='k',zorder=100)
#    my_map.contour(LonP,LatP,CountPvS,[5e4],colors='k',zorder=100)
    title=['Integration dt=30 min','Integration dt=15 min']
    plt.title(title[typ],fontsize=24,fontweight='bold')
    return density
#    cbar=my_map.colorbar(density)
#    cbar.ax.tick_params(labelsize=12)
#    cbar.set_label("Plastic Counts ($10^{-3}$ # km$^{-2}$)", rotation=90,fontsize=12)
                                    
location='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al/Densities/'
File=['AtlanticIntegration_30m','AtlanticIntegration_15m']
fig,axes=plt.subplots(nrows=2, ncols=1,figsize=(8*1.5,8*1.5))
for i in range(len(File)):
    density=np.load(location+File[i])
    density[np.isnan(density)]=0
    meanFinalYear=np.sum(density[-183:,:,:]/density[-183:,:,:].shape[0],axis=0)
    meanFinalYear[meanFinalYear==0]=np.nan
    latD=np.linspace(-80,80,160)
    lonD=np.linspace(-180,180,360)
    plt.subplot(2,1,i+1)
    density=plotDensity(i,lonD,latD,meanFinalYear)
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.85, 0.06, 0.02, 0.86])
cbar=fig.colorbar(density,cax=cbar_ax)
cbar.ax.tick_params(labelsize=22)
cbar.set_label("Plastic Counts ($10^{-3}$ # km$^{-2}$)", rotation=90,fontsize=24)
plt.tight_layout(pad=5)
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\Figures\NorthAtlanticIntegrationTimeStep.jpg')
