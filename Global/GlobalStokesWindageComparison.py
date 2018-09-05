# -*- coding: utf-8 -*-
"""
Created on Wed Jul 04 10:21:55 2018

@author: Victor Onink
Here we create the plots showing the RMSE error of the windage and correlations with the individual components
"""
from netCDF4 import Dataset
import numpy as np
from matplotlib import colors
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from datetime import datetime, timedelta
import scipy
from matplotlib.ticker import LogFormatter,LogLocator 

#%%
FileRMS='D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al/RMSE_WindageStokesMagnitudeUV02_14.nc'
Filecor='D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al/CorrelationsWindageStokesMagnitudeUV02_14.nc'


datasetRMS=Dataset(FileRMS)
RMS=datasetRMS.variables['Field'][:]
lon=datasetRMS.variables['longitude'][:]
lat=datasetRMS.variables['latitude'][:]
Lon,Lat=np.meshgrid(lon,lat)
datasetCOR=Dataset(Filecor)
Cor=datasetCOR.variables['Field'][:]
#%%
latmin,latmax=-90,90
lonmin,lonmax=-180,180
axeslabelsize=14
textsize=12
fig=plt.figure(figsize=(10*1,8*1.2))
for i in range(2):
    if i==0:
        ax1 = fig.add_subplot(211)    
    if i==1:
        ax2= fig.add_subplot(212)
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
    #my_map.drawcoastlines()
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    if i==0:
        my_map.drawmeridians(np.arange(0, 360, 30),zorder=10,fontsize=textsize)
    else:
        my_map.drawmeridians(np.arange(0, 360, 30),labels=[0,0,0,1],zorder=10,fontsize=textsize)
    my_map.drawparallels(np.arange(-90, 91, 30),labels=[1,0,0,0],zorder=10,fontsize=textsize)
    correlation=my_map.pcolormesh(Lon,Lat,np.square(Cor[i,:,:]),
                      vmin=0.2,vmax=1,
                      cmap='rainbow')
ax1.set_title('Zonal Velocity U',fontsize=axeslabelsize,fontweight='bold')
ax2.set_title('Meridional Velocity V',fontsize=axeslabelsize,fontweight='bold')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.11, 0.035, 0.77])
cbar=fig.colorbar(correlation,cax=cbar_ax)
cbar.ax.tick_params(labelsize=textsize)
cbar.ax.set_yticklabels(['<0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'])
cbar.set_label(r"Coefficient of Determination (r$^{2}$)", rotation=90,fontsize=axeslabelsize)
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\Global\Figures\GlobalWindageStokesCorrelation.jpg',
            bbox_inches='tight')

#%%
def plotDensityRMS(typ,lon,lat,correlation):
    Lon,Lat=np.meshgrid(lon,lat)
    latmin,latmax=-90,90
    lonmin,lonmax=-180,180
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 30),labels=[0,0,0,1],fontsize=11)
    my_map.drawparallels(np.arange(-90, 91, 30),labels=[1,0,0,0],fontsize=11)
    sizeRMS=plt.pcolormesh(Lon,Lat,correlation,zorder=1,cmap='rainbow',
                        vmin=0.01,vmax=0.08)
    return sizeRMS
fig=plt.figure(1,figsize=(10*1.5,8*1.5))
ax1=fig.add_subplot(111)
plotDensityRMS(0,lon,lat,RMS[0,:,:])
sizeRMS=plotDensityRMS(0,lon,lat,RMS[0,:,:])
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.93, 0.25, 0.02, 0.49])
cbar=fig.colorbar(sizeRMS,cax=cbar_ax)
cbar.ax.set_yticklabels(['<0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08<'])
cbar.ax.tick_params(labelsize=12)
cbar.set_label("RMSE (m s$^{-1}$)", rotation=90,fontsize=13)
ax1.set_title('RMSE 1% Windage', fontsize=14, fontweight='heavy')
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\Global\Figures\GlobalWindage1perRMSE.jpg',
            bbox_inches='tight')

#%% I want to scale the RMS by the mean Stokes drift, so I can in relative terms where
# the biggest deviations are.
def plotDensityRMSscale(typ,lon,lat,correlation):
    Lon,Lat=np.meshgrid(lon,lat)
    latmin,latmax=-90,90
    lonmin,lonmax=-180,180
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 30),labels=[0,0,0,1],fontsize=11)
    my_map.drawparallels(np.arange(-90, 90, 30),labels=[1,0,0,0],fontsize=11)
    sizeRMS=plt.pcolormesh(Lon,Lat,correlation,zorder=1,cmap='rainbow',
                        vmin=20,vmax=200)
    return sizeRMS
fig=plt.figure(1,figsize=(10*1.5,8*1.5))
ax1=fig.add_subplot(111)
plotDensityRMSscale(0,lon,lat,np.divide(RMS[0,:,:],RMS[3,:,:])*100)
sizeRMS=plotDensityRMSscale(0,lon,lat,np.divide(RMS[0,:,:],RMS[3,:,:])*100)
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.93, 0.25, 0.02, 0.49])
cbar=fig.colorbar(sizeRMS,cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label("Percentage Difference (%)", rotation=90,fontsize=13)
ax1.set_title('RMSE 1% Windage', fontsize=14, fontweight='heavy')
#plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\Global\Figures\GlobalWindage1perRMSEscale.jpg',
#            bbox_inches='tight')
