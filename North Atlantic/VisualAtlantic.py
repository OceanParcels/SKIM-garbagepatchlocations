# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:08:27 2018

@author: Victor Onink
We are going to see what we can get out of the output files
"""
from netCDF4 import Dataset
from parcels import plotTrajectoriesFile,ParticleSet,JITParticle,FieldSet
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from datetime import datetime, timedelta

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
#%%
#lon=np.load('D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\InputArray/LonsEastTestgrid0_3.npy')
#lat=np.load('D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\InputArray/LatsWestTestgrid0_3.npy')
#filenames = {'U': "D:\Desktop\Thesis\Data sets\GlobCurrent/2005/AllTogether/20050101*.nc",
#             'V': "D:\Desktop\Thesis\Data sets\GlobCurrent/2005/AllTogether/20050101*.nc"}
#variables = {'U': 'eastward_eulerian_current_velocity',
#             'V': 'northward_eulerian_current_velocity'}
#dimensions = {'lat': 'lat',
#              'lon': 'lon',
#              'time': 'time'}
#fieldset = FieldSet.from_netcdf(filenames, variables, dimensions,allow_time_extrapolation=True)

#File='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\AtlanticGeoTestRun.nc'
#File='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\AtlanticStokeTotalRunNoBeach05x05.nc'
#File='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\AtlanticStokeTotalRunNoBeach05x05.nc'
#File='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\AtlanticRun05x05.nc'
#File='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\NorthAtlanticStokeRunNoBeach1x1.nc'
File='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\AtlanticWindage0.05.nc'

dataset=Dataset(File)
trajectory=dataset.variables['trajectory'][:]
timeM=dataset.variables['time'][:]
lat=dataset.variables['lat'][:]
lon=dataset.variables['lon'][:]
lon[lon>180]-=360
#distance=dataset.variables['distance'][:]
#EKE=dataset.variables['EKE'][:]
#age=dataset.variables['Age'][:]
#plotTrajectoriesFile(File,mode='2d')
#pset=ParticleSet(fieldset=fieldset,pclass=JITParticle,lon=lon[:,-1],lat=lat[:,-1])
#density=AreaWeighAverage(pset.density(),4)

##Now we try to create basemaps
##set the background of the animation
#m=Basemap(projection='cyl', llcrnrlon=lonmin, 
#                  urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
#                  resolution='c')
#m.drawmapboundary()
#m.drawcoastlines()
#m.fillcontinents(color='grey')
#scat=m.plot([],[])
##scat=m.plot(lon[:,1],lat[:,1],'r.')
#def init():
#    del scat
#
#def animate(i):
#    longitude=lon[:,i]
#    latitude=lat[:,i]
#    scat=m.plot(longitude,latitude,'r.')
#    return scat
##Now to actually run the animation
#animation.FuncAnimation(plt.gcf(), animate,
#                               frames=20, interval=20, blit=True)
#
###anim.save('D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\Animations/Trial2005.mp4', 
###          fps=10, extra_args=['-vcodec', 'libx264'])
#Pick the dimensions of the basin we want to cover
#%%
latmin,latmax=-5,75
lonmin,lonmax=-100,20
plt.figure(figsize=(10*1.5,8*1.5))
my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                  urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                  resolution='l')
#my_map.drawcoastlines()
my_map.fillcontinents(color = 'gray')
my_map.drawmapboundary()
#my_map.drawmeridians(np.arange(0, 360, 30))
#my_map.drawparallels(np.arange(-90, 90, 30))
plt.title('Stokes Drift', fontsize=14,weight='bold')
plt.tight_layout()
#text=plt.annotate('hi', xy=(0.5, 0.9), xycoords='axes fraction')
text=plt.text(-98.75, -3.1,'',
                     ha='center',va='center',fontsize=12,
                     bbox={'facecolor':'white', 'alpha':1,'pad':10},zorder=200,
                     horizontalalignment='left')

#x,y=my_map(0.5,0.5)
#text=plt.text(x,y,'hi')
x,y = my_map(0, 0)
point = my_map.plot(x, y, 'r.', markersize=3)[0]
def init():
    point.set_data([], [])
    text.set_text('')
#    plt.title('blah')
    return point,#,text

# animation function.  This is called sequentially
def animate(i):
    lons, lats =  lon[:,i],lat[:,i]
    Time=datetime(2002,1,1)+timedelta(seconds=timeM[8474,i])
    x, y = my_map(lons, lats)
    point.set_data(x, y)
    text.set_text('Date='+Time.strftime('%d-%m-%y'))
#    text.set_text('Date='+Time.strftime('%m-%Y'))
    return point,text
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(plt.gcf(), animate, init_func=init,
                               frames=lon.shape[1], interval=60, blit=True)
#anim.save('D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\Animations/StokesNorthAtlantic.mp4', 
#          fps=40, extra_args=['-vcodec', 'libx264'])
#anim.save('D:\Desktop\Thesis\Azores Cruise\TestSimulation.mp4', 
#          fps=50, extra_args=['-vcodec', 'libx264'])
plt.show()
