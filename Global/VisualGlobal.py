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
#File='D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\GlobalStokeTotalRunNoBeach1x1.nc'
File='D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al\GlobalStokeTotal.nc'
dataset=Dataset(File)
trajectory=dataset.variables['trajectory'][:]
timeM=dataset.variables['time'][:]
lat=dataset.variables['lat'][:]
lon=dataset.variables['lon'][:]
lon[lon>180]-=360
#%%
latmin,latmax=-90,90
lonmin,lonmax=-180,180
plt.figure(figsize=(10*1.5,8*1.05))
my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                  urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                  resolution='l')
#my_map.drawcoastlines()
my_map.fillcontinents(color = 'gray')
my_map.drawmapboundary()
plt.tight_layout(pad=0)
#my_map.drawmeridians(np.arange(0, 360, 30))
#my_map.drawparallels(np.arange(-90, 90, 30))
plt.title('Ekman + Geostrophic Currents + Stokes Drift',weight='bold',fontsize=14)
text=plt.text(-176.3, -84.5,'hi',
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
    Time=datetime(2002,1,1)+timedelta(seconds=timeM[822,i])
    x, y = my_map(lons, lats)
    point.set_data(x, y)
    text.set_text('Date='+Time.strftime('%d-%m-%y'))
#    text.set_text('Date='+Time.strftime('%m-%Y'))
    return point,text
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(plt.gcf(), animate, init_func=init,
                               frames=lon.shape[1], interval=20, blit=True)
anim.save('D:\Desktop\Thesis\ParcelsFigData\Data\Global\Animations/GlobalStokeTotal.mov', 
          fps=50, extra_args=['-vcodec', 'libx264'])
plt.show()
