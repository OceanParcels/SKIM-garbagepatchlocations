# -*- coding: utf-8 -*-
"""
Created on Wed May 23 11:40:04 2018

@author: Victor Onink
"""
from netCDF4 import Dataset
from parcels import plotTrajectoriesFile,ParticleSet,JITParticle,FieldSet,Field
import numpy as np
from matplotlib import colors
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from datetime import datetime, timedelta
from matplotlib.patches import Polygon

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

def draw_screen_poly( lats, lons, m,color):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, edgecolor=color, linewidth=1.5,facecolor='none' )
    plt.gca().add_patch(poly)

#%% With boxes showing the garbage patch and extended garbage patch region
latD=np.linspace(-80,80,160)
lonD=np.linspace(-180,180,360)
NPlon,NPlonex=[-130,-130,-150,-150],[-125,-125,-155,-155]
NPlat,NPlatex=[30,40,40,30],[25,45,45,25]
NAlon,NAlonex=[-40,-40,-70,-70],[-35,-35,-75,-75]
NAlat,NAlatex=[25,35,35,25],[20,40,40,20]
lonPatch=[NAlon,NAlonex,NPlon,NPlonex]
latPatch=[NAlat,NAlatex,NPlat,NPlatex]
color=['red','black','red','black']

def plotDensity(typ,lon,lat,dens):
    Lat,Lon=np.meshgrid(lat,lon)
    latmin,latmax=-90,90
    lonmin,lonmax=-180,180
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    my_map.drawmeridians(np.arange(0, 360, 30),labels=[0,0,0,1],fontsize=11)
    if (typ+1)%2==1:
        my_map.drawparallels(np.arange(-90, 91, 30),labels=[1,0,0,0],fontsize=11)
    else:
        my_map.drawparallels(np.arange(-90, 91, 30),fontsize=11)
    density=my_map.contourf(Lon,Lat,dens/1e-9,np.linspace(1e-1,2e0,20),
                              #norm=colors.LogNorm(1e-10,1e-9),
                              cmap='rainbow',extend='both')
    title=['(a) Average Year 9','(b) Average Year 10','(c) Average Year 11','(d) Average Year 12']
    plt.title(title[typ],fontsize=14,fontweight='bold')    
    return density
#    cbar=my_map.colorbar(density)
#    cbar.ax.tick_params(labelsize=12)
#    cbar.set_label("Plastic Counts ($10^{-3}$ # km$^{-2}$)", rotation=90,fontsize=12)
location='D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al/Densities/'
File=['GlobalGeostrophicDensity','GlobalGeostrophicDensity',
      'GlobalGeostrophicDensity','GlobalGeostrophicDensity']
fig,axes=plt.subplots(nrows=2, ncols=2,figsize=(10*2,8*1.5))
for i in range(len(File)):
    density=np.load(location+File[i])
    density[np.isnan(density)]=0
    meanFinalYear=np.sum(density[183*(i+1):183*(i+2),:,:]/
                                 density[183*(i+1):183*(i+2),:,:].shape[0],axis=0)
    meanFinalYear[meanFinalYear==0]=np.nan
    plt.subplot(2,2,i+1)
    density=plotDensity(i,lonD,latD,meanFinalYear)
fig.subplots_adjust(right=0.90)
cbar_ax = fig.add_axes([0.93, 0.133, 0.02, 0.72])
cbar=fig.colorbar(density,cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_yticklabels(['<0.1','0.3','0.5','0.7','0.9','1.1','1.3','1.5','1.7','1.9<'])
cbar.set_label("Plastic Counts ($10^{-3}$ # km$^{-2}$)", rotation=90,fontsize=13)
plt.subplots_adjust(wspace=0.1)

plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\Global\Figures\GlobalGeostrophicStability.jpg',
            bbox_inches='tight')