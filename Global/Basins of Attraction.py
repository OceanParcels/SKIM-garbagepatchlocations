# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:15:13 2018

@author: Victor Onink
We are going to do basins of attraction, wohoooo
Basically, color the particles according to the end location where they end up
First, we will do the total currents
"""
from netCDF4 import Dataset
from parcels import plotTrajectoriesFile,ParticleSet,JITParticle,FieldSet
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from matplotlib import colors as c
from datetime import datetime, timedelta
from matplotlib.patches import Polygon

def BasinShaper(finalLon,finalLat,firstLon,firstLat):
    lonGrid=np.linspace(-180,179,360)
    latGrid=np.linspace(-90,90,181)
    identifier=np.zeros((len(latGrid),len(lonGrid)))
    for i in range(len(lon)):
        if 0<finalLat[i]<50:
            #so, first north Atlantic
            if -80<finalLon[i]<30:
                Endspot=7
            if -100<finalLon[i]<-80:
                if finalLat[i]>10:
                    Endspot=7
                else:
                    Endspot=5
            #North Pacific
            if -100>finalLon[i]>-180:
                Endspot=5
            if 180>finalLon[i]>100:
                Endspot=5
            #Indian Ocean
            if 40<finalLon[i]<100:
                Endspot=3
            #Red Sea, which we won't plot
            if 10<finalLat[i]<30:
                if 20<finalLon[i]<43:
                    Endspot=np.nan
        if 50<finalLat[i]<60:
            #North Atlantic
            if -80<finalLon[i]<80:
                Endspot=7
        	#North Pacific
            if -100>finalLon[i]>-180:
                Endspot=5
            if 180>finalLon[i]>100:
                Endspot=5
        if 60<finalLat[i]:
            if -120<finalLon[i]<90:
                Endspot=1 #arctic Ocean
        if -50<finalLat[i]<0:
            #South Atlantic
            if 20>finalLon[i]>-70:
                Endspot=6
            #South Pacific
            if -70>finalLon[i]>-180:
                Endspot=4
            if 180>finalLon[i]>150:
                Endspot=4
            #indian Ocean
            if 150>finalLon[i]>20:
                Endspot=3
            #southern ocean
        if finalLat[i]<-50:
            Endspot=2
        identifier[np.argmin(np.abs(firstLat[i]-latGrid)),np.argmin(np.abs(firstLon[i]-lonGrid))]=Endspot
    LonG,LatG=np.meshgrid(lonGrid,latGrid)
    identifier[identifier==0]=np.nan
    return [identifier,LonG,LatG]
def draw_screen_poly( lats, lons, m):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, edgecolor='black', linewidth=2,facecolor='none' )
    plt.gca().add_patch(poly)
#%%
File=['D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al\GlobalTotal24h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al\GlobalStokeTotal.nc']
fig,axes=plt.subplots(nrows=2, ncols=1,figsize=(10*1.2,8*1.5))
latCoordBasins=[[0,-50,-50,0],#South Atlantic
                [0,0,-50,-50],#south pacific East
                [-30,0,0,-50,-50,-30],#south pacific West
                [-50,-50,-30,0,0,15,30,30,30],#indian ocean
                [-50,-50],#Southern Ocean
                [60,60]#arctic ocean
                ]
lonCoordBasins=[[20,20,-69,-69],#South Atlantic
                [-69,-180,-180,-69],#south pacific East
                [150,120,180,180,150,150],#south pacific West
                [20,150,150,120,100,100,105,105,20],#indian ocean
                [-180,20],#Southern Ocean
                [-120,90]#arctic ocean
                ]
textlabel=['1','2','3','3','4','4','5','6','7']
lonText=[-50,-35,-175,151,-175,151,45,-175,-15]
latText=[50,-10,40,40,-10,-10,-10,-60,67]
for i in range(len(File)):
    dataset=Dataset(File[i])
    lat=dataset.variables['lat'][:]
    lon=dataset.variables['lon'][:]
    lon[lon>180]-=360
    lon[lon<-180]+=360
    finalLat,finalLon=lat[:,-1],lon[:,-1]
    finalLon[finalLon>180]-=360
    firstLat,firstLon=lat[:,0],lon[:,0]
    grid,lonG,latG=BasinShaper(finalLon,finalLat,firstLon,firstLat)
    latmin,latmax=-90,90
    lonmin,lonmax=-180,180
    plt.subplot(2,1,i+1)
    my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
                      urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
                      resolution='l')
    #draw the polygons to mark the different ocean basins
    for j in range(len(latCoordBasins)):
        draw_screen_poly(latCoordBasins[j],lonCoordBasins[j],my_map)
    #my_map.drawcoastlines()
    my_map.fillcontinents(color = 'gray')
    my_map.drawmapboundary()
    if i==1:
        my_map.drawmeridians(np.arange(0, 360, 30),labels=[0,0,0,1],zorder=10,fontsize=11)
    else:
        my_map.drawmeridians(np.arange(0, 360, 30),zorder=10)
    my_map.drawparallels(np.arange(-90, 91, 30),labels=[1,0,0,0],zorder=10,fontsize=11)
    for k in range(len(lonText)):
        plt.text(lonText[k],latText[k],textlabel[k],fontsize=14,weight='bold')
    #colors=['red','darkorange','forestgreen','lightgreen','gold','cyan']
    colors=['deepskyblue','mediumblue','gold','lightgreen','forestgreen','darkorange','red']
    cmap_c=c.ListedColormap(colors)
    label=['Free','North Atlantic','North Pacific','South Pacific','South Atlantic','Indian Ocean','Southern Ocean','Arctic Ocean']
    basin=my_map.pcolormesh(lonG,latG,grid,cmap=cmap_c)
    title=['(a) Total Currents','(b) Total Currents +  Stokes Drift']
    plt.title(title[i],fontweight='bold',fontsize=14)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.115, 0.035, 0.76])
cbar=fig.colorbar(basin,cax=cbar_ax)
cbar.set_ticks(np.linspace(1.4,6.6,7))
cbar.ax.set_yticklabels(['Arctic (7)','Southern (6)', 'Indian (5)', 'South Pacific (4)','North Pacific (3)','South Atlantic (2)','North Atlantic (1)'])
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\Global\Figures\BasinsOfAttractionLongitudeTotalStokesCurrentsGrids.png',
            bbox_inches='tight')