# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 12:21:06 2018

@author: Victor Onink
Statistics of the basins of attraction. I want to be able to see changes in the origin
of particles statistically.
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
            #remove baltic issue
            
            #so, first north Atlantic
            if -77<finalLon[i]<45:
                Endspot=7
            if -85<finalLon[i]<-77:
                if finalLat[i]>9:
                    Endspot=7
                else:
                    Endspot=5
            if -100<finalLon[i]<-85:
                if finalLat[i]>16:
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
        if 50<=finalLat[i]<60:
            #North Atlantic
            if -80<finalLon[i]<100:
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
        if finalLat[i]<=-50:
            Endspot=2
        identifier[np.argmin(np.abs(firstLat[i]-latGrid)),np.argmin(np.abs(firstLon[i]-lonGrid))]=Endspot
    LonG,LatG=np.meshgrid(lonGrid,latGrid)
    identifier[identifier==0]=np.nan
    return [identifier,LonG,LatG]
#%%
"""
1=arctic, 2=southern,3=Indian,4=South Pacific,5=North Pacific,6=South Atlantic,7=North Atlantic
"""
File=['D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al\GlobalTotal24h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al\GlobalStokeTotal.nc']
k=1
dataset=Dataset(File[k])
lat=dataset.variables['lat'][:]
lon=dataset.variables['lon'][:]
lon[lon>180]-=360
lon[lon<-180]+=360
finalLat,finalLon=lat[:,-1],lon[:,-1]
finalLon[finalLon>180]-=360
firstLat,firstLon=lat[:,0],lon[:,0]
Finalgrid,lonG,latG=BasinShaper(finalLon,finalLat,firstLon,firstLat)
Initialgrid,_,_=BasinShaper(firstLon,firstLat,firstLon,firstLat)

##%%
#latmin,latmax=-90,90
#lonmin,lonmax=-180,180
#my_map = Basemap(projection='cyl', llcrnrlon=lonmin, 
#                  urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax, 
#                  resolution='l')
##draw the polygons to mark the different ocean basins
##for j in range(len(latCoordBasins)):
##    draw_screen_poly(latCoordBasins[j],lonCoordBasins[j],my_map)
##my_map.drawcoastlines()
#my_map.fillcontinents(color = 'gray')
#my_map.drawmapboundary()
#
#my_map.drawmeridians(np.arange(0, 360, 15),labels=[0,0,0,1],zorder=10,fontsize=11)
#
#my_map.drawparallels(np.arange(-90, 91, 15),labels=[1,0,0,0],zorder=10,fontsize=11)
##for k in range(len(lonText)):
##    plt.text(lonText[k],latText[k],textlabel[k],fontsize=14,weight='bold')
##colors=['red','darkorange','forestgreen','lightgreen','gold','cyan']
#colors=['deepskyblue','cyan','gold','lightgreen','forestgreen','darkorange','red']
#cmap_c=c.ListedColormap(colors)
#label=['Free','North Atlantic','North Pacific','South Pacific','South Atlantic','Indian Ocean','Southern Ocean','Arctic Ocean']
#basin=my_map.pcolormesh(lonG,latG,Initialgrid,cmap=cmap_c)
#title=['Total Currents','Total + Stokes Currents']
#plt.title(title[i],fontweight='bold',fontsize=14)
#%%
"""
Now, for each of the basins, I want to know where the particles originate from, with 
each basin being given as a percentage of the total amount of particles within that
basin at the end
"""
fractions=[]
for i in range(1,8):
    flatFinal,flatInitial=Finalgrid.flatten(),Initialgrid.flatten()
    basin=flatFinal==i
    N=len(basin[basin==True])
    fractionssub=[N]
    for j in range(1,8):
        fractionssub.append(len(flatInitial[basin==True][flatInitial[basin==True]==j])*100./N)
    fractions.append(fractionssub)
    