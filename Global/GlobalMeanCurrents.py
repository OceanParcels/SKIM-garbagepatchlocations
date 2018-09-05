# -*- coding: utf-8 -*-
"""
Created on Wed Jul 04 10:21:55 2018

@author: Victor Onink
Global mean currents
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

##%% Loading all the data!
folder='D:\Desktop\Thesis\ParcelsFigData\Data\Average Currents/'
#dataset=Dataset(folder+'AverageEkmanCurrents.nc')
#latEk=dataset.variables['lat'][:]
#lonEk=dataset.variables['lon'][:]
#lonEk[lonEk>180]-=360
#uEk=dataset.variables['Uavg'][:]
#vEk=dataset.variables['Vavg'][:]
#magEk=np.sqrt(np.square(uEk)+np.square(vEk))
#
#dataset=Dataset(folder+'AverageGeostrophicCurrents.nc')
#latGeo=dataset.variables['lat'][:]
#lonGeo=lonEk
##lonGeo[lonGeo>180]-=360
#uGeo=dataset.variables['Uavg'][:]
#uGeo=uGeo.reshape(720,1440)
#uGeo=np.hstack((uGeo[:,720:],uGeo[:,:720]))
#vGeo=dataset.variables['Vavg'][:]
#vGeo=vGeo.reshape(720,1440)
#vGeo=np.hstack((vGeo[:,720:],vGeo[:,:720]))
#magGeo=np.sqrt(np.square(uGeo)+np.square(vGeo))
#
#dataset=Dataset(folder+'AverageTotalCurrentsNoTime.nc')
#latTot=dataset.variables['lat'][:]
#lonTot=dataset.variables['lon'][:]
##lonTot[lonTot>180]-=360
#uTot=dataset.variables['Uavg'][:]
#vTot=dataset.variables['Vavg'][:]
#magTot=np.sqrt(np.square(uTot)+np.square(vTot))

#%%
fileTot='D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al/AverageTotalCurrents.nc'
fileEk='D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al/AverageEkmanCurrents.nc'
fileGeo='D:\Desktop\Thesis\ParcelsFigData\Data\Global\OutputFiles\Onink et al/AverageGeostrophicCurrents.nc'

datasetTot=Dataset(fileTot)
latTot=datasetTot.variables['latitude'][:]
lonTot=datasetTot.variables['longitude'][:]
lonTot[lonTot>180]-=360
uTot=datasetTot.variables['Uavg'][:]
vTot=datasetTot.variables['Vavg'][:]
magTot=np.sqrt(np.square(uTot)+np.square(vTot))
magTot[magTot==0]=np.nan

datasetEk=Dataset(fileEk)
latEk=datasetEk.variables['latitude'][:]
lonEk=datasetEk.variables['longitude'][:]
lonEk[lonEk>180]-=360
uEk=datasetEk.variables['Uavg'][:]
vEk=datasetEk.variables['Vavg'][:]
magEk=np.sqrt(np.square(uEk)+np.square(vEk))
magEk[magEk==0]=np.nan

datasetTot=Dataset(fileGeo)
latGeo=datasetTot.variables['latitude'][:]
lonGeo=datasetTot.variables['longitude'][:]
lonGeo[lonGeo>180]-=360
uGeo=datasetTot.variables['Uavg'][:]
vGeo=datasetTot.variables['Vavg'][:]
magGeo=np.sqrt(np.square(uGeo)+np.square(vGeo))
magGeo[magGeo==0]=np.nan

dataset=scipy.io.loadmat(folder+'StokesDriftMean.mat')
lonS=np.arange(-180,180,.5)
latS=dataset['latStoke']
uS=dataset['avgUss']
uS=np.hstack((uS[:,360:],uS[:,:360]))
vS=dataset['avgVss']
vS=np.hstack((vS[:,360:],vS[:,:360]))
magS=np.sqrt(np.square(uS)+np.square(vS))

#%% Now we plot everything
def plotDensity(typ,lon,lat,u,v,mag):
    Lon,Lat=np.meshgrid(lon,lat)
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
    normU,normV=np.divide(u,mag),np.divide(v,mag)
    if typ!=3:
        step=20
    else:
        step=10
    plt.quiver(Lon[::step,::step],Lat[::step,::step],normU[::step,::step],normV[::step,::step],
                      width=1e-3,color='k',zorder=5)
    size=plt.pcolormesh(Lon,Lat,mag,zorder=1,cmap='rainbow',
                        norm=colors.LogNorm(1e-2,1))
#                        vmin=0,vmax=.6)
    title=['(a) Total Currents','(b) Ekman Currents','(c) Geostrophic Currents','(d) Stokes Drift']
    plt.title(title[typ],fontsize=14,fontweight='bold')
    return size
fig,axes=plt.subplots(nrows=2, ncols=2,figsize=(10*2,8*1.5))
plt.subplot(2,2,1)
size=plotDensity(0,lonTot,latTot,uTot[0,:,:],vTot[0,:,:],magTot[0,:,:])
plt.subplot(2,2,2)
plotDensity(1,lonEk,latEk,uEk[0,:,:],vEk[0,:,:],magEk[0,:,:])
plt.subplot(2,2,3)
plotDensity(2,lonGeo,latGeo,uGeo[0,:,:],vGeo[0,:,:],magGeo[0,:,:])
plt.subplot(2,2,4)
plotDensity(3,lonS,latS,uS,vS,magS)
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.93, 0.14, 0.02, 0.715])
formatter = LogFormatter(10, labelOnlyBase=True) 
cbar=fig.colorbar(size,cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
labels=['<0.01','0.1','1<']
#labelsM=["%.2f" % z for z in labels]
actualLabels=[]
k=0
for i in range(20):
    if i%9==0:
        actualLabels.append(labels[k])
        k+=1
    else:
        actualLabels.append('')
cbar.ax.set_yticklabels(actualLabels)
cbar.set_label("Mean Current Velocity (m s$^{-1}$)", rotation=90,fontsize=13)
plt.subplots_adjust(wspace=0.1)
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\Global\Figures\MeanGlobalCurrents.jpg',
            bbox_inches='tight')
