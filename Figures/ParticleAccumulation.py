# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 15:39:19 2018

@author: Victor Onink
Figure showing the fraction of particles within the garbage patches in the North
Atlantic and the North Pacific
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.dates as mdates



#%%
NPlon=[-125,-155]
NPlat=[25,45]
NAlon=[-30,-75]
NAlat=[20,40]
color=['red','forestgreen','dodgerblue']
label=['Total','Ekman','Geostrophic']
#formatting the figures
axeslabelsize=14
textsize=12
fig=plt.figure(figsize=(10*1.5,8*1.5))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
yearsFmt = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_locator(years)
ax1.xaxis.set_major_formatter(yearsFmt)
ax1.xaxis.set_minor_locator(months)
ax2.xaxis.set_major_locator(years)
ax2.xaxis.set_major_formatter(yearsFmt)
ax2.xaxis.set_minor_locator(months)
datemin = datetime(2002, 1, 1)
datemax = datetime(2015, 1, 1)
ax1.set_xlim(datemin, datemax)
ax2.set_xlim(datemin, datemax)
#ax1.set_xlabel('Time (yr)',fontsize=axeslabelsize)
ax2.set_xlabel('Time (yr)',fontsize=axeslabelsize)
ax1.set_ylabel('Fraction (%)',fontsize=axeslabelsize)
ax2.set_ylabel('Fraction (%)',fontsize=axeslabelsize)

ax2.set_title('North Atlantic',fontsize=axeslabelsize,fontweight='bold')
ax1.set_title('North Pacific',fontsize=axeslabelsize,fontweight='bold')
ax1.set_ylim([0,80])
ax2.set_ylim([0,80])



File=['D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al/NorthPacificTotal24h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al/NorthPacificEkman.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al/NorthPacificGeostrophic.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al/AtlanticTotal24h.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al/NorthAtlanticEkman.nc',
      'D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al/AtlanticGeostrophic.nc'
      ]

for i in range(len(File)):
    dataset=Dataset(File[i])
    lon=dataset.variables['lon'][:]
    lon[lon>180]-=360
    lat=dataset.variables['lat'][:]
    time=dataset.variables['time'][:]
    time_axis=[]
    for z in range(lat.shape[1]):
        time_axis.append(datetime(2002,1,1,0)+timedelta(seconds=time[0,:][z]))
    
    patch,lonmatch,latmatch=np.zeros(lon.shape),np.zeros(lon.shape),np.zeros(lon.shape)
    if i==0 or i==1 or i==2: #Pacific
        lonmatch[(lon<NPlon[0])&(lon>NPlon[1])]=1
        latmatch[(lat>NPlat[0])&(lat<NPlat[1])]=1
        patch[(lonmatch==1)&(latmatch==1)]=1
        fraction=np.sum(patch,axis=0)/lon.shape[0]*100
        ax1.plot(time_axis,fraction,color[i%3],label=label[i%3])
    else:#Atlantic
        lonmatch[(lon<NAlon[0])&(lon>NAlon[1])]=1
        latmatch[(lat>NAlat[0])&(lat<NAlat[1])]=1
        patch[(lonmatch==1)&(latmatch==1)]=1
        fraction=np.sum(patch,axis=0)/lon.shape[0]*100
        ax2.plot(time_axis,fraction,color[i%3],label=label[i%3])
ax2.legend(fontsize=16)
plt.tight_layout()
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data\AtlanticPacificConvergencePatchExt.jpg')

