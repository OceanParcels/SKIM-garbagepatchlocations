# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 09:38:54 2018

@author: Victor Onink
"""

import numpy as np
from matplotlib import colors
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
from datetime import datetime, timedelta
from scipy import io
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
import scipy.integrate as integrate

def BasinShaperAtlantic(inputDistribution,inputLon,inputLat):
    Endspot=np.zeros(inputDistribution.shape)
    for i in range(len(inputLon)):
        
        for j in range(len(inputLat)):
            if 0<inputLat[j]<50:
                #so, first north Atlantic
                if -80<inputLon[i]<30:
                    Endspot[j,i]=7
                if -100<inputLon[i]<-80:
                    if inputLat[j]>10:
                        Endspot[j,i]=7
                    else:
                        Endspot[j,i]=np.nan
                #North Pacific
                elif -100>inputLon[i]>-180:
                    Endspot[j,i]=np.nan
                if 180>inputLon[i]>100:
                    Endspot[j,i]=np.nan
                #Indian Ocean
                if 40<inputLon[i]<100:
                    Endspot[j,i]=np.nan
                #Red Sea, which we won't plot
                if 10<inputLat[j]<30:
                    if 20<inputLon[i]<43:
                        Endspot[j,i]=np.nan
            if 50<inputLat[j]<60:
                #North Atlantic
                if -80<inputLon[i]<80:
                    Endspot[j,i]=7
            	#North Pacific
                if -100>inputLon[i]>-180:
                    Endspot[j,i]=np.nan
                if 180>inputLon[i]>100:
                    Endspot[j,i]=np.nan
            if 60<inputLat[j]:
                Endspot[j,i]=np.nan
            if inputLat[j]<0:
                Endspot[j,i]=np.nan
    return Endspot
def BasinShaperPacific(inputDistribution,inputLon,inputLat):
    Endspot=np.zeros(inputDistribution.shape)
    for i in range(inputLon.shape[1]):
        for j in range(inputLat.shape[1]):
            if 0<inputLat[0,j]<50:
                if 260<inputLon[0,i]<280:
                    if inputLat[0,j]>10:
                        Endspot[j,i]=np.nan
                    else:
                        Endspot[j,i]=5
                #North Pacific
                elif 260>inputLon[0,i]>180:
                    Endspot[j,i]=5
                elif 180>inputLon[0,i]>100:
                    Endspot[j,i]=5
                else:
                    Endspot[j,i]=np.nan
            elif 50<inputLat[0,j]<60:
                if 260>inputLon[0,i]>180:
                    Endspot[j,i]=5
                elif 180>inputLon[0,i]>100:
                    Endspot[j,i]=5
                else:
                    Endspot[j,i]=np.nan
            else:
                Endspot[j,i]=np.nan
    inputDistribution[np.isnan(Endspot)==True]==np.nan
    return inputDistribution,Endspot
#%% Atlantic Data
latDA=np.linspace(-80,80,160)
#latDA=np.linspace(-90,90,180)
lonDA=np.linspace(-180,180,360)

location='D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\OutputFiles\Onink et al\Densities/'
File=['NorthAtlanticTotalDensity24h','NorthAtlanticEkmanDensity','NorthAtlanticGeostrophicDensity',
      'NorthAtlanticStokesDensity','NorthAtlanticStokesTotal24h']
zonalMeanA,meriMeanA=[],[]
for i in range(len(File)):
    density=np.load(location+File[i])
    density[np.isnan(density)]=0
    meanFinalYearA=np.sum(density[-183:,:,:]/density[-183:,:,:].shape[0],axis=0)
    #Now the zonal and meridional means between 0-50 N and 90-30W
    zonalMeanA.append(np.mean(meanFinalYearA,axis=0))
    meriMeanA.append(np.mean(meanFinalYearA,axis=1))
    meriMeanA[i]=meriMeanA[i]/integrate.simps(meriMeanA[i],lonDA)
    meriMeanA[i]=meriMeanA[i]*np.mean(meriMeanA[0])/np.mean(meriMeanA[i])
    zonalMeanA[i]=zonalMeanA[i]/integrate.simps(zonalMeanA[i],latDA)
    zonalMeanA[i]=zonalMeanA[i]*np.mean(zonalMeanA[0])/np.mean(zonalMeanA[i])
e=io.loadmat('D:\Desktop\Thesis\Data sets\Standardized Data Sets/SamplesMaxiGrid.mat')
sampleA=e['MaxiSamples']
sampleA=np.hstack((sampleA[:,180:],sampleA[:,:180]))
spotA=BasinShaperAtlantic(sampleA,lonDA,latDA)
sampleA[np.isnan(spotA)==True]=np.nan
sampleA[np.isnan(sampleA)==True]=0
zonalMeanA.append(np.mean(sampleA[79:171,89:240],axis=1))#79:151,89:151
meriMeanA.append(np.mean(sampleA[79:171,89:240],axis=0))
meriMeanA[-1]=meriMeanA[-1]/integrate.simps(meriMeanA[-1],lonDA[89:240])
zonalMeanA[-1]=zonalMeanA[-1]/integrate.simps(zonalMeanA[-1],latDA[79:171])
#%% Pacific Data
location='D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\OutputFiles\Onink et al\Densities/'
File=['NorthPacificTotalDensity24h','NorthPacificEkmanDensity','NorthPacificGeostrophicDensity',
      'NorthPacificStokesDensity','NorthPacificStokesTotalDensity24h']
latDP=np.linspace(-80,80,160)
lonDP=np.linspace(0,359,360)
zonalMeanP,meriMeanP=[],[]
for i in range(len(File)):
    density=np.load(location+File[i])
    density[np.isnan(density)]=0
    meanFinalYear=np.sum(density[-183:,:,:]/density[-183:,:,:].shape[0],axis=0)
    #Means between 0-60N latitude and east of 120 E
    zonalMeanP.append(np.mean(meanFinalYear,axis=0))
    meriMeanP.append(np.mean(meanFinalYear,axis=1))
    meriMeanP[i]=meriMeanP[i]/integrate.simps(meriMeanP[i],lonDP)
    zonalMeanP[i]=zonalMeanP[i]/integrate.simps(zonalMeanP[i],latDP)
#Observations
e=io.loadmat('D:\Desktop\Thesis\Data sets\Standardized Data Sets/SamplesMaxiGrid.mat')
sampleP=e['MaxiSamples']
lonSamP=e['lonSam']
lonSamP[lonSamP<0]+=360
latSamP=e['latSam']
LonSamP,LatSamP=np.meshgrid(lonSamP,latSamP)
_,spotP=BasinShaperPacific(sampleP,lonSamP,latSamP)
sampleP[np.isnan(spotP)==True]=np.nan
sampleP[np.isnan(sampleP)==True]=0
zonalMeanP.append(np.mean(sampleP[79:151,120:280],axis=1))
meriMeanP.append(np.mean(sampleP[79:151,120:280],axis=0))
meriMeanP[-1]=meriMeanP[-1]/integrate.simps(meriMeanP[-1],lonDP[120:280])
zonalMeanP[-1]=zonalMeanP[-1]/integrate.simps(zonalMeanP[-1],latDP[79:151])

#%%
labels=['Total','Ekman','Geostrophic','Stokes','Total + Stokes','Observations']
Colors=['red','forestgreen','dodgerblue','gold','darkmagenta','black']
axeslabelsize=14
textsize=12
fig=plt.figure(figsize=(10*1.5,8*1.5))
ax1 = fig.add_subplot(223)
for i in range(len(zonalMeanA)-1):
    ax1.plot(latDA,zonalMeanA[i],label=labels[i],color=Colors[i],linewidth=2)
#ax1.set_xlabel('Latitude ($^{\circ}$N)',fontsize=axeslabelsize)
ax1.set_ylabel('Normalized Density (# km$^{-1}$)',fontsize=axeslabelsize)
ax1.tick_params(labelsize=textsize)
ym1=0.3
ax1.set_xlim([0,50])
ax1.set_ylim([0,ym1])
ax1.set_xticks([0,10,20,30,40,50,60,70,80])
ax1.set_yticks(np.linspace(0,ym1,6))
ax1.set_title('(c) North Atlantic Zonal Mean',fontsize=axeslabelsize,fontweight='bold')
#ax2= ax1.twinx()
ax1.plot(latDA[79:171],zonalMeanA[-1],label=labels[-1],color=Colors[-1],linewidth=2)
#ax2.set_xlabel('Latitude ($^{\circ}$N)',fontsize=axeslabelsize)
#ax2.set_ylabel('Sampled Density\n'+r'($\times 10^{5}$ # km$^{-2}$)',fontsize=axeslabelsize)
#ax2.tick_params(labelsize=textsize)
#ax2.set_xlim([0,50])
#ym2=0.25
#ax2.set_ylim([0,ym1])
#ax2.set_xticks([0,10,20,30,40,50,60])
#ax2.set_yticks(np.linspace(0,ym2,6))
h1, l1 = ax1.get_legend_handles_labels()
#h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1, l1, fontsize=axeslabelsize)

ax3 = fig.add_subplot(224)
for j in range(len(meriMeanA)-1):
    ax3.plot(lonDA,meriMeanA[j],label=labels[j],color=Colors[j],linewidth=2)
ax3.set_xlim([-90,-30])
ym3=0.15
ax3.set_ylim([0,ym3])
ax3.set_yticks(np.linspace(0,ym3,6))
ax3.set_xticks([-90,-80,-70,-60,-50,-40,-30,-20,-10])
ax3.set_xticklabels([-90,-80,-70,-60,-50,-40,-30,-20,-10])
#ax3.set_xlabel('Longitude ($^{\circ}$W)',fontsize=axeslabelsize)
ax3.tick_params(labelsize=textsize)
#ax3.set_ylabel('Modelled Density\n'+r'($\times 10^{-3}$ # km$^{-2}$)',fontsize=axeslabelsize)
ax3.set_title('(d) North Atlantic Meridional Mean',fontsize=axeslabelsize,fontweight='bold')
#ax4 = ax3.twinx()
ax3.plot(lonDA[89:240],meriMeanA[-1],label=labels[-1],color=Colors[-1],linewidth=2)
#ax4.set_xlim([-90,-30])
#ax4.set_ylim([0,2])
#ax4.set_yticks(np.linspace(0,2,6))
#ax4.set_xticks([-90,-80,-70,-60,-50,-40,-30])
#ax4.set_xticklabels([-90,-80,-70,-60,-50,-40,-30])
#ax4.set_xlabel('Longitude ($^{\circ}$E)',fontsize=axeslabelsize)
#ax4.tick_params(labelsize=textsize)
#ax4.set_ylabel('Observed Density\n'+r'($\times 10^{5}$ # km$^{-2}$)',fontsize=axeslabelsize)

################################################################################
# Here we now add all the pacific plot stuff
################################################################################
ax5 = fig.add_subplot(221)
for i in range(len(zonalMeanP)-1):
    ax5.plot(latDP,zonalMeanP[i],label=labels[i],color=Colors[i],linewidth=2)
#plt.title('Zonal Mean')
ax1.set_xlabel('Latitude ($^{\circ}$N)',fontsize=axeslabelsize)
ax5.set_ylabel('Normalized Density (# km$^{-1}$)',fontsize=axeslabelsize)
ax5.tick_params(labelsize=textsize)
ym5=0.2
ax5.set_xlim([0,60])
ax5.set_ylim([0,ym5])
ax5.set_xticks([0,10,20,30,40,50,60,70,80])
ax5.set_yticks(np.linspace(0,ym5,6))
ax5.set_title('(a) North Pacific Zonal Mean',fontsize=axeslabelsize,fontweight='bold')
#ax7= ax5.twinx()
ax5.plot(latSamP[0,79:151],zonalMeanP[-1],label=labels[-1],color=Colors[-1],linewidth=2)
#ax7.set_xlabel('Latitude ($^{\circ}$N)',fontsize=axeslabelsize)
#ax7.set_ylabel('Sampled Density\n'+r'($\times 10^{5}$ # km$^{-2}$)',fontsize=axeslabelsize)
#ax7.tick_params(labelsize=textsize)
#ax7.set_xlim([0,60])
#ax7.set_ylim([0,2])
#ax7.set_xticks([0,10,20,30,40,50,60])
#ax7.set_yticks(np.linspace(0,2,6))
#h5, l5 = ax5.get_legend_handles_labels()
#h7, l7 = ax7.get_legend_handles_labels()
#ax5.legend(h5+h7, l5+l7, fontsize=16)

ax6 = fig.add_subplot(222)
for j in range(len(meriMeanP)-1):
    ax6.plot(lonDP,meriMeanP[j],label=labels[j],color=Colors[j],linewidth=2)
ax6.set_xlim([120,280])
ax6.set_xticks([120,140,160,180,200,220,240,260,280])
ax6.set_xticklabels([120,140,160,180,-160,-140,-120,-100,-80])
ym6=0.15
ax6.set_ylim([0,ym6])
ax6.set_yticks(np.linspace(0,ym6,6))
ax3.set_xlabel('Longitude ($^{\circ}$)',fontsize=axeslabelsize)
ax6.tick_params(labelsize=textsize)
ax6.set_title('(b) North Pacific Meridional Mean',fontsize=axeslabelsize,fontweight='bold')
#ax6.set_ylabel('Modelled Density\n'+r'($\times 10^{-3}$ # km$^{-2}$)',fontsize=axeslabelsize)
#ax8 = ax6.twinx()
ax6.plot(lonSamP[0,120:280],meriMeanP[-1],label=labels[-1],color=Colors[-1],linewidth=2)
#ax8.set_xlim([120,280])
#ax8.set_ylim([0,2.5])
#ax8.set_yticks(np.linspace(0,2.5,6))
#ax8.set_xticks([120,140,160,180,200,220,240,260,280])
#ax8.set_xticklabels([120,140,160,180,-160,-140,-120,-100,-80])
#ax8.set_xlabel('Longitude ($^{\circ}$)',fontsize=axeslabelsize)
#ax8.tick_params(labelsize=textsize)
#ax8.set_ylabel('Observed Density\n'+r'($\times 10^{5}$ # km$^{-2}$)',fontsize=axeslabelsize)
plt.tight_layout()
plt.savefig('D:\Desktop\Thesis\ParcelsFigData\Data/AtlanticPacificComponentComparisonNormalFull.jpg')
