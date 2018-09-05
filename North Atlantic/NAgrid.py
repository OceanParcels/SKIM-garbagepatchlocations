#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 17:21:28 2018

Function to generate a uniform (in degree) grid of particles for initial positions
for the North Atlantic Basin
@author: David Wichmann, modified by Victor Onink
"""
import numpy as np
from parcels import Field
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def GetGlobCurrentLandArray(filename, fieldname):
    """
    Function to return a Field with 1's at land points and 0's at ocean points, based on python basemap
    :param f: a field of the .nc file, but not a parcels field object! For OFES, land points are masked. This is used here!
    """
    pfile = Dataset(filename, 'r')
    Lon = pfile.variables['lon'][:]
    Lat = pfile.variables['lat'][:]
    f = pfile.variables[fieldname][:]
    f = f[0]
    L= np.ma.getmask(f)
    Land=Field('Land',L,transpose=False,lon=Lon,lat=Lat)
    return Land

def Plotdistr():
#    lons=np.load('D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\InputArray/NA_TM_LonsInitial_ddeg02.npy')
#    lats=np.load('D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\InputArray/NA_TM_LatsInitial_ddeg02.npy')
#    lons= np.load('/Users/wichmann/sources/inertialparticles_David/paper/1_Distribution/ParticleGrid/Lons_N1275_NA_deg2.npy')
#    lats = np.load('/Users/wichmann/sources/inertialparticles_David/paper/1_Distribution/ParticleGrid/Lats_N1275_NA_deg2.npy')   
#    slons= np.load('/Users/wichmann/sources/inertialparticles_David/paper/1_Distribution/ParticleGrid/Lons_N1275_NA_deg2_shift.npy')
#    slats = np.load('/Users/wichmann/sources/inertialparticles_David/paper/1_Distribution/ParticleGrid/Lats_N1275_NA_deg2_shift.npy')   
    
    print lons[lons<0.]
#    print slons[slons<0.]
#    
#    m = Basemap(projection='merc', llcrnrlat=0.,  urcrnrlat=75., llcrnrlon=-60.0, urcrnrlon=360., resolution='l')
#    m.drawcoastlines()
#    xs, ys = m(lons, lats)
#    m.scatter(xs,ys)
    
#Plotdistr()

def GenerateP(landfilename, lonmin, lonmax, latmin, latmax, spacing, name):
    Land = GetGlobCurrentLandArray(landfilename, 'northward_eulerian_current_velocity')
    grid=np.mgrid[lonmin:lonmax:spacing,latmin:latmax:spacing]
    n=grid[0].size;
    lons=np.reshape(grid[0],n)
    lats=np.reshape(grid[1],n)

    print 'number of particles beginning: ', len(lats)
    #Remove land particles and mdeiterranean particles
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,lo,la,0]==0.0  ]),np.array([la for lo, la in zip(lons,lats) if Land[0,lo,la,0]==0.0 ])]
    
#    #Remove pacific part
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la< - lo - 75)  ]),np.array([la for lo, la in zip(lons,lats) if not (la< - lo - 75)])]
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la< 10 and lo<-70)  ]),np.array([la for lo, la in zip(lons,lats) if not (la< 10 and lo<-70)])]    
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la< 15 and lo<-83)  ]),np.array([la for lo, la in zip(lons,lats) if not (la< 15 and lo<-83)])]    
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la< 12 and lo<-82)  ]),np.array([la for lo, la in zip(lons,lats) if not (la< 12 and lo<-82)])]    
#
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la> 50 and -70<lo<-20)  ]),np.array([la for lo, la in zip(lons,lats) if not (la> 50 and -70<lo<-20)])]    

#    #Remove mediterranian
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la< 48 and la > 30 and lo >0)  ]),np.array([la for lo, la in zip(lons,lats) if not (la< 48 and la > 30 and lo >0)])]    
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la< 39 and la > 35 and lo >-7)  ]),np.array([la for lo, la in zip(lons,lats) if not (la< 39 and la > 35 and lo >-7)])]        
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (lo<5. and lo>-5 and la<41 and la>35)  ]),np.array([la for lo, la in zip(lons,lats) if not (lo<5. and lo>-5 and la<41 and la>35)])]        
#
#    #remove some parts in canada
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la> 50 and la<72 and lo < -70)  ]),np.array([la for lo, la in zip(lons,lats) if not (la> 50  and la<72 and lo < -70)])]            
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (lo < -78 and la > 71)  ]),np.array([la for lo, la in zip(lons,lats) if not (lo < -78 and la > 71)])]            
#
#    #Remove Ostsee
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la<65 and lo > 10)  ]),np.array([la for lo, la in zip(lons,lats) if not (la<65 and lo > 10)])]            
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la>60 and la<67 and lo > 20)  ]),np.array([la for lo, la in zip(lons,lats) if not (la>60 and la<67 and lo > 20)])]            
#    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,lo,la,0]==0.0 and not (lo>354. and la<40. and la>30.) and not (lo<288. and la>50.)  and not (lo<295.5 and lo>290. and la>57.) and not (lo<283. and la<8.6) and not (lo>357. and la> 52.)]),np.array([la for lo, la in zip(lons,lats) if Land[0,lo,la,0]==0.0 and not (lo>354. and la<40. and la>30.) and not (lo<288. and la>50.) and not (lo<295.5 and lo>290. and la>57.) and not (lo<283. and la<8.6) and not (lo>357. and la>52.)])]    

    print 'number of particles end: ', len(lats)
    
    #Plot particles and density to check if it is correct
    fig = plt.figure(figsize=(22, 16))
    ax = fig.add_subplot(211)
    ax.set_title("Particles")
    m = Basemap(projection='merc', llcrnrlat=latmin,  urcrnrlat=latmax, llcrnrlon=lonmin, urcrnrlon=lonmax, resolution='l')
    m.drawparallels(np.array([20,40,60]), labels=[True, False, False, True])
    m.drawmeridians(np.array([280,300,320,340,360,380]), labels=[False, False, False, True])
    m.drawcoastlines()
    xs, ys = m(lons, lats)
    m.scatter(xs,ys)
    
    ax = fig.add_subplot(212)
    ax.set_title("Particles per bin. Should be 1 everywhere but on land.")
    m = Basemap(projection='merc', llcrnrlat=latmin,  urcrnrlat=latmax, llcrnrlon=lonmin, urcrnrlon=lonmax, resolution='l')
    m.drawcoastlines()
    lon_bins = np.arange(lonmin, lonmax, spacing)
    lat_bins = np.arange(latmin, latmax, spacing)
    density, _, _ = np.histogram2d(lats, lons, [lat_bins, lon_bins])
    lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)
    xs, ys = m(lon_bins_2d, lat_bins_2d)
    plt.pcolormesh(xs, ys, density,cmap=plt.cm.RdBu_r)
    cbar = plt.colorbar(orientation='vertical', shrink=0.625, aspect=20, fraction=0.2,pad=0.02)
    cbar.set_label('Particles per bin',size=8)
        
#    np.save('D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\InputArray/Lons' + name,lons)
#    np.save('D:\Desktop\Thesis\ParcelsFigData\Data\North Atlantic\InputArray/Lats' + name,lats)
    
    plt.show()
    


landfilename = "D:\Desktop\Thesis\Data sets\GlobCurrent/2005/AllTogether/20050101000000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc"
lonmin, lonmax = -100., 30.
latmin, latmax = 0., 75.
spacing = 0.5
name = 'GeoTestgrid0_5'
GenerateP(landfilename, lonmin, lonmax, latmin, latmax, spacing, name)

