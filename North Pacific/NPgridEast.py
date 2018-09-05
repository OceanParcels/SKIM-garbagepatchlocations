#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 17:21:28 2018

Function to generate a uniform (in degree) grid of particles for initial positions

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
    #Remove land particles 
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,lo,la,0]==0.0  ]),np.array([la for lo, la in zip(lons,lats) if Land[0,lo,la,0]==0.0 ])]
    
    #Remove the carribean
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la> 18 and lo>-100)  ]),np.array([la for lo, la in zip(lons,lats) if not (la> 18 and lo>-100)])]    
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la> 14 and lo>-90)  ]),np.array([la for lo, la in zip(lons,lats) if not (la> 14 and lo>-90)])]    
    [lons,lats]=[np.array([lo for lo, la in zip(lons,lats) if not (la> 9 and lo>-85)  ]),np.array([la for lo, la in zip(lons,lats) if not (la> 9 and lo>-85)])]    

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
        
    np.save('D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\InputArray/Lons' + name,lons)
    np.save('D:\Desktop\Thesis\ParcelsFigData\Data\North Pacific\InputArray/Lats' + name,lats)
    
    plt.show()
    


landfilename = "D:\Desktop\Thesis\Data sets\GlobCurrent/2005/AllTogether/20050101000000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc"
lonmin, lonmax = -155, -125.
latmin, latmax = 25, 35.
spacing = 0.1
name = 'FreqFourier'
GenerateP(landfilename, lonmin, lonmax, latmin, latmax, spacing, name)

