# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 12:17:27 2018

@author: Erik van Sebille, modified by Victor Onink
Code is from Erik and it is used to create a current at the land boundaries which pushes
particles away from the coast, preventing the high amount of beaching that we observe
Now the goal is to adapt this to globcurrent...
"""
from parcels import FieldSet
import numpy as np
import xarray
from progressbar import ProgressBar


def set_globcurrent_fieldset():
    datadir = 'D:\Desktop\Thesis\Data sets/20141231000000-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v03.0-fv01.0.nc'
    filenames = {'U': datadir,
                 'V': datadir}
    variables = {'U': 'eastward_geostrophic_current_velocity', 'V': 'northward_geostrophic_current_velocity'}
    dimensions = {'lat': 'lat', 'lon': 'lon',
                  'time': 'time'}
    return FieldSet.from_netcdf(filenames, variables, dimensions)

#def set_globcurrent_fieldset():
#    datadir = 'D:\Desktop\Thesis\Data sets/20141231-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v03.0-fv01.0.nc*'
#    filenames = {'U': datadir,
#                 'V': datadir}
#    variables = {'U': 'eastward_eulerian_current_velocity', 'V': 'northward_eulerian_current_velocity'}
#    dimensions = {'lat': 'lat', 'lon': 'lon',
#                  'time': 'time'}
#    return FieldSet.from_netcdf(filenames, variables, dimensions)

def isocean(u, v):
    return True if u == 0 and v == 0 else False


def landborders(u, v, J, I, nx):
    mask = np.ones((3, 3), dtype=bool)
    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            mask[j+1, i+1] = isocean(u[ j+J, (i+I) % nx], v[ j+J, (i+I) % nx])
    return mask


fset = set_globcurrent_fieldset()
nx = fset.U.lon.size
ny = fset.U.lat.size
#nz = fset.U.depth.size

mask_uvel = np.zeros_like(fset.U.data)
mask_vvel = np.zeros_like(fset.U.data)
#for k in range(0, nz):
#    print k
pbar = ProgressBar()
for i in pbar(range(0, nx)):
    for j in range(1, ny-1):
        if isocean(fset.U.data[0, j, i], fset.V.data[0, j, i]):
            mask = landborders(fset.U.data[0, :, :], fset.V.data[0, :, :], j, i, nx)
            if not mask.all():
                mask_uvel[0, j, i] = sum(mask[:, 2]) - sum(mask[:, 0])
                mask_vvel[0, j, i] = sum(mask[2, :]) - sum(mask[0, :])


direc = 'D:\Desktop\Thesis\ParcelsFigData\Data\BorderCurrents/'
coords = [('time', [0]),
          ('lat', fset.U.lat), ('lon', fset.U.lon)]
uvel = xarray.DataArray(mask_uvel, coords=coords)
vvel = xarray.DataArray(mask_vvel, coords=coords)
dcoo = {'lon':  fset.U.lon, 'lat': fset.U.lat,
        'time': [0]}
dset = xarray.Dataset({'MaskUvel': uvel, 'MaskVvel': vvel}, coords=dcoo)
dset.to_netcdf(direc+"boundary_velocities_v3Geo.nc")