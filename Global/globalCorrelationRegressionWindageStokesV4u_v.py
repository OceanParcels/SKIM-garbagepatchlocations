# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 13:27:01 2018

@author: Victor Onink
Going to see if we can create a correlation/regression map for between Stokes drift
and the various windage velocity fields
We also consider RMSE between Stokes and windage scenarios 
We will consider 1, 3 and 5 % windages
"""
import numpy as np
from netCDF4 import Dataset
import datetime
from scipy import stats

def rmse(predictions, targets):
    return np.sqrt(np.nanmean(((predictions - targets) ** 2)))

def rmseAngle(predictions, targets):
    return np.sqrt(np.nanmean((((predictions - targets)%360) ** 2)))

def netCDFcreator(filename,corField,lon,lat,label):
    root_grp = Dataset(filename, 'w', format='NETCDF4')
    root_grp.description = label
    #Create Dimensions
    latSize=lat.size
    lonSize=lon.size
    root_grp.createDimension('windages', corField.shape[0])
    root_grp.createDimension('latitude', latSize)
    root_grp.createDimension('longitude', lonSize)
    # variables
    levels = root_grp.createVariable('windages', 'f8', ('windages',))
    latitude = root_grp.createVariable('latitude', 'f8', ('latitude',))
    longitude = root_grp.createVariable('longitude', 'f8', ('longitude',))
    Field = root_grp.createVariable('Field', 'f8', ('windages', 'latitude', 'longitude',))
    #Now assigning values to the variables
    latitude[:]=lat
    longitude[:]=lon
    if corField.shape[0]==3:
	levels[:]=[0.01,0.03,0.05,20]
    if corField.shape[0]==2:
	levels[:]=[0,0]
    Field[:,:,:]=corField
    #Finally, close the new field
    root_grp.close()
#%%

filenameWind='/scratch/Victor/WaveWatchWind24h/WaveWatchWindDailyMean'
filenameStokes='/scratch/Victor/StokeData24h/StokesDriftDailyMean'

#Fields within which we save the correlation and regression values
CorrelationFieldMag=np.ones((2,317,720))
RMSEFieldmag=np.ones((4,317,720))
CorrelationFieldAngle=np.ones((1,317,720))
RMSEFieldAngle=np.ones((1,317,720))

#for i in range(5):#range(317):
for j in range(720):
    if j%50==0:
	print j
    stokesMagFull=np.zeros((317,366)) 
    stokesMagFull[:]=np.nan
    stokesU=np.zeros((317,366))
    stokesV=np.zeros((317,366))
    stokesU[:],stokesV[:]=np.nan,np.nan
    #stokesAngleFull=np.zeros((317,366))
    #stokesAngleFull[:]=np.nan
    windMagFull=np.zeros((317,366))
    windMagFull[:]=np.nan
    windU=np.zeros((317,366))
    windV=np.zeros((317,366))
    windU[:],windV[:]=np.nan,np.nan
    #windAngleFull=np.zeros((317,366))
    #windAngleFull[:]=np.nan
    #first we need to determine the datestamp for the files
    for k in range(2002,2015):
        if k==2004 or k==2008 or k==2012:
            days=366
        else:
            days=365 
        for d in range(1,days+1):
            #if d%100==0:
            #    print d
            if d<10:
                day='00'+str(d)
            elif d>=10 and d<100:
                day='0'+str(d)
            else:
                day=str(d)
            Date=datetime.datetime(k, 1, 1) + datetime.timedelta(d - 1)
            [year, mon, day]=[Date.year,Date.month,Date.day]
            if mon<10:
                if day<10:
                    datestamp=str(year)+'0'+str(mon)+'0'+str(day)
                else:
                    datestamp=str(year)+'0'+str(mon)+str(day)
            else:
                if day<10:
                    datestamp=str(year)+str(mon)+'0'+str(day)
                else:
                    datestamp=str(year)+str(mon)+str(day)
            windFile=filenameWind+datestamp+'.nc'
            dataWind=Dataset(windFile)
            windMag=np.sqrt(np.square(dataWind.variables['uwnd'][0,:,j]) \
                            +np.square(dataWind.variables['vwnd'][0,:,j]))
            windMagFull[:,d-1]=windMag
            windU[:,d-1]=dataWind.variables['uwnd'][0,:,j]
	    windV[:,d-1]=dataWind.variables['vwnd'][0,:,j]
            #windAngleFull[:,d-1]=np.arctan2(dataWind.variables['uwnd'][0,:,j], \
            #                                dataWind.variables['vwnd'][0,:,j])/np.pi*180
            
            stokesFile=filenameStokes+datestamp+'.nc'
            dataStokes=Dataset(stokesFile)
            stokesMag=np.sqrt(np.square(dataStokes.variables['uuss'][0,:,j])\
                              +np.square(dataStokes.variables['vuss'][0,:,j]))
            stokesMagFull[:,d-1]=stokesMag
	    stokesU[:,d-1]=dataStokes.variables['uuss'][0,:,j]
	    stokesV[:,d-1]=dataStokes.variables['vuss'][0,:,j]
            #stokesAngleFull[:,d-1]=np.arctan2(dataStokes.variables['uuss'][0,:,j],\
            #                                  dataStokes.variables['vuss'][0,:,j])/np.pi*180
                           
    #Here we then compute the correlation coefficient
    for i in range(317):
        if np.isfinite(stokesMagFull[i,:]).any()==True:
            r_1 = np.ma.corrcoef(stokesU[i,:][np.isnan(stokesU[i,:])==False],windU[i,:][np.isnan(stokesU[i,:])==False])[0,1]
	    #r_1= np.ma.corrcoef(stokesU[i,:],windU[i,:])[0,1]
	    r_2 = np.ma.corrcoef(stokesV[i,:][np.isnan(stokesV[i,:])==False],windV[i,:][np.isnan(stokesV[i,:])==False])[0,1]
	    #r_2= np.ma.corrcoef(stokesV[i,:],windV[i,:])[0,1]
	    #print r_1, r_2
        else:
            r_1,r_2=np.nan,np.nan
	    #print 'we have nothing'
        CorrelationFieldMag[0,i,j]=r_1
	CorrelationFieldMag[1,i,j]=r_2
        #And now the root mean error between the time series
        if np.isfinite(stokesMagFull[i,:]).any()==True:
            RMSE_1=rmse(windMagFull[i,:]*0.01,stokesMagFull[i,:])
            RMSE_3=rmse(windMagFull[i,:]*0.03,stokesMagFull[i,:])
            RMSE_5=rmse(windMagFull[i,:]*0.05,stokesMagFull[i,:])
	    meanStokes=np.nanmean(stokesMagFull[i,:])
            errors=[RMSE_1,RMSE_3,RMSE_5,meanStokes]
        else:
            errors=[np.nan,np.nan,np.nan,np.nan]
        for z in range(len(errors)):
            RMSEFieldmag[z,i,j]=errors[z]
        #Now the correlation coefficient of the current direction
        #if np.isfinite(stokesAngleFull[i,:]).any()==True:
        #    r_1 = np.ma.corrcoef(stokesAngleFull[i,:][np.isnan(stokesAngleFull[i,:])==False],0.01*windAngleFull[i,:][np.isnan(windMagFull[i,:])==False])[0,1]
        #else:
        #    r_1=np.nan
        #CorrelationFieldAngle[0,i,j]=r_1
        #Now the RMSE of the angle between Stokes Drift and Windage
        #if np.isfinite(stokesAngleFull[i,:]).any()==True:
        #    error = rmseAngle(stokesAngleFull[i,:],windAngleFull[i,:])
        #else:
        #    error=np.nan
        #RMSEFieldAngle[0,i,j]=error
        
labelCor='Correlation of magnitude for time series at each 0.5x0.5 grid cel between Stokes drift and windages, 0=U, 1=V'
labelRMS='RMSE of magnitude for time series at each 0.5x0.5 grid cel between Stokes drift and windages'
#labelCorAngle='Correlation of angle for time series at each 0.5x0.5 grid cel between Stokes drift and windages'
#labelRMSangle='RMSE of angle for time series at each 0.5x0.5 grid cel between Stokes drift and windages'
lon=dataWind.variables['longitude'][:]
lat=dataWind.variables['latitude'][:]
filenameCor='/scratch/Victor/CorrelationsWindageStokesMagnitudeUV02_14.nc'
filenameRMSE='/scratch/Victor/RMSE_WindageStokesMagnitudeUV02_14.nc'
#filenameCorAngle='/scratch/Victor/CorrelationsWindageStokesAngle.nc'
#filenameRMSEangle='/scratch/Victor/RMSE_WindageStokesAngle.nc'
netCDFcreator(filenameCor,CorrelationFieldMag,lon,lat,labelCor)
netCDFcreator(filenameRMSE,RMSEFieldmag,lon,lat,labelRMS)
#netCDFcreator(filenameCorAngle,CorrelationFieldAngle,lon,lat,labelCorAngle)
#netCDFcreator(filenameRMSEangle,RMSEFieldAngle,lon,lat,labelRMSangle)


