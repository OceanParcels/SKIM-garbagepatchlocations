# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:43:49 2018

@author: Victor Onink
"""

from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4,ErrorCode, plotTrajectoriesFile,Variable
from datetime import timedelta, datetime
import numpy as np
from operator import attrgetter
import math
#We can add or remove all the zeros according to preference. In case that they are left there, we only get daily data for the currents which will end up with the code running faster, but we do lose time resolution. Tests will determine if this loss in time resolution is actually important
filenames = {'U': "/scratch/Victor/GeostrophicData/20*000000*.nc",
             'V': "/scratch/Victor/GeostrophicData/20*000000*.nc",
             #'avgU':"/scratch/Victor/AvgTotCur/AverageG*.nc",
             #'avgV':"/scratch/Victor/AvgTotCur/AverageG*.nc",
	     'borU':"/scratch/Victor/AvgTotCur/boundary_velocitiesG*",
	     'borV':"/scratch/Victor/AvgTotCur/boundary_velocitiesG*"}
variables = {'U': 'eastward_geostrophic_current_velocity',
             'V': 'northward_geostrophic_current_velocity',
             #'avgU':'Uavg',
             #'avgV':'Vavg',
	     'borU':'MaskUvel',
	     'borV':'MaskVvel'}
dimensions = {'lat': 'lat',
              'lon': 'lon',
              'time': 'time'}
#%%
#Create the fieldset with the periodic halo and time extrapolation for the EKE
print 'Creating the fieldset'
fieldset = FieldSet.from_netcdf(filenames, variables, dimensions,allow_time_extrapolation=True)
fieldset.add_periodic_halo(zonal=True)

#The starting coordinates of the Particles, for the North Pacific. They are generated
#by the code NAgrid.py, graciously send to me by David.
lons=np.load('/home/students/4056094/Desktop/Thesis/ParcelsOutput/North Atlantic/InputDistribution/LonsTestgrid0_5.npy')
lats=np.load('/home/students/4056094/Desktop/Thesis/ParcelsOutput/North Atlantic/InputDistribution/LatsTestgrid0_5.npy')
#lons, lats = np.meshgrid(lon,lat)


#And now we define what sort of particles we are actually dealing with
class SampleParticle(JITParticle):
   # EKE=Variable('EKE', initial=0., dtype=np.float32)
#    #Now the part to determine the age of the particle
    Age=Variable('Age',initial=0.,dtype=np.float32)#agr is gonna be in seconds
    prev_time=Variable('prev_time',initial=attrgetter('time'),to_write=False)
    #Now the part to track the distance covered
#    distance = Variable('distance', initial=0., dtype=np.float32)
#    prev_lon = Variable('prev_lon', dtype=np.float32, to_write=False,
#                        initial=attrgetter('lon'))
#    prev_lat = Variable('prev_lat', dtype=np.float32, to_write=False,
#                        initial=attrgetter('lat'))
#    #Now I also want the particle to be deleted if it is on land (so it won't move)
#    count=Variable('count',initial=0,to_write=False)
#    init_lon = Variable('init_lon', dtype=np.float32, to_write=False,
#                        initial=attrgetter('lon'))
#    init_lat = Variable('init_lat', dtype=np.float32, to_write=False,
#                        initial=attrgetter('lat'))
#The starting point of the similation and the endtime
print 'Creating the pset'
starttime=datetime(2002,1,1,0,0)
endtime=datetime(2014,12,31,21,0)
pset = ParticleSet(fieldset=fieldset, pclass=SampleParticle, lon=lons, lat=lats,time=starttime)
#%% All the different functions/kernels we want to have
def moveOver(particle,fieldset,time,dt):
    if particle.lon<0:
	particle.lon+=360.0
    elif particle.lon>360:
	particle.lon-=360.0
def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()
    print 'we deleted it at '+str(particle.lon)+' and '+str(particle.lat)
#def EKEsample(particle, fieldset, time, dt):
#    u=fieldset.U[time, particle.lon, particle.lat, particle.depth]
#    v=fieldset.V[time, particle.lon, particle.lat, particle.depth]
#    u_avg=fieldset.avgU[time, particle.lon, particle.lat, particle.depth]
#    v_avg=fieldset.avgV[time, particle.lon, particle.lat, particle.depth]
#    u_p=u-u_avg
#    v_p=v-v_avg
#    particle.EKE=(u_p*u_p+v_p*v_p)/2
def AgeSample(particle, fiedset,time,dt):
    current_time=particle.time
    timedifference=current_time-particle.prev_time
    particle.Age+=timedifference
    particle.prev_time=current_time
#def TotalDistance(particle, fieldset, time, dt):
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
#    lat_dist = (particle.lat - particle.prev_lat) * 1.11e2
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
#    lon_dist = (particle.lon - particle.prev_lon) * 1.11e2 * math.cos(particle.lat * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
#    particle.distance += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))
#    particle.prev_lon = particle.lon  # Set the stored values for next iteration.
#    particle.prev_lat = particle.lat
def antiBeach(particle,fieldset,time,dt):
    bu=fieldset.borU[time,particle.lon,particle.lat,particle.depth]
    bv=fieldset.borV[time,particle.lon,particle.lat,particle.depth]
    particle.lon-=bu*dt*0.00001
    particle.lat-=bv*dt*0.00001
#EKEsam=pset.Kernel(EKEsample)
Agesam=pset.Kernel(AgeSample)
#Distsam=pset.Kernel(TotalDistance) 
Beach=pset.Kernel(antiBeach)
#Landsam=pset.Kernel(antiInitialLand)
move=pset.Kernel(moveOver)
totalKernal=move+AdvectionRK4+Beach+Agesam
#%%
pfile = pset.ParticleFile(name="/scratch/Victor/AtlanticGeostrophic",
                          outputdt=timedelta(hours=48))

Time=starttime
steps=0
while Time<=endtime:
    steps+=1
    Time+=timedelta(hours=48)
print 'now we start advecting them for how many steps? '+str(steps)
for i in range(steps-1):
#    if (i+1)%steps==0:
#        pset.show(field='vector',particles=True,land=True,domain=[80,0,20,-120])
#    if i%20==0:
#	print i
    pset.execute(totalKernal,
                 runtime=timedelta(hours=48),  # runtime controls the interval of the plots
                 dt=timedelta(minutes=30),
                 recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
                 output_file=pfile
                 )  # the recovery kernel
#%%

