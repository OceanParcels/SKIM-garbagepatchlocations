# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 09:47:10 2018

@author: Victor Onink
Downloading the globcurrent geostrophic currents
"""
import urllib
import numpy as np
import datetime

#testfile = urllib.URLopener()
#testfile.retrieve(url, "D:\Desktop\Thesis/works.nc")

direc='http://www.ifremer.fr/opendap/cerdap1/globcurrent/v3.0/global_025_deg/geostrophic' 
savedirec='/scratch/Victor/GeostrophicData'
s='/' #this marks the directories, not \ since that fails with numbers in the directory names
standard='-GLOBCURRENT-L4-CURgeo_0m-ALT_OI-v03.0-fv01.0.nc'
for i in range(2002,2015):
        direcyear=direc+s+str(i)
        if i==2004 or i==2008 or i==2012:
            days=366
        else:
            days=365 
	if i==2011:
	    start=1
	else:
	    start=1
        for j in range(start,start+days):
	    if j%5==0:
		print j
            if j<10:
                day='00'+str(j)
            elif j>=10 and j<100:
                day='0'+str(j)
            else:
                day=str(j)
            direcday=direcyear+s+day+s
            Date=datetime.datetime(i, 1, 1) + datetime.timedelta(j - 1)
            [y, m, d]=[Date.year,Date.month,Date.day]
            if m<10:
                if d<10:
                    datestamp=str(y)+'0'+str(m)+'0'+str(d)
                else:
                    datestamp=str(y)+'0'+str(m)+str(d)
            else:
                if d<10:
                    datestamp=str(y)+str(m)+'0'+str(d)
                else:
                    datestamp=str(y)+str(m)+str(d)
            time=['000000']
            for k in range(len(time)):
                File=direcday+datestamp+time[k]+standard
                SaveName=savedirec+s+datestamp+time[k]+standard
                testfile = urllib.URLopener()
                testfile.retrieve(File, SaveName)

