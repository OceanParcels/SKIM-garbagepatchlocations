# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 09:47:10 2018

@author: Victor Onink
Downloading daily mean total currents
"""
import urllib
import numpy as np
import datetime

#testfile = urllib.URLopener()
#testfile.retrieve(url, "D:\Desktop\Thesis/works.nc")

direc='http://www.ifremer.fr/opendap/cerdap1/globcurrent/v3.0/global_025_deg/total_hs' 
savedirec='/scratch/Victor/TotalData24h'
s='/' #this marks the directories, not \ since that fails with numbers in the directory names
standard='-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v03.0-fv01.0.nc'
for i in range(2002,2015):
        direcyear=direc+s+str(i)
        if i==2004 or i==2008 or i==2012:
            days=366
        else:
            days=365
	if i==2006:
	    start=1
	else:
	    start=1 
        for j in range(start,days+1):
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
	    File=direcday+datestamp+standard
		#if i==2013 and j==160 and k==6:
		 #   File='http://www.ifremer.fr/opendap/cerdap1/globcurrent/v2.0/global_025_deg/total_hs/20130609180000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v02.0-fv01.0.nc'
            SaveName=savedirec+s+datestamp+standard
	#	print time[k]
            testfile = urllib.URLopener()
            testfile.retrieve(File, SaveName)
