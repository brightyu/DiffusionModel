# -*- coding: utf-8 -*-
"""
Created on Wed Dec 02 11:15:20 2015

@author: khan
"""

import csv
import os
from numpy import array
import numpy
import datetime
#from datetime import datetime
import time
import matplotlib.pyplot as plt
import dateutil
from pylab import *
import bisect
import random

rootdir= 'C:\\Users\\khan\\Documents\\GitHub\\Python_model_output\\'
filenamelocation = rootdir+'Florida hourly data from Mukul simplified.csv'
timestamp = time.time()
first_time_seconds = str(math.floor(timestamp))

data = recfromcsv(filenamelocation)
numbercolumns = 36
numberrows = len(data)

weatherarray = numpy.zeros((numberrows,numbercolumns))

for i in range(numbercolumns):
    for j in range(numberrows):
        weatherarray[j,i] = data[j][i]
# column 22 is temperature, 24 is humidity (C, %RH)        
           
TimeStep = 60 # minutes (keep this at one hour because weather data is every hour)
SimulationLength = 1 # years
TimePoints = SimulationLength * 365 * 24 * 60 / TimeStep

#Store T, RH data and calculate leakage current density based on leakage current model for Al with mitsui front gating
#from leakage current model
current_inv_tmep_param = -2531
current_hum_param = 0.009016
current_int_param = 6.964
weathercurrentarray = numpy.zeros((len(data),12))

weathercurrentarray[:,0] = weatherarray[:,22] # logs the termperature at each time point for 25 years (or simulation length)
weathercurrentarray[:,1] = weatherarray[i,24] # %RH for each tiem point
weathercurrentarray[:,2] = exp(current_int_param+\
current_inv_tmep_param/(weathercurrentarray[:,0]+273+20)\
+current_hum_param*weathercurrentarray[:,1]) #current in nanoamps
weathercurrentarray[:,3] = weatherarray[:,33] # precipitation in mm
weathercurrentarray[:,4] = weatherarray[:,23] # dew point in C
weathercurrentarray[:,6] = weatherarray[:,4] # illumination

XPET = 200e-6 # m
XEVA = 450e-6 # m

delt = .01 # seconds
xSteps = 650
tSteps = len(data)*100*3600
delx = (XPET+XEVA)/xSteps # m
persecond = 1/delt
print delt*tSteps/60


#DAPET = 1.71e-11 # m^2/s at 85C diffusion coefficient
#DAEVA = 6.33e-10 # m^2/s at 85C

k=8.617e-5 #eV/K
D0PET = 2.89e-5 #(m^2/s)
D0Mitsui = 7.57e-9    #m^2/sec
EaPET = 0.46 #eV
EaMitsui = 0.145 #eV
CAarray = numpy.zeros((xSteps))
#DAPETarray = numpy.zeros((tSteps)) # use this if we want to change temperature with time
#DAEVAarray = numpy.zeros((tSteps))

CAarray_past = numpy.zeros((xSteps))
CAarray_output = numpy.zeros((xSteps,11))

#tStepsarray = np.linspace(0,tSteps*delt/60.0,tSteps)
xStepsarray = np.linspace(0,xSteps*delx*1e6,xSteps)	

#DAPETarray[:] = D0PET*numpy.exp(-(EaPET)/k/(TRHarray[:,0]+273))
#DAEVAarray[:] = D0Mitsui*numpy.exp(-(EaMitsui)/k/(TRHarray[:,0]+273))

intPETstop = int(xSteps*XPET/(XPET+XEVA)) # location of PET-EVA boundary in xSteps

CAarray_past[:] = 1 # initial condition g/cm^3 water in encap  materials at beginning

#font = {'family' : 'sans-serif',
#        'weight' : 'bold',
#        'size'   : 22}
#title_font = {'fontname':'Arial', 'size':'12', 'color':'black', 'weight':'normal',
#          'verticalalignment':'bottom'} # Bottom vertical alignment for more space
#legend_font = {'fontname':'Arial', 'size':'18', 'color':'black', 'weight':'normal',
#          'verticalalignment':'bottom'} # Bottom vertical alignment for more space
#matplotlib.rc('font', **font)

#figure(num=None, figsize=(12, 8), dpi=480, facecolor='w', edgecolor='k')
#xlabel('Distance across domain (um)')
#ylabel('Water Concentration')
Temp = 10
Humidity = 50
DAPET = 0
DAEVA =  0
j=1

#for j in xrange(1,tSteps):
while (j<tSteps):
    Temp = weatherarray[int(j/3600.0/persecond),22]
    Humidity = weatherarray[int(j/3600.0/persecond),24]
    DAPET = D0PET*numpy.exp(-(EaPET)/k/(Temp+273))
    DAEVA = D0Mitsui*numpy.exp(-(EaMitsui)/k/(Temp+273))
    CAarray[0]=0.01325*Humidity #g/cm^3 moisture in the outer most layer of PET
    CAarray[1:intPETstop]=CAarray_past[1:intPETstop]+DAPET*delt*(CAarray_past[1-1:intPETstop-1]-2*CAarray_past[1:intPETstop]+CAarray_past[1+1:intPETstop+1])/(delx**2)
    CAarray[intPETstop] = CAarray[intPETstop-1] # derivitive = 0 boundary condition at boundary of two materials
    CAarray[intPETstop+1:xSteps-1]=CAarray_past[intPETstop+1:xSteps-1]+DAEVA*delt*(CAarray_past[intPETstop:xSteps-2]-2*CAarray_past[intPETstop+1:xSteps-1]+CAarray_past[intPETstop+2:xSteps])/(delx**2)
    CAarray[xSteps-1] = CAarray[xSteps-2] # derivitive = 0 boundary condition at silicon
    CAarray_past[0:xSteps]=CAarray[0:xSteps]
#    TRHarray[j,3]= CAarray[xSteps-1] # moisture content right next to cell
    if j % int(tSteps/3600.0/100.0)==0:
        weathercurrentarray[int(j/360000.0),7] = CAarray[0] # boundary condition
        weathercurrentarray[int(j/360000.0),8] = CAarray[xSteps-1] #moisture near cell
        weathercurrentarray[int(j/360000.0),9] = numpy.round(j/360000.0) #time in hours
        weathercurrentarray[int(j/360000.0),10] = DAPETarray[int(j/360000.0)] #diffusion coefficient for  PET
        weathercurrentarray[int(j/360000.0),11] = DAEVAarray[int(j/360000.0)] #diffusion coefficient for encap (Mitsui)
    if j % int(tSteps/10000.0)==0:
        print 'Progress ',100.0*j/tSteps
#        plt.plot(xStepsarray[:], CAarray[:])
    j = j+1
        

numpy.savetxt('C:\\Users\\khan\\Documents\\GitHub\\DiffusionModel\\Diffusion'+first_time_seconds+'.txt',weathercurrentarray)