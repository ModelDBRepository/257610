#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pylab as plt
import numpy as np
import seaborn as sns
sns.set_style("darkgrid")
sns.set(font_scale=2)


# microcircuit incorporating DBCs
x=0
I_GABA_PYR0_list=[]  # I_GABA_PYR0 values, total inhibitory input current received by pyramidal cells in MC0 
ts_GABA_PYR0_list=[]
PYR_MC0=30 #number of pyramidal cells of MC0
iterations=100
for x in range(PYR_MC0*iterations):  # for every neuron at each iteration
    I_GABA_PYR0_list.append(np.loadtxt('I_GABA_PYR0_list{}.txt'.format(x), dtype=float))
    ts_GABA_PYR0_list.append(np.loadtxt('ts_GABA_PYR0_list{}.txt'.format(x), dtype=float))
    
I_GABA_3to4_sec=[]
for i in range(len(I_GABA_PYR0_list)):
    a=np.array_split(I_GABA_PYR0_list[i], 2)  # split list and keep the values from 2500ms to 5000ms
    I_GABA_3to4_sec.append(a[1])  # a[1] includes values from 2500ms to 5000ms and we concentrate between 3000ms and 4000ms later on 

timeWindow=10  #the total inhibitory input current is collected every 10ms 
down=0
up=timeWindow

data=[]
keepValues=[]
itera=0
for itera in range(2000/timeWindow):  # from 2500ms to 4500ms divided according to a time window e.g 50th iteration corresponds to 3000ms and 200th to 4500ms 
    for kappa in range(len(I_GABA_3to4_sec)):
        sample = I_GABA_3to4_sec[kappa][down:up]
        keepValues.extend(sample)  #collect values every 10ms 
    data.append(keepValues) # stored values from 2500ms to 4500ms 
    keepValues=[]
    down=down+timeWindow
    up=up+timeWindow


data=np.array(data)
data=data.transpose() # for plotting purposes 
    

# microcircuit without DBCs
# the following code is dublicated for the model without DBCs 
PYR_MC0=30
x=0
I_GABA_PYR0_NoDBC_list=[]
ts_GABA_PYR0_NoDBC_list=[]
for x in range(PYR_MC0*iterations):
    I_GABA_PYR0_NoDBC_list.append(np.loadtxt('I_GABA_PYR0_NoDBC_list{}.txt'.format(x), dtype=float))
    ts_GABA_PYR0_NoDBC_list.append(np.loadtxt('ts_GABA_PYR0_NoDBC_list{}.txt'.format(x), dtype=float))
    
I_GABA_3to4_NoDBC_sec=[]
lamda=0
for lamda in range(len(I_GABA_PYR0_NoDBC_list)):
    a_NoDBC=np.array_split(I_GABA_PYR0_NoDBC_list[lamda], 2)
    I_GABA_3to4_NoDBC_sec.append(a_NoDBC[1])
    
    

down=0
up=timeWindow
data_NoDBC=[]
keepValues_NoDBC=[]
itera=0
for itera in range(2000/timeWindow):
    for kappa_NoDBC in range(len(I_GABA_3to4_NoDBC_sec)):
        sample_NoDBC = I_GABA_3to4_NoDBC_sec[kappa_NoDBC][down:up]
        keepValues_NoDBC.extend(sample_NoDBC)
    data_NoDBC.append(keepValues_NoDBC)
    keepValues_NoDBC=[]
    down=down+timeWindow
    up=up+timeWindow


data_NoDBC=np.array(data_NoDBC)
data_NoDBC=data_NoDBC.transpose()

#The diagram shows the population averaged total inhibitory input current received by pyramidal cells in MC0 in both architectures between 2500ms and 4500ms
plt.figure(1)
plt.figure(figsize=(12,10))
sns.tsplot(data=data,time=range(2500,4500,10), ci="sd",interpolate=False,condition="New cortical model",ms=7,color='darkmagenta',legend='WTA with DBCs')
sns.tsplot(data=data_NoDBC,time=range(2500,4500,10),ci="sd",interpolate=False,condition="Previous cortical model",ms=7,color='royalblue',legend='WTA')
plt.ylabel('$I_{GABA}$ [pA]')
plt.xlabel('time [ms]')
plt.ylim(-500,300)
plt.title('Functionality verification')
plt.savefig('Fig.3B.png')
