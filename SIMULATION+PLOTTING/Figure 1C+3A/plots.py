#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import pylab as plt
import numpy as np
import seaborn as sns
sns.set(font_scale=2.0)

Tsim=5000

#load data from multimeter
#DBC0
VmsDBC0 = np.loadtxt('VmsDBC0.txt', dtype=float)
ts_mul_DBC0 = np.loadtxt('ts_mul_DBC0.txt', dtype=float)

#load data from spike detectors
#HYPERCOLUMN0
evsDBC0 = np.loadtxt('evsDBC0.txt', dtype=float)
tsDBC0 = np.loadtxt('tsDBC0.txt', dtype=float)
spikesDBC0 = np.loadtxt('spikesDBC0.txt', dtype=float)

evsPYR0 = np.loadtxt('evsPYR0.txt', dtype=float)
tsPYR0 = np.loadtxt('tsPYR0.txt', dtype=float)
spikesPYR0 = np.loadtxt('spikesPYR0.txt', dtype=float)

evsBS_HC0 = np.loadtxt('evsBS_HC0.txt', dtype=float)
tsBS_HC0 = np.loadtxt('tsBS_HC0.txt', dtype=float)
spikesPYR1 = np.loadtxt('spikesPYR1.txt', dtype=float)

evsDBC1 = np.loadtxt('evsDBC1.txt', dtype=float)
tsDBC1 = np.loadtxt('tsDBC1.txt', dtype=float)

evsPYR1 = np.loadtxt('evsPYR1.txt', dtype=float)
tsPYR1 = np.loadtxt('tsPYR1.txt', dtype=float)

#HYPERCOLUMN1
evsDBC2 = np.loadtxt('evsDBC2.txt', dtype=float)
tsDBC2 = np.loadtxt('tsDBC2.txt', dtype=float)

evsPYR2 = np.loadtxt('evsPYR2.txt', dtype=float)
tsPYR2 = np.loadtxt('tsPYR2.txt', dtype=float)
#spikesPYR2 = np.loadtxt('spikesPYR2.txt', dtype=float)

evsBS_HC1 = np.loadtxt('evsBS_HC1.txt', dtype=float)
tsBS_HC1 = np.loadtxt('tsBS_HC1.txt', dtype=float)

evsDBC3 = np.loadtxt('evsDBC3.txt', dtype=float)
tsDBC3 = np.loadtxt('tsDBC3.txt', dtype=float)

evsPYR3 = np.loadtxt('evsPYR3.txt', dtype=float)
tsPYR3 = np.loadtxt('tsPYR3.txt', dtype=float)



# Spiking activity in HC0 (Fig. 3A)
#Neuron ID -> y axes
#### HYPERCOLUMN 0 
### MINICOLUMN 0
#DBC_MC0       Neuron ID 1
#PYR_MC0       Neurons ID 2-31
### SHARED BASKETCELLS BETWEEN MC0 AND MC1]
#BS_HC0        Neurons ID 32-35
### MINICOLUMN 1    
#DBC_MC1       Neuron ID 36
#PYR_MC1       Neurons ID 37-66  
#### HYPERCOLUMN 1
### MINICOLUMN 2
#DBC_MC2       Neuron ID 67
#PYR_MC2       Neurons ID 68-97
###SHARED BASKETCELLS
#BS_HC1        Neurons ID 98-101
### MINICOLUMN 3        
#DBC_MC3       Neuron ID 102
#PYR_MC3       Neurons ID 103-132


sns.set_style("darkgrid")
plt.figure(1)
plt.figure(figsize=(12,10))
plt.title('Spike raster of neurons in HC0')
sns.regplot(x=tsDBC0, y=evsDBC0,fit_reg=False,scatter_kws={"s": 12},color='limegreen')
sns.regplot(x=tsPYR0, y=evsPYR0,fit_reg=False,scatter_kws={"s": 15},color='tomato')
sns.regplot(x=tsBS_HC0, y=evsBS_HC0,fit_reg=False,scatter_kws={"s": 12})
sns.regplot(x=tsDBC1, y=evsDBC1,fit_reg=False,scatter_kws={"s": 12},color='limegreen')
sns.regplot(x=tsPYR1, y=evsPYR1,fit_reg=False,scatter_kws={"s": 15},color='tomato')
plt.xlabel('time [ms]')

plt.ylabel('Neuron ID')
plt.savefig('Fig.3A.png')
plt.show()



### Spiking activity in HC1
#plt.figure(2)
#plt.figure(figsize=(15,13))
#plt.title('Spiking activity in HC1')
#sns.regplot(x=tsDBC2, y=evsDBC2,fit_reg=False,scatter_kws={"s": 10},color='limegreen')
#sns.regplot(x=tsPYR2, y=evsPYR2,fit_reg=False,scatter_kws={"s": 10},color='tomato')
#sns.regplot(x=tsBS_HC1, y=evsBS_HC1,fit_reg=False,scatter_kws={"s": 10})
#sns.regplot(x=tsDBC3, y=evsDBC3,fit_reg=False,scatter_kws={"s": 10},color='limegreen')
#sns.regplot(x=tsPYR3, y=evsPYR3,fit_reg=False,scatter_kws={"s": 10},color='tomato')
#plt.xlabel('time [ms]')
#plt.ylabel('Neuron ID')
#plt.show()


params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

# DBC0 membrane voltage (Fig. 1C)
plt.figure(2)
plt.figure(figsize=(12,7))
plt.plot(ts_mul_DBC0,VmsDBC0,color = "g")
plt.ylabel("$V_m$ $(DBC_{MC0}^{HC0})$ [mV]")
plt.xlabel("time [ms]")
plt.title("Membrane voltage of a stimulated DBC")
plt.savefig('Fig.1C.png')
plt.show()


