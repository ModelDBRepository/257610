#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import pylab as plt
import numpy as np
import seaborn as sns

##AMPA_weights
##PYR0_PYR0
w_AMPA_PYR0_PYR0_1s = np.loadtxt('w_AMPA_PYR0_PYR0_1s.txt', dtype=float)
w_AMPA_PYR0_PYR0_4s = np.loadtxt('w_AMPA_PYR0_PYR0_4s.txt', dtype=float)

##PYR0_PYR2
w_AMPA_PYR0_PYR2_1s = np.loadtxt('w_AMPA_PYR0_PYR2_1s.txt', dtype=float)
w_AMPA_PYR0_PYR2_4s = np.loadtxt('w_AMPA_PYR0_PYR2_4s.txt', dtype=float)


#PYR0_DBC1
w_AMPA_PYR0_DBC1_1s = np.loadtxt('w_AMPA_PYR0_DBC1_1s.txt', dtype=float)
w_AMPA_PYR0_DBC1_4s = np.loadtxt('w_AMPA_PYR0_DBC1_4s.txt', dtype=float)

#PYR0_DBC3
w_AMPA_PYR0_DBC3_1s = np.loadtxt('w_AMPA_PYR0_DBC3_1s.txt', dtype=float)
w_AMPA_PYR0_DBC3_4s = np.loadtxt('w_AMPA_PYR0_DBC3_4s.txt', dtype=float)


### load weights from microcircuit without DBCs
##PYR0_PYR0
w_AMPA_PYR0_PYR0_1s_NoDBC = np.loadtxt('w_AMPA_PYR0_PYR0_1s_NoDBC.txt', dtype=float)
w_AMPA_PYR0_PYR0_4s_NoDBC = np.loadtxt('w_AMPA_PYR0_PYR0_4s_NoDBC.txt', dtype=float)

##PYR0_PYR2
w_AMPA_PYR0_PYR2_1s_NoDBC = np.loadtxt('w_AMPA_PYR0_PYR2_1s_NoDBC.txt', dtype=float)
w_AMPA_PYR0_PYR2_4s_NoDBC = np.loadtxt('w_AMPA_PYR0_PYR2_4s_NoDBC.txt', dtype=float)


#PYR0_DBC1
w_AMPA_PYR0_PYR1_1s_NoDBC = np.loadtxt('w_AMPA_PYR0_PYR1_1s_NoDBC.txt', dtype=float)
w_AMPA_PYR0_PYR1_4s_NoDBC = np.loadtxt('w_AMPA_PYR0_PYR1_4s_NoDBC.txt', dtype=float)

#PYR0_DBC3
w_AMPA_PYR0_PYR3_1s_NoDBC = np.loadtxt('w_AMPA_PYR0_PYR3_1s_NoDBC.txt', dtype=float)
w_AMPA_PYR0_PYR3_4s_NoDBC = np.loadtxt('w_AMPA_PYR0_PYR3_4s_NoDBC.txt', dtype=float)

sns.set_style("darkgrid")
#sns.set_style("whitegrid")

sns.set(font_scale=2.0)
def plot_weightsSeaBornTest2(var_1000ms,var_4000ms,var_1000msNoDBC,var_4000msNoDBC,title,figNum):
    f, (ax1, ax2,ax3) = plt.subplots(3,figsize=(11,6),sharex=True,gridspec_kw={"height_ratios": (.10, .10, .85)})
    ax1.set_title(title)
    ax1.set_title(title, pad=20,fontsize=20)
    
    
    royalblue = dict(markerfacecolor='royalblue',marker='o',markeredgecolor='royalblue', markersize=5,linestyle='none')
    lime = dict(markerfacecolor='lime', marker='o',markeredgecolor='lime',markersize=5,linestyle='none')
    sns.boxplot(var_1000ms,color='royalblue',flierprops=royalblue,ax=ax1)
    sns.boxplot(var_4000ms,color='lime',flierprops=lime,ax=ax1)
    
    darkmagenta = dict(markerfacecolor='darkmagenta',marker='o',markeredgecolor='darkmagenta', markersize=5,linestyle='none')
    tomato = dict(markerfacecolor='tomato',marker='o',markeredgecolor='tomato', markersize=5,linestyle='none')
    sns.boxplot(var_1000msNoDBC,color='orange',flierprops=tomato,ax=ax2)
    sns.boxplot(var_4000msNoDBC,color='darkmagenta',flierprops=darkmagenta,ax=ax2)

    
    sns.kdeplot(var_1000ms, shade=True,kernel='triw',label='IWD-new model',ax=ax3,legend='weights') 
    sns.kdeplot(var_4000ms, shade=True,kernel='triw',label='LWD-new model',color='lime',ax=ax3) 

    plt.ylabel('Number of synaptic weights')
    plt.xlabel('Plastic weights')
    sns.kdeplot(var_1000msNoDBC, shade=True,kernel='triw',label='IWD-previous model',color='orange',ax=ax3,legend='weights') 
    sns.kdeplot(var_4000msNoDBC, shade=True,kernel='triw',label='LWD-previous model',color='darkmagenta',ax=ax3) 


    plt.subplots_adjust(wspace=None, hspace=0.05)
    plt.savefig('{}.png'.format(figNum))

    plt.show()

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

# new model (NM), previous model (PM)
plot_weightsSeaBornTest2(w_AMPA_PYR0_PYR0_1s,w_AMPA_PYR0_PYR0_4s,w_AMPA_PYR0_PYR0_1s_NoDBC,w_AMPA_PYR0_PYR0_4s_NoDBC,"$PYR_{MC0}^{HC0}$ to $PYR_{MC0}^{HC0}$","Fig. 2A")
plot_weightsSeaBornTest2(w_AMPA_PYR0_PYR2_1s,w_AMPA_PYR0_PYR2_4s,w_AMPA_PYR0_PYR2_1s_NoDBC,w_AMPA_PYR0_PYR2_4s_NoDBC,"$PYR_{MC0}^{HC0}$ to $PYR_{MC2}^{HC1}$","Fig. 2B")
plot_weightsSeaBornTest2(w_AMPA_PYR0_DBC1_1s,w_AMPA_PYR0_DBC1_4s,w_AMPA_PYR0_PYR1_1s_NoDBC,w_AMPA_PYR0_PYR1_4s_NoDBC,"$PYR_{MC0}^{HC0}$ to $DBC_{MC1}^{HC0}$ (new model), $PYR_{MC0}^{HC0}$ to $PYR_{MC1}^{HC0}$ (previous model)","Fig. 2C")
plot_weightsSeaBornTest2(w_AMPA_PYR0_DBC3_1s,w_AMPA_PYR0_DBC3_4s,w_AMPA_PYR0_PYR3_1s_NoDBC,w_AMPA_PYR0_PYR3_4s_NoDBC,"$PYR_{MC0}^{HC0}$ to $DBC_{MC3}^{HC1}$ (new model), $PYR_{MC0}^{HC0}$ to $PYR_{MC3}^{HC1}$ (previous model) ","Fig. 2D")
