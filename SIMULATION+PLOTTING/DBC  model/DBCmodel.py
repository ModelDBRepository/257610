#!/usr/bin/env python2
# coding: utf-8

import numpy as np, pylab as plt
import nest
import seaborn as sns; sns.set()
import sys
sys.path.insert(0,'/home/nik/Documents/BCPNN_NEST_Module') #Python checks and inserts the new directory
import BCPNN # 'pt_module'

nest.ResetKernel()
BCPNN.InstallBCPNN()
nest.SetKernelStatus({'resolution':0.001})
sns.set(font_scale=1.7)


syn_ports = {'AMPA':1,'NMDA':2,'GABA':3} #receptor types
f_desired=7.5
f_max=55.

#In the DBCmodel.py, the parameters we used to achieve satisfactory electrophysiological fidelity are included.
#The simulations aim at reproducing spike patterns under sweeps of increasing  suprathreshold current 
#steps (10 pA each) and other reported activity. The range of the stimulation input current is on the same 
#level with the one reported in the paper below.

#The spike patterns produced (figure DBC_ActivityPatterns) can be directly compared with the findings of fig.4B appeared in Cluster 
#analysis–Based Physiological Classification and Morphological Properties of Inhibitory 
#Neurons in Layers 2–3 of Monkey Dorsolateral Prefrontal Cortex (Krimer et al., 2005).


# DBC parameters based on reported results and tuning 
NRN={
         'cell_model': 'aeif_cond_exp_multisynapse',
         'neuron_params': {
          'AMPA_NEG_E_rev': -75.0,#pseudo-negative reversal potential used for negative BCPNN weights
          'AMPA_Tau_decay': 5.0,#synaptic time constant
          'Delta_T': 1.0,
          'E_L': -76.0,#Leak Reversal Potention
          'GABA_E_rev': -80.0,
          'GABA_Tau_decay': 10.0,
          'NMDA_NEG_E_rev': -75.0,
          'NMDA_Tau_decay': 100.0,
          'V_reset': -52.0,#Reset Potential
          'V_th': -44.0, #Spike Threshold
          'a': 0.0,  #subthreshold  adaptation
          'b': 3.0, #spike adaptation in [pA]
          'bias': np.log(f_desired/f_max), #initial BCPNN bias 
          'epsilon': 0.01,#BCPNN epsilon
          'fmax': 20.0, #BCPNN fmax
          'g_L':1.52, #leak conductance
          'gain': 0., #BCPNN bias gain. Should be set such that noise activity matches f_desired. Leads to zero mean weights
          'gsl_error_tol': 1e-12,
          'kappa': 0.0,#BCPNN plasticity switch
          'p_j': f_desired/f_max, #BCPNN pj trance
          't_ref': 2.0,
          'tau_e': 0.5,#BCPNN time constant
          'tau_j': 5.0,#BCPNN time constant
          'tau_p': 5000.0,#BCPNN learning time constant
          'tau_w': 200.0,#adaptation time constant
          'w': 0.0}}


DBC_params=NRN['neuron_params']
if 'DBC' not in nest.Models():
    nest.CopyModel(NRN['cell_model'],'DBC',DBC_params)  #create parameterized L23e_cell for use later on


#The simulations aim at reproducing spike patterns under sweeps of increasing 
#suprathreshold current steps (10 pA each) and other reported activity. The 
#range of the stimulation input current is on the same level with the one 
#reported in the paper below.
    
#The spike patterns produced can be directly compared with the findings of fig.4B appeared in Cluster 
#analysis–Based Physiological Classification and Morphological Properties of Inhibitory 
#Neurons in Layers 2–3 of Monkey Dorsolateral Prefrontal Cortex (Krimer et al., 2005).


I=[47.5,57.5,67.5,77.5]

# DBC membrane capacitance tuned to 15 nS
C_m = 15.0

# DBCs receive different input current steps
DBC0=nest.Create("DBC")
nest.SetStatus(DBC0, {"I_e": I[0]})
nest.SetStatus(DBC0, {"C_m": C_m})

DBC1=nest.Create("DBC")
nest.SetStatus(DBC1, {"I_e": I[1]})
nest.SetStatus(DBC1, {"C_m": C_m})

DBC2=nest.Create("DBC")
nest.SetStatus(DBC2, {"I_e": I[2]})
nest.SetStatus(DBC2, {"C_m": C_m})

DBC3=nest.Create("DBC")
nest.SetStatus(DBC3, {"I_e": I[3]})
nest.SetStatus(DBC3, {"C_m": C_m})


# Create and connect each multimeter with the corresponding DBC
multimeter0 = nest.Create("multimeter")
nest.SetStatus(multimeter0,{"withtime":True,"interval":0.001,"record_from":["V_m"]})
nest.Connect(multimeter0,DBC0)

multimeter1 = nest.Create("multimeter")
nest.SetStatus(multimeter1,{"withtime":True,"interval":0.001,"record_from":["V_m"]})
nest.Connect(multimeter1,DBC1)

multimeter2 = nest.Create("multimeter")
nest.SetStatus(multimeter2,{"withtime":True,"interval":0.001,"record_from":["V_m"]})
nest.Connect(multimeter2,DBC2)

multimeter3 = nest.Create("multimeter")
nest.SetStatus(multimeter3,{"withtime":True,"interval":0.001,"record_from":["V_m"]})
nest.Connect(multimeter3,DBC3)

# Simulation in the afformentioned paper (fig.4B) lasts roughly 470 ms
Tsim=470
nest.Simulate(float(Tsim))

# Returns membrane voltage with 0.001 resolution
def membVolt(multimeter):
    dmm = nest.GetStatus(multimeter)[0]
    Vms = dmm["events"]["V_m"]
    ts = dmm["events"]["times"]
    return ts,Vms

# Membrane voltage diagrams according to the 4 inpute current steps.
fig, axs = plt.subplots(2, 2)
fig.set_size_inches(14.5, 8.5)
fig.subplots_adjust(hspace=0.45, wspace=0.2)
axs[0, 0].plot(membVolt(multimeter0)[0], membVolt(multimeter0)[1])
axs[0, 0].set_title('Stimulation current, 47.5 pA')
axs[0, 1].plot(membVolt(multimeter1)[0], membVolt(multimeter1)[1], 'tab:orange')
axs[0, 1].set_title('Stimulation current, 57.5 pA')
axs[1, 0].plot(membVolt(multimeter2)[0], membVolt(multimeter2)[1], 'tab:green')
axs[1, 0].set_title('Stimulation current, 67.5 pA')
axs[1, 1].plot(membVolt(multimeter3)[0], membVolt(multimeter3)[1], 'tab:red')
axs[1, 1].set_title('Stimulation current, 77.5 pA')

for ax in axs.flat:
    ax.set(xlabel='time [ms]', ylabel='Vm (DBC) [mV]')
fig.savefig('DBC_ActivityPatterns.png')

