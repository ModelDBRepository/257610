#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import time
import numpy as np
import numpy as np, pylab as plt
import nest
import sys
sys.path.insert(0,'/home/nik/Documents/BCPNN_NEST_Module') #Python checks and inserts the new directory
import BCPNN # 'pt_module'



# PYR0_PYR1 , PYR0_PYR3
w_AMPA_PYR0_PYR1_1s_NoDBC,w_AMPA_PYR0_PYR3_1s_NoDBC=([] for i in range(2))
i=0
w_AMPA_PYR0_PYR1_4s_NoDBC,w_AMPA_PYR0_PYR3_4s_NoDBC=([] for i in range(2))
i=0
       
#PYR0_PYR0  , PYR0_PYR2
w_AMPA_PYR0_PYR0_1s_NoDBC,w_AMPA_PYR0_PYR2_1s_NoDBC=([] for i in range(2))
i=0
w_AMPA_PYR0_PYR0_4s_NoDBC,w_AMPA_PYR0_PYR2_4s_NoDBC=([] for i in range(2))
i=0   

#multiple rounds 
for iter in range(10):
    nest.ResetKernel()
    nest.SetKernelStatus({'resolution':0.1})
    seed= int( time.time() * 1000.0 )
    nest.SetKernelStatus({'rng_seeds': [seed]})

    BCPNN.InstallBCPNN()
    
    syn_ports = {'AMPA':1,'NMDA':2,'GABA':3} #receptor types
    f_desired=1.
    f_max=20.
    
    
    
    NRN_L23e={
             'cell_model': 'aeif_cond_exp_multisynapse',
             'neuron_params': {
              'AMPA_NEG_E_rev': -75.0,#pseudo-negative reversal potential used for negative BCPNN weights
              'AMPA_Tau_decay': 5.0,#synaptic time constant
              'Delta_T': 3.0,
              'E_L': -70.0,#Leak Reversal Potention
              'GABA_E_rev': -75.0,
              'GABA_Tau_decay': 5.0,
              'NMDA_NEG_E_rev': -75.0,
              'NMDA_Tau_decay': 100.0,
              'V_reset': -80.0,#Reset Potential
              'V_th': -55.0, #Spike Threshold
              'a': 0.0,  #subthreshold  adaptation
              'b': 86.0, #spike adaptation in [pA]
              'bias': np.log(f_desired/f_max), #initial BCPNN bias 
              'epsilon': 0.01,#BCPNN epsilon
              'fmax': 20.0, #BCPNN fmax
              'g_L': 14.0, #leak conductance
              'gain': 31.7,#65.,#31.7, #BCPNN bias gain. Should be set such that noise activity matches f_desired. Leads to zero mean weights
              'gsl_error_tol': 1e-12,
              'kappa': 1.0,#BCPNN plasticity switch
              'p_j': f_desired/f_max, #BCPNN pj trance
              't_ref': 5.0,
              'tau_e': 0.1,#BCPNN time constant
              'tau_j': 5.0,#BCPNN time constant
              'tau_p': 5000.0,#BCPNN learning time constant
              'tau_w': 500.0,#adaptation time constant
              'w': 0.0}}
    
    
    
    ST={
             'L23e_zmn_rate': 750.0, 
             'STIM0_rate':1700.0,
             'STIM1_rate':1700.0,
             'stim_length': 250.0,
             'stim_rate': 1700.0,
             'stim_weight': 2.0,
             'zmn_rate': 950.0,
             'zmn_delay': 0.1,
             'zmn_weight': 1.5}
    
    SYN={
             'AMPA_synapse_param': {
                'K': 1.0, #BCPNN plasticity switch
                'U': 0.25,#vescicle depletion/synaptic depression (markram_tsodyks type)
                'bias': np.log(f_desired/f_max),#initial BCPNN bias (computed at runtime)
                'delay': 1.0,#synaptic conductance delay
                'epsilon': 0.01,#BCPNN epsilon
                'fmax': 20.0,#BCPNN fmax
                'gain': 3.92,#5.3,#3.92,#BCPNN synaptic gain
                'p_i':  f_desired/f_max,#BCPNN p trace (presynaptic)
                'p_ij':(f_desired/f_max)**2,#BCPNN p trace (joint)
                'p_j':  f_desired/f_max, #BCPNN p trace (postsynaptic)
                'receptor_type': 1, #1=AMPA
                'stp_flag': 1.0, #STP enabled
                't_k': 0.0, #reset this when switching plasticity via k/kappa
                'tau_e': 0.1, #BCPNN time constant 
                'tau_fac': 0.0,#Facilitation time constant (0->no facilitation)
                'tau_i': 5.0,#BCPNN time constant
                'tau_j': 5.0,#BCPNN time constant
                'tau_p': 5000.0,#BCPNN time constant
                'tau_rec': 500.0,#depression time constant (markram_tsodyks type)
                'u': 0.25,#vescicle depletion/synaptic depression (markram_tsodyks type)
                'weight': 0.0,#default weight (computed at runtime)
                'x': 0.25},#depression variable (markram_tsodyks type)
             'delay_min': 1.,#minimum synaptic delay in the grid
             'delay_max': 7.,#minimum synaptic delay across the grid
             'delay_eie': 1.5,#feedback inhibition connection delay
             'e2i_weight': 3.5,#feedback inhibition connection weight
             'fmax': 20.0,#BCPNN fmax
             'gain': 3.92,    #BCPNN gain
             'i2e_weight': -35.0,
             'prob_e2e': 0.2,#recurrent connection probability
             'prob_e2i': 0.7,#feedback inhibition connection probability
             'prob_i2e': 0.7,#feedback inhibition connection probability
             'receptors': ['GABA'],#receptors used (no GABA/NMDA in this simplified version)
             'synapse_model': 'bcpnn_synapse',
             'tau_AMPA': 5.0,#synaptic time constant
             'tau_NMDA': 100.0,#synaptic time constant
             'tau_p': 5000.0,#BCPNN learning time constant
             'tau_rec': 500.0,#adaptation time constant
             'tau_w': 500.0}
    
    
    
    L23e_cell_params=NRN_L23e['neuron_params']
    if 'PYR' not in nest.Models():
        nest.CopyModel(NRN_L23e['cell_model'],'PYR',L23e_cell_params)  #create parameterized L23e_cell (pyramidal cell) for use later on
    
    
    basket_cell_params=NRN_L23e['neuron_params']
    basket_cell_params.update(b = 0.0, gain= 0.0) #basket cells have no neural plasticity.
    if 'basket_cell' not in nest.Models(): 
        nest.CopyModel(NRN_L23e['cell_model'],'basket_cell',basket_cell_params) #basket cell
    
    if 'stim_synapse' not in nest.Models():
        nest.CopyModel('static_synapse','stim_synapse',{'weight':ST['stim_weight'],'delay': 0.1,'receptor_type': syn_ports['AMPA']})
    if 'zmn_synapse' not in nest.Models():
        nest.CopyModel('static_synapse','zmn_synapse',{'weight':ST['zmn_weight'],'delay': 0.1,'receptor_type': syn_ports['AMPA']})       
    if 'e2i_synapse' not in nest.Models():
        nest.CopyModel('static_synapse','e2i_synapse',{'weight':SYN['e2i_weight'],'delay': SYN['delay_eie'],'receptor_type':syn_ports['AMPA']})
    if 'i2e_synapse' not in nest.Models():
        nest.CopyModel('static_synapse','i2e_synapse',{'weight':SYN['i2e_weight'],'delay': SYN['delay_eie'],'receptor_type':syn_ports['GABA']})
    #BCPNN synapse (AMPA)
    if 'AMPA_synapse' not in nest.Models():
        nest.CopyModel(SYN['synapse_model'],'AMPA_synapse',SYN['AMPA_synapse_param'])
    
    
    nest.CopyModel("e2i_synapse","e2i_PYR_BS",{"delay": 1.5})
    nest.CopyModel("AMPA_synapse","AMPA_synapse_short_delay",{"delay": 1.5}) #delay within HCs
    nest.CopyModel("AMPA_synapse","AMPA_synapse_larger_delay",{"delay": 4.5}) #delay between HCs
    
    
    
    conn_dict_PYR_TO_PYR = {'rule': 'fixed_total_number', 'N': 180}  #20% cp (connection probability)
    conn_dict_PYR_TO_BS = {'rule': 'fixed_total_number', 'N': 84}  #70% cp
    conn_dict_PYR_TO_DBC={'rule': 'fixed_total_number', 'N': 6} #20% cp
    conn_dict_PYR_TO_PYR_DIFF_HC={'rule': 'fixed_total_number', 'N': 180} #20% cp
    
    ### HYPERCOLUMN 0 
    ### MINICOLUMN 0
    PYR_MC0=nest.Create("PYR", 30)             #nodes 2-31
    ### SHARED BASKETCELLS BETWEEN MC0 AND MC1]
    BS_HC0 =nest.Create("basket_cell", 4)      #nodes 32-35
    ### MINICOLUMN 1    
    
    PYR_MC1=nest.Create("PYR", 30)             #nodes 37-66  
    ###### HYPERCOLUMN 1
    ### MINICOLUMN 2
    
    PYR_MC2=nest.Create("PYR", 30)             #68-97
    ###SHARED BASKETCELLS
    BS_HC1 =nest.Create("basket_cell", 4)      #nodes 98-101
    ### MINICOLUMN 3        
    
    PYR_MC3=nest.Create("PYR", 30)             #nodes 103-132
    
    
    ### CONNECTIONS 
    #MINICOLUMN 0
    # DBC_TO_PYR
    #PYR_TO_BS
    nest.Connect(PYR_MC0, BS_HC0, conn_dict_PYR_TO_BS,syn_spec = {'model':'e2i_synapse'})
    #BS_TO_PYR
    nest.Connect(BS_HC0, PYR_MC0, conn_dict_PYR_TO_BS,syn_spec = {'model':'i2e_synapse'})
    #PYR_TO_PYR(connection to itself)
    nest.Connect(PYR_MC0, PYR_MC0, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_short_delay'})
    #PYR0_TO_DBC1_DBC3
    nest.Connect(PYR_MC0, PYR_MC1, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_short_delay'})
    nest.Connect(PYR_MC0, PYR_MC3, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_larger_delay'})
    #PYR_TO_PYR(DIFFERENT HYPERCOLUMNS)
    nest.Connect(PYR_MC0, PYR_MC2, conn_dict_PYR_TO_PYR_DIFF_HC,syn_spec = {'model':'AMPA_synapse_larger_delay'})
    
    
    
    #MINICOLUMN 1
    # DBC_TO_PYR
    #PYR_TO_BS
    nest.Connect(PYR_MC1, BS_HC0, conn_dict_PYR_TO_BS,syn_spec = {'model':'e2i_synapse'})
    #BS_TO_PYR
    nest.Connect(BS_HC0, PYR_MC1, conn_dict_PYR_TO_BS,syn_spec = {'model':'i2e_synapse'})
    #PYR_TO_PYR(connection to itself)
    nest.Connect(PYR_MC1, PYR_MC1, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_short_delay'})
    #PYR1_TO_DBC0_DBC2
    nest.Connect(PYR_MC1, PYR_MC0, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_short_delay'})
    nest.Connect(PYR_MC1, PYR_MC2, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_larger_delay'})
    #PYR_TO_PYR(DIFFERENT HYPERCOLUMNS)
    nest.Connect(PYR_MC1, PYR_MC3, conn_dict_PYR_TO_PYR_DIFF_HC,syn_spec = {'model':'AMPA_synapse_larger_delay'})
    
    
    #MINICOLUMN 2
    # DBC_TO_PYR
    #PYR_TO_BS
    nest.Connect(PYR_MC2, BS_HC1, conn_dict_PYR_TO_BS,syn_spec = {'model':'e2i_synapse'})
    #BS_TO_PYR
    nest.Connect(BS_HC1, PYR_MC2, conn_dict_PYR_TO_BS,syn_spec = {'model':'i2e_synapse'})
    #PYR_TO_PYR(connection to itself)
    nest.Connect(PYR_MC2, PYR_MC2, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_short_delay'})
    #PYR2_TO_DBC1_DBC3
    nest.Connect(PYR_MC2, PYR_MC1, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_larger_delay'})
    nest.Connect(PYR_MC2, PYR_MC3, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_short_delay'})
    #PYR_TO_PYR(DIFFERENT HYPERCOLUMNS)
    nest.Connect(PYR_MC2, PYR_MC0, conn_dict_PYR_TO_PYR_DIFF_HC,syn_spec = {'model':'AMPA_synapse_larger_delay'})
    
    
    #MINICOLUMN 3
    # DBC_TO_PYR
    #PYR_TO_BS
    nest.Connect(PYR_MC3, BS_HC1, conn_dict_PYR_TO_BS,syn_spec = {'model':'e2i_synapse'})
    #BS_TO_PYR
    nest.Connect(BS_HC1, PYR_MC3, conn_dict_PYR_TO_BS,syn_spec = {'model':'i2e_synapse'})
    #PYR_TO_PYR(connection to itself)
    nest.Connect(PYR_MC3, PYR_MC3, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_short_delay'})
    #PYR3_TO_DBC0_DBC2
    nest.Connect(PYR_MC3, PYR_MC0, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_larger_delay'})
    nest.Connect(PYR_MC3, PYR_MC2, conn_dict_PYR_TO_PYR,syn_spec = {'model':'AMPA_synapse_short_delay'})
    #PYR_TO_PYR(DIFFERENT HYPERCOLUMNS)
    nest.Connect(PYR_MC3, PYR_MC1, conn_dict_PYR_TO_PYR_DIFF_HC,syn_spec = {'model':'AMPA_synapse_larger_delay'})
    
    
    
    # STIMULATIONS 
    # ZERO MEAN NOISE 
    zmn_nodes_L23e=nest.Create('poisson_generator', params={'rate'  : ST['L23e_zmn_rate']})
    zmn_nodes_L23i=nest.Create('poisson_generator', params={'rate'  : ST['L23e_zmn_rate']})
    nest.SetStatus(zmn_nodes_L23e, {"start": 0.0})
    nest.SetStatus(zmn_nodes_L23e, {"stop": 5000.0})
    
    nest.SetStatus(zmn_nodes_L23i, {"start": 0.0})
    nest.SetStatus(zmn_nodes_L23i, {"stop":5000.0})
    
    syn_dict_e = {'model': 'zmn_synapse', 'weight': +ST['zmn_weight'], 'delay': ST['zmn_delay']}
    syn_dict_i = {'model': 'zmn_synapse', 'weight': -ST['zmn_weight'], 'delay': ST['zmn_delay']}
    nest.DivergentConnect(zmn_nodes_L23e, PYR_MC0, model=syn_dict_e['model'],weight=syn_dict_e['weight'],delay=syn_dict_e['delay'])  #adjusted for nest 2.2.2 
    nest.DivergentConnect(zmn_nodes_L23i, PYR_MC0, model=syn_dict_i['model'],weight=syn_dict_i['weight'],delay=syn_dict_i['delay'])  #adjusted for nest 2.2.2 
    
    
    nest.DivergentConnect(zmn_nodes_L23e, PYR_MC1, model=syn_dict_e['model'],weight=syn_dict_e['weight'],delay=syn_dict_e['delay'])  #adjusted for nest 2.2.2 
    nest.DivergentConnect(zmn_nodes_L23i, PYR_MC1, model=syn_dict_i['model'],weight=syn_dict_i['weight'],delay=syn_dict_i['delay'])  #adjusted for nest 2.2.2 
    
    
    nest.DivergentConnect(zmn_nodes_L23e, PYR_MC2, model=syn_dict_e['model'],weight=syn_dict_e['weight'],delay=syn_dict_e['delay'])  #adjusted for nest 2.2.2 
    nest.DivergentConnect(zmn_nodes_L23i, PYR_MC2, model=syn_dict_i['model'],weight=syn_dict_i['weight'],delay=syn_dict_i['delay'])  #adjusted for nest 2.2.2 
    
    
    nest.DivergentConnect(zmn_nodes_L23e, PYR_MC3, model=syn_dict_e['model'],weight=syn_dict_e['weight'],delay=syn_dict_e['delay'])  #adjusted for nest 2.2.2 
    nest.DivergentConnect(zmn_nodes_L23i, PYR_MC3, model=syn_dict_i['model'],weight=syn_dict_i['weight'],delay=syn_dict_i['delay'])  #adjusted for nest 2.2.2 
    
    
    ## STIMO TO PYR0,PYR2
    STIM0=nest.Create('poisson_generator', params={'rate'  : ST['STIM0_rate']})
    nest.SetStatus(STIM0, {"start": 1000.0})
    nest.SetStatus(STIM0, {"stop": 2000.0})
    syn_dict_ex0 = {'model': 'zmn_synapse', 'weight': +ST['zmn_weight'], 'delay': ST['zmn_delay']}
    nest.DivergentConnect(STIM0, PYR_MC0, model=syn_dict_ex0['model'],weight=syn_dict_ex0['weight'],delay=syn_dict_ex0['delay'])  #adjusted for nest 2.2.2 
    nest.DivergentConnect(STIM0, PYR_MC2, model=syn_dict_ex0['model'],weight=syn_dict_ex0['weight'],delay=syn_dict_ex0['delay'])  #adjusted for nest 2.2.2 
    
    ## STIM1 TO PYR1,PYR3
    STIM1=nest.Create('poisson_generator', params={'rate'  : ST['STIM1_rate']})
    nest.SetStatus(STIM1, {"start": 3000.0})
    nest.SetStatus(STIM1, {"stop": 4000.0})
    
    syn_dict_ex1 = {'model': 'zmn_synapse', 'weight': +ST['zmn_weight'], 'delay': ST['zmn_delay']}
    
    nest.DivergentConnect(STIM1, PYR_MC1, model=syn_dict_ex1['model'],weight=syn_dict_ex1['weight'],delay=syn_dict_ex1['delay'])  #adjusted for nest 2.2.2 
    nest.DivergentConnect(STIM1, PYR_MC3, model=syn_dict_ex1['model'],weight=syn_dict_ex1['weight'],delay=syn_dict_ex1['delay'])  #adjusted for nest 2.2.2 
    
    
    
    #MULTIMETERS
    #HYPERCOLUMN0
    #DBC0
    interv=0.1
    
    #PYR0
    multimeter_PYR0 = nest.Create("multimeter")
    nest.SetStatus(multimeter_PYR0,{"withtime":True,"interval":interv,"record_from":["V_m"]})
    nest.ConvergentConnect(multimeter_PYR0,[PYR_MC0[0]])
    
    #BS_HC0
    multimeter_BS_HC0 = nest.Create("multimeter")
    nest.SetStatus(multimeter_BS_HC0,{"withtime":True,"interval":interv,"record_from":["V_m"]})
    nest.ConvergentConnect(multimeter_BS_HC0,[BS_HC0[0]])
    
        
    #SPIKEDETECTORS
    #HYPERCOLUMN0
    PYR_MC0_spikedetector = nest.Create("spike_detector",params={"withgid":True,"withtime":True})
    BS_HC0_spikedetector = nest.Create("spike_detector",params={"withgid":True,"withtime":True})
    PYR_MC1_spikedetector = nest.Create("spike_detector",params={"withgid":True,"withtime":True})
    
    #HYPERCOLUMN1
    PYR_MC2_spikedetector = nest.Create("spike_detector",params={"withgid":True,"withtime":True})
    BS_HC1_spikedetector = nest.Create("spike_detector",params={"withgid":True,"withtime":True})
    PYR_MC3_spikedetector = nest.Create("spike_detector",params={"withgid":True,"withtime":True})
    
    
    nest.ConvergentConnect(PYR_MC0,PYR_MC0_spikedetector)
    nest.ConvergentConnect(BS_HC0,BS_HC0_spikedetector)
    
    nest.ConvergentConnect(PYR_MC1,PYR_MC1_spikedetector)
    
    nest.ConvergentConnect(PYR_MC2,PYR_MC2_spikedetector)
    nest.ConvergentConnect(BS_HC1,BS_HC1_spikedetector)
    
    nest.ConvergentConnect(PYR_MC3,PYR_MC3_spikedetector)
    
    
    
    
    
    #PYR_MC0
    AMPA_PYR0_PYR0=nest.GetConnections(PYR_MC0,target=PYR_MC0,synapse_model="AMPA_synapse_short_delay")
    AMPA_PYR0_PYR2=nest.GetConnections(PYR_MC0,target=PYR_MC2,synapse_model="AMPA_synapse_larger_delay")
    AMPA_PYR0_PYR1=nest.GetConnections(PYR_MC0,target=PYR_MC1,synapse_model="AMPA_synapse_short_delay")
    AMPA_PYR0_PYR3=nest.GetConnections(PYR_MC0,target=PYR_MC3,synapse_model="AMPA_synapse_larger_delay")
    
    # PYR_MC1
    AMPA_PYR1_PYR1=nest.GetConnections(PYR_MC1,target=PYR_MC1,synapse_model="AMPA_synapse_short_delay")
    AMPA_PYR1_PYR3=nest.GetConnections(PYR_MC1,target=PYR_MC3,synapse_model="AMPA_synapse_larger_delay")
    AMPA_PYR1_PYR0=nest.GetConnections(PYR_MC1,target=PYR_MC0,synapse_model="AMPA_synapse_short_delay")
    AMPA_PYR1_PYR2=nest.GetConnections(PYR_MC1,target=PYR_MC2,synapse_model="AMPA_synapse_larger_delay")
    
    # PYR_MC2
    AMPA_PYR2_PYR2=nest.GetConnections(PYR_MC2,target=PYR_MC2,synapse_model="AMPA_synapse_short_delay")
    AMPA_PYR2_PYR0=nest.GetConnections(PYR_MC2,target=PYR_MC0,synapse_model="AMPA_synapse_larger_delay")
    AMPA_PYR2_PYR1=nest.GetConnections(PYR_MC2,target=PYR_MC1,synapse_model="AMPA_synapse_larger_delay")
    AMPA_PYR2_PYR3=nest.GetConnections(PYR_MC2,target=PYR_MC3,synapse_model="AMPA_synapse_short_delay")
    
    # PYR_MC3
    AMPA_PYR3_PYR3=nest.GetConnections(PYR_MC3,target=PYR_MC3,synapse_model="AMPA_synapse_short_delay")
    AMPA_PYR3_PYR1=nest.GetConnections(PYR_MC3,target=PYR_MC1,synapse_model="AMPA_synapse_larger_delay")
    AMPA_PYR3_PYR2=nest.GetConnections(PYR_MC3,target=PYR_MC2,synapse_model="AMPA_synapse_short_delay")
    AMPA_PYR3_PYR0=nest.GetConnections(PYR_MC3,target=PYR_MC0,synapse_model="AMPA_synapse_larger_delay")
    
    
    
    
    Tsim=5000 # 5 seconds
    for x in range (Tsim):
        nest.Simulate(1.)
        
        if (x==1000): # keep the weights at 1000 ms (before training)
            w_AMPA_PYR0_PYR1_1s_NoDBC.extend([np.log(a['p_ij']/(a['p_i']*a['p_j'])) for a in nest.GetStatus(AMPA_PYR0_PYR1)])
            w_AMPA_PYR0_PYR3_1s_NoDBC.extend([np.log(a['p_ij']/(a['p_i']*a['p_j'])) for a in nest.GetStatus(AMPA_PYR0_PYR3)])
    
            w_AMPA_PYR0_PYR0_1s_NoDBC.extend([np.log(a['p_ij']/(a['p_i']*a['p_j'])) for a in nest.GetStatus(AMPA_PYR0_PYR0)])
            w_AMPA_PYR0_PYR2_1s_NoDBC.extend([np.log(a['p_ij']/(a['p_i']*a['p_j'])) for a in nest.GetStatus(AMPA_PYR0_PYR2)])
    
        if (x==4000): # keep the weights at 4000 ms (after training)
           
            w_AMPA_PYR0_PYR1_4s_NoDBC.extend([np.log(a['p_ij']/(a['p_i']*a['p_j'])) for a in nest.GetStatus(AMPA_PYR0_PYR1)])
            w_AMPA_PYR0_PYR3_4s_NoDBC.extend([np.log(a['p_ij']/(a['p_i']*a['p_j'])) for a in nest.GetStatus(AMPA_PYR0_PYR3)])
    
            w_AMPA_PYR0_PYR0_4s_NoDBC.extend([np.log(a['p_ij']/(a['p_i']*a['p_j'])) for a in nest.GetStatus(AMPA_PYR0_PYR0)])
            w_AMPA_PYR0_PYR2_4s_NoDBC.extend([np.log(a['p_ij']/(a['p_i']*a['p_j'])) for a in nest.GetStatus(AMPA_PYR0_PYR2)])
    

#save weights before training
np.savetxt("w_AMPA_PYR0_PYR1_1s_NoDBC.txt", w_AMPA_PYR0_PYR1_1s_NoDBC, delimiter=",")
np.savetxt("w_AMPA_PYR0_PYR3_1s_NoDBC.txt", w_AMPA_PYR0_PYR3_1s_NoDBC, delimiter=",")
np.savetxt("w_AMPA_PYR0_PYR0_1s_NoDBC.txt", w_AMPA_PYR0_PYR0_1s_NoDBC, delimiter=",")
np.savetxt("w_AMPA_PYR0_PYR2_1s_NoDBC.txt", w_AMPA_PYR0_PYR2_1s_NoDBC, delimiter=",")




#save weights after training
np.savetxt("w_AMPA_PYR0_PYR1_4s_NoDBC.txt", w_AMPA_PYR0_PYR1_4s_NoDBC, delimiter=",")
np.savetxt("w_AMPA_PYR0_PYR3_4s_NoDBC.txt", w_AMPA_PYR0_PYR3_4s_NoDBC, delimiter=",")
np.savetxt("w_AMPA_PYR0_PYR0_4s_NoDBC.txt", w_AMPA_PYR0_PYR0_4s_NoDBC, delimiter=",")
np.savetxt("w_AMPA_PYR0_PYR2_4s_NoDBC.txt", w_AMPA_PYR0_PYR2_4s_NoDBC, delimiter=",")

    
