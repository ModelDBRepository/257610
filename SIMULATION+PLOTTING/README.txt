Simulation code accompanying the manuscript:
"Introducing double bouquet cells into a modular cortical associative memory model"
By Nikolaos Chrysanthidis, Florian Fiebig, Anders Lansner
Manuscript submitted to Springer, Journal of Computational Neuroscience (JCNS)


We use NEST (Neural Simulation Tool) version 2.4.2 along with Python 2.7
To install NEST 2.4.2 please follow the instructions at the following link: NEST 2.4.2 https://nest-simulator.readthedocs.io/en/latest/installation/oldvers_install.html

In particular and as far as the simulation code is concerned, 
In the DBCmodel.py, the parameters we used to achieve satisfactory electrophysiological fidelity are included.
The simulations aim at reproducing spike patterns under sweeps of increasing  suprathreshold current 
steps (10 pA each) and other reported activity. The range of the stimulation input current is on the same level with the one reported in the paper below.
    
The spike patterns produced (figure DBC_ActivityPatterns) can be directly compared with the findings of fig.4B appeared in Cluster 
analysis–Based Physiological Classification and Morphological Properties of Inhibitory Neurons in Layers 2–3 of Monkey Dorsolateral Prefrontal Cortex (Krimer et al., 2005).

After installing BCPNN module (see. BCPNN_NEST_Module), in order for the figures related to the publication (Introducing double bouquet cells into a modular cortical associative
memory model) to be produced, the following scripts should be compiled: 

mainNewModel.py, plots.py in Figure 1C+3A folder --> Fig.1C (Membrane voltage of a stimulated DBC), Fig.3A (Spike raster of neurons in HC0). For plotting purposes (Fig.1C) the resolution is set to 0.001 for high spikes' amplitude quality.

multipleRoundsNewModel.py, multipleRoundsPreviousModel.py, WeightsDistribution.py in Figure 2 folder --> Fig.2A, Fig.2B, Fig.2C, Fig.2D (Weights distribution, multi-trial averages). 

mainNewModel.py, mainPreviousModel.py, I_GABA.py in Figure 3B folder --> Fig.3B (Functionality verification - Population averaged total inhibitory input current received by pyramidal cells in MC0 in both architectures)

