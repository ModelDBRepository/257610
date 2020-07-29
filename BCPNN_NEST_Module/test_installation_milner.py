import nest
print nest.Models()                     # 'bcpnn_synapse' should NOT show up

# ON LINDGREN: instead of nest.Install('pt_module'), do:
#if (not 'bcpnn_synapse' in nest.Models('synapses')):
#    nest.sr('(/home/bernhard/Downloads/nest/nest-2.2.1-build/share/nest/sli) addpath')
#    nest.Install('/home/bernhard/workspace/BCPNN-module/build-module-100725/pt_module')

on_milner = True
if (not 'bcpnn_synapse' in nest.Models('synapses')):
    if on_milner:
	nest.sr('(/cfs/milner/scratch/b/bkaplan/BCPNN-Module/share/nest/sli) addpath')
	nest.Install('/cfs/milner/scratch/b/bkaplan/BCPNN-Module/lib/nest/pt_module')
    else:
	nest.Install('pt_module')

nest.Models()                     # now, 'bcpnn_synapse' is available
nest.GetDefaults('bcpnn_synapse')
