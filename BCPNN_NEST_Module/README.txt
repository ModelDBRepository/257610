Bernhard Kaplan, 2014


For NEST 2.2.2 [for NEST 2.4.2 some changes to the .cpp and .sli files need to be done, see below]

--First, set the environment variables
export PYTHONPATH=/home/bernhard/Downloads/nest-2.2.2-build/lib/python2.7/site-packages
export PATH=$PATH:/home/bernhard/Downloads/nest-2.2.2-build/bin
export LD_LIBRARY_PATH=/home/bernhard/workspace/BCPNN-Module/lib/nest


--Second, to add this module ('pt_module') to nest.Models():
rm -rf bootstrap-module-100725/
rm -rf build-module-100725/
./compile-module.sh module-100725/

--then, 
ipython
import nest
nest.Models()                     // 'bcpnn_synapse', 'iaf_cond_alpha_bias' should NOT show up
try:
    nest.Install('pt_module')
except:
    nest.Install('pt_module')
nest.Models()                     // now, 'bcpnn_synapse' is available

--use the new synapse model!


For Nest 2.4.2, some changes need to be done before compiling
- In your module's main .cpp file, in the function ::commandstring(),
change the returned string to just "(ml_module) run".

- In your module's .sli file, remove the line with "provide-component"
at the end.

Also removed the line with require-component at the end
in the ml:module.sli file to make it work.
