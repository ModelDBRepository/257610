#!/bin/bash -l
# The name of the script is myjob
#SBATCH -J bcpnn_test_job

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 0:10:00

# Number of MPI tasks.
# always ask for complete nodes (i.e. mppwidth should normally
# be a multiple of 20)
#SBATCH -n 20

#SBATCH -e error_file.e
#SBATCH -o output_file.o


#load the nest module
module swap PrgEnv-cray PrgEnv-gnu
module add nest

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cfs/milner/scratch/b/bkaplan/BCPNN-Module/build-module-100725
export PYTHONPATH=/pdc/vol/nest/2.2.2/lib/python2.7/site-packages:/pdc/vol/python/2.7.6-gnu/lib/python2.7/site-packages

echo "Starting test 1 at `date`"
# Run and write the output into my_output_file
aprun -n 20 python test_installation_milner.py > delme_test_installation 2>&1
echo "Stopping test 1 at `date`"

echo "Starting test 2 at `date`"
aprun -n 20 python test_installation_milner.py > delme_test_connect 2>&1
echo "Stopping test 2 at `date`"




