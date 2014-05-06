#!/bin/sh
#PBS -S /bin/sh
#PBS -N SELF-STACK
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@/
#PBS -j oe
#PBS -t 0-@@job_array@@
#

# Go to Working Directory
cd /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@

# Run Code with run_num dynamically defined
python /glusterfs/users/caustics1/nkern/OSDCStacking/caustic_mass_stack2D.py $PBS_ARRAYID


