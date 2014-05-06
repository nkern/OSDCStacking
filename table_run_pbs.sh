#!/bin/sh
#PBS -S /bin/sh
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@/
#PBS -j oe
#PBS -t 0-@@job_array@@
#

# Go to Working Directory
cd /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@

# Replace @@run_num@@ in caustic_params.py with $PBS_ARRAYID
sed -i "s:@@run_num@@:$PBS_ARRAYID:g" caustic_params.py

python caustic_mass_stack2D.py 


