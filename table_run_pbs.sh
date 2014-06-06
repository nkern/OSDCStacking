#!/bin/sh
#PBS -S /bin/sh
#PBS -N @@job_name@@
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o @@root@@/OSDCStacking/@@data_loc@@/@@write_loc@@/
#PBS -j oe
#PBS -t 0-@@job_array@@
#

# Go to Working Directory
cd @@root@@/OSDCStacking/@@data_loc@@/@@write_loc@@

# Run Code with run_num dynamically defined
python @@root@@/OSDCStacking/@@filename@@.py $PBS_ARRAYID


