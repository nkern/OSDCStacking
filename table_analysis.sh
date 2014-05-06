#!/bin/sh
#PBS -S /bin/sh
#PBS -N table_analysis
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=45:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/mass_mix/mm_0.05_run_table1
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking 
python flux_stack_recovery.py 5


