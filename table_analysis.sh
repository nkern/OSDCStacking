#!/bin/sh
#PBS -S /bin/sh
#PBS -N table_analysis
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=1:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/binstack/
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking 
python flux_stack_recovery.py 1


