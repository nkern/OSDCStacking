#!/bin/sh
#PBS -S /bin/sh
#PBS -N table_analysis
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/selfstack/ss_run_table1
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking/selfstack/ss_run_table1

python /glusterfs/users/caustics1/nkern/OSDCStacking/caustic_stack_recovery.py 1


