#!/bin/sh
#PBS -S /bin/sh
#PBS -N table_analysis
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=3:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/binstack/bs_run_table5
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking/binstack/bs_run_table5

python /glusterfs/users/caustics1/nkern/OSDCStacking/millennium_recover.py 5


