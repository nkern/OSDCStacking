#!/bin/sh
#PBS -S /bin/sh
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/binstack_run_table4/bs_m0_run5/
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking
python caustic_mass_stack2D.py 6 6 5 25 0 5 4 0




