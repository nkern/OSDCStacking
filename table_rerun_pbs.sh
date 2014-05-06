#!/bin/sh
#PBS -S /bin/sh
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@/
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@

python caustic_mass_stack2D.py @@run_num@@ 


