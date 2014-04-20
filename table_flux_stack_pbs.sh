#!/bin/sh
#PBS -S /bin/sh
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@/
#PBS -j oe
#PBS -t 0-@@job_array@@
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking

python caustic_mass_stack2D.py $PBS_ARRAYID @@clus_num@@ @@gal_num@@ @@line_num@@ @@method_num@@ @@cell_num@@ @@table_num@@ @@run_los@@


