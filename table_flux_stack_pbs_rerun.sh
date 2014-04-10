#!/bin/sh
#PBS -S /bin/sh
#PBS -N BIN-STACK
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/mass_mix/mm_0.50_run_table@@table_num@@/@@write_loc@@/
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking

python caustic_mass_stack2D.py @@run_num@@ @@clus_num@@ @@gal_num@@ @@line_num@@ @@method_num@@ @@cell_num@@ @@table_num@@ @@run_los@@


