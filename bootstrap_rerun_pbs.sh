#!/bin/sh
#PBS -S /bin/sh
#PBS -N BOOTSTRAP
#PBS -l nodes=1:ppn=1,pmem=2000mb,walltime=2:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@/
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@

python /glusterfs/users/caustics1/nkern/OSDCStacking/caustic_mass_stack2D.py @@run_num@@ @@clus_num@@ @@gal_num@@ @@line_num@@ @@method_num@@ @@cell_num@@ @@bootstrap_num@@ @@run_los@@
