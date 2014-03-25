#!/bin/sh
#PBS -S /bin/sh
#PBS -A christoq_flux
#PBS -l qos=flux
#PBS -q flux
#PBS -N table_analysis
#PBS -l nodes=1:ppn=1,pmem=30000mb,walltime=3:00:00
#PBS -V
#PBS -o /nfs/christoq_ls/nkern/Stacking/
#PBS -j oe
#

cd /nfs/christoq_ls/nkern/Stacking 
python flux_stack_recovery.py


