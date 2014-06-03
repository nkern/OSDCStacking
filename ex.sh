#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=10mb,walltime=1:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking/
#PBS -j oe
#

cd /glusterfs/users/caustics1/nkern/OSDCStacking/
echo $SHELL
python ex.py

