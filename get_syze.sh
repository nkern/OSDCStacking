#!/bin/sh
#PBS -S /bin/sh
#PBS -N SYZE
#PBS -l nodes=1:ppn=1,pmem=1000mb,walltime=1:00:00
#PBS -V
#PBS -o /glusterfs/users/caustics1/nkern/OSDCStacking
#PBS -j oe
#

TZ=America/Detroit date

cd /glusterfs/users/caustics1/nkern
du -hcs *

