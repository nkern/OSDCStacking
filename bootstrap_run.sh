#!/bin/bash
# This program iteratively submits qsub bootstrap_flux_stack_pbs.sh

echo -n "Are you sure you want to run the script bootstrap_run.sh? (y/n):"
read accept
if [ $accept != 'y' ];
	then echo 'Quitting...'
	exit
fi 

##################
## Begin Script ##
##################

## Initialize Configuration Arrays and Other Constants
cell_num=($(seq 1 49))				# Number of Cells
line_num=(2 5 10 15 25 50 100)			# Line of Sight Number 
gal_num=(5 10 15 25 50 100 150)			# Ngal number
halo_num=2100					# Number of Halos in Sample
method_num=0					# Ensemble Build Method
table_num=5					# Version of entire run table
write_stem="mm_m0_run"				# Stem of write_loc directory
data_loc="mass_mix/mm_0.05_run_table$table_num"	# Highest Directory for Data








