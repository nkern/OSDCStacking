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
cell_num=($(seq 1 49))						# Number of Cells
line_num=(2 5 10 15 25 50 100)					# Line of Sight Number 
gal_num=(5 10 15 25 50 100 150)					# Ngal number
halo_num=2100							# Number of Halos in Sample
method_num=0							# Ensemble Build Method
bootstrap_num=2							# Version of entire bootstrap directory
repetitions=($(seq 1 5))					# Array of rep directories to loop bootstrap over
write_stem="bo_m0_run"						# Stem of write_loc directory
data_loc="binstack/bootstrap$bootstrap_num"			# Highest Directory for Data
base_dir="/glusterfs/users/caustics1/nkern/OSDCStacking"	# Root directory

# Go to OSDCStacking Directory
cd $base_dir

## Double Check if directory chain exists
# Parent directories check
echo ""
echo "...Checking Parent Directories:"
if [ -d $data_loc ]
        then
        echo "Directory: $data_loc , exists..."
        echo -n "Do you want to continue?(y/n):"
        read accept
        if [ $accept != 'y' ]
                then echo 'Quitting...'
                exit
        fi
        else
        mkdir $data_loc
        echo "Created Directory: $data_loc"
        echo -n "Continue?(y/n):"
        read accept
        if [ $accept != 'y' ]
                then echo 'Quitting...'
                exit
        fi
fi

# Sub Directory Check
echo ""
echo "...Checking Subdirectories:"
for i in ${repetitions[*]}
do
	dir="$data_loc/rep$i"
	if [ -d $dir ]
		then
		echo -n
		else
		mkdir -v $dir
	fi
	for i in ${cell_num[*]}
	do
		sub_dir=$write_stem$i
		if [ -d $dir/$sub_dir ]
			then
			echo -n
			else
			mkdir -v $dir/$sub_dir
		fi
	done
done


# Begin Looping
for a in ${repetitions[*]}
do
	echo ""
	echo -e "#######################################"
	echo -e "## Working on Bootstrap Repetition "$a" ##"
	echo -e "#######################################"
	for i in $(seq 1 6)		# For now don't do 1st column and row of run table
	do
		for j in $(seq 1 6)
		do
			let "k=($i*7)+$j"
			echo '----------------------------------------------------------'
			echo -e "cell_num=${cell_num[$k]}\tgal_num=${gal_num[$i]}\tline_num=${line_num[$j]}"
			# Submit Job Array To FLUX by feeding "bootstrap_flux_stack_pbs.sh" job parameters
			# $1 : clus_num                 (first caustic_mass_stack2D.py positional parameter)
			# $2 : gal_num                  (second positional parameter)
			# $3 : line_num                 (third positional parameter)
			# $4 : method_num               (fourth positional parameter)
			# $5 : cell_num                 (fifth positional parameter)
			# $6 : job_array                (number of runs per job)
			# $7 : Submission Number        (appended to data_loc, ex: binstack_run_table3)
			# $8 : write_loc                (2nd level directory to write .pkl files in)

			# Define Constants to be used in job submission
			job_array=(13 13 13 13 13 13 20)
			let "clus_num=$halo_num/(${job_array[$j]}+1)/${line_num[$j]}"
			write_loc=$write_stem${cell_num[$k]}

			_clus_num="$clus_num"
			_gal_num="${gal_num[$i]}"
			_line_num="${line_num[$j]}"
			_method_num="$method_num"
			_cell_num="${cell_num[$k]}"
			_job_array="${job_array[$j]}"
			_data_loc="$data_loc/rep$a"
			_write_loc="$write_loc"

			# create bootstrap_params.py file in $data_loc/$write_loc directory
			sed -e "s:@@bootstrap_num@@:$bootstrap_num:g;s:@@bootstrap_rep@@:$a:g;s:@@file_location@@:$_data_loc/$_write_loc:g" < $base_dir/bootstrap_pbs_params.py > $base_dir/$_data_loc/$_write_loc/bootstrap_params.py

			# create script.sh file in $data_loc/$write_loc
			sed -e "s:@@write_loc@@:$_write_loc:g;s:@@data_loc@@:$_data_loc:g;s:@@job_array@@:$_job_array:g;s:@@clus_num@@:$_clus_num:g;s:@@gal_num@@:$_gal_num:g;s:@@line_num@@:$_line_num:g;s:@@method_num@@:$_method_num:g;s:@@cell_num@@:$_cell_num:g;s:@@bootstrap_num@@:$bootstrap_num:g;s:@@run_los@@:0:g" < $base_dir/bootstrap_pbs.sh > $base_dir/$_data_loc/$_write_loc/script.sh

			# Submit script.sh to PBS Scheduler
			qsub $_data_loc/$_write_loc/script.sh

		done
	done
done	













