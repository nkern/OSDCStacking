#!/bin/bash
# This program iteratively finds failed jobs and re-submits qsub bootstrap_pbs_rerun.sh

echo -n "Are you sure you want to run the script bootstrap_run_rerun.sh? (y/n):"
read accept
if [ $accept != 'y' ];
        then echo 'Quitting...'
        exit
fi

##################
## Begin Script ##
##################

## Initialize Configuration Arrays and Other Constants
cell_num=($(seq 1 49))                                          # Number of Cells
line_num=(2 5 10 15 25 50 100)                                  # Line of Sight Number 
gal_num=(5 10 15 25 50 100 150)                                 # Ngal number
clus_num=(75 30 15 10 6 3 1)					# Number of Ens Clusters done per instance
job_num=(14 14 14 14 14 14 21)					# Number of Jobs Submitted
halo_num=2100                                                   # Number of Halos in Sample
method_num=0                                                    # Ensemble Build Method
bootstrap_num=2                                                 # Version of entire bootstrap directory
repetitions=($(seq 1 5))                                       # Array of rep directories to loop bootstrap over
write_stem="bo_m0_run"                                          # Stem of write_loc directory
data_stem="binstack/bootstrap$bootstrap_num"                    # Highest Directory for Data
base_dir="/glusterfs/users/caustics1/nkern/OSDCStacking"        # Root directory


## Check Directory ##
echo "Loaded Directory is: $data_stem"
echo -n "Is this the desired directory? (y/n): "
read accept
if [ $accept != 'y' ];
	then echo 'Quitting...'
	exit
fi 


## Iterate over Repetitions, find failed jobs, re-submit them ##
for a in ${repetitions[*]}
do
	echo ""
	echo "Working on Repetition: rep$a"	
	echo "-------------------------------"
	data_loc="$data_stem/rep$a"

	## Find failed jobs ##
	DIRS=()
	NUMS=()
	# i loops over Ngal, j loops over LOS
	for i in {1..6}		# For now ignore first column and row of each run table
	do
		for j in {1..6}
		do
			let "k=$i*7+($j+1)"
			dir=$write_stem$k
			echo "Working on Directory: $dir..."
			m=0
			for n in $(seq 1 $(($halo_num/${line_num[$j]}/${job_num[$j]})) $(($halo_num/${line_num[$j]})))
			do	
				n=$((n-1))
				if [ -a $data_loc/$dir/Ensemble_$n\_Data.pkl ]
					then
					echo -n
				else
					DIRS+=($k)
					NUMS+=($m)	
				fi
				m=$((m+1))
			done
		done
	done

	## Check this worked
	echo "DIRS=${DIRS[*]}"
	echo "NUMS=${NUMS[*]}"
	echo -n "Does this seem right? (y/n):"
	read accept
	if [ $accept != 'y' ];
		then echo 'Quitting...'
		exit
	fi 

	## Begin PBS Job Submission
	echo "Beginning PBS job submission..."
	echo "-------------------------------------------"
	echo ""

	## Re create arrays
	line_num=(${line_num[*]} ${line_num[*]} ${line_num[*]} ${line_num[*]} ${line_num[*]} ${line_num[*]} ${line_num[*]})
	Job_Array=(13 13 13 13 13 13 20)
	job_array=(${Job_Array[*]} ${Job_Array[*]} ${Job_Array[*]} ${Job_Array[*]} ${Job_Array[*]} ${Job_Array[*]} ${Job_Array[*]})
	Gal_Num=(${gal_num[*]})
	gal_num=()
	for i in {0..6}
	do
		gal_num+=("${Gal_Num[$i]}")
		gal_num+=("${Gal_Num[$i]}")
		gal_num+=("${Gal_Num[$i]}")
		gal_num+=("${Gal_Num[$i]}")
		gal_num+=("${Gal_Num[$i]}")
		gal_num+=("${Gal_Num[$i]}")
		gal_num+=("${Gal_Num[$i]}")
	done

	## Go to Stacking Directory
	cd /glusterfs/users/caustics1/nkern/OSDCStacking

	m=0
	for k in ${DIRS[*]}
	do
		k=$((k-1))
		run_num=${NUMS[$m]}

		echo '----------------------------------------------------------'
		echo -e "cell_num=${cell_num[$k]}\tgal_num=${gal_num[$k]}\tline_num=${line_num[$k]}"
		# Submit Job Array To FLUX by feeding "mr_flux_stack_pbs.sh" job parameters
		# $1 : clus_num			(first caustic_mass_stack2D.py positional parameter)
		# $2 : gal_num			(second positional parameter)
		# $3 : line_num			(third positional parameter)
		# $4 : method_num		(fourth positional parameter)
		# $5 : cell_num			(fifth positional parameter)
		# $6 : job_array		(number of runs per job)
		# $7 : Submission Number	(appended to data_loc, ex: binstack_run_table3)
		# $8 : write_loc		(2nd level directory to write .pkl files in)

		# Define Constants
		let "clus_num=$halo_num/(${job_array[$k]}+1)/${line_num[$k]}"

		# Submit FLUX JOB for Ensembles
		_run_num="$run_num"
		_clus_num="$clus_num"
		_gal_num="${gal_num[$k]}"
		_line_num="${line_num[$k]}"
		_method_num="$method_num"
		_cell_num="${cell_num[$k]}"
		_data_loc="$data_loc"
		_write_loc="$write_stem${cell_num[$k]}"


		# Create bootstrap_params.py file in $data_loc/$write_loc directory
		sed -e "s:@@bootstrap_num@@:$bootstrap_num:g;s:@@bootstrap_rep@@:$a:g;s:@@file_location@@:$_data_loc/$_write_loc:g" < $base_dir/bootstrap_pbs_params.py > $base_dir/$_data_loc/$_write_loc/bootstrap_params.py

		# Create script_rerun.sh script
		sed -e "s:@@write_loc@@:$_write_loc:g;s:@@data_loc@@:$_data_loc:g;s:@@run_num@@:$_run_num:g;s:@@clus_num@@:$_clus_num:g;s:@@gal_num@@:$_gal_num:g;s:@@line_num@@:$_line_num:g;s:@@method_num@@:$_method_num:g;s:@@cell_num@@:$_cell_num:g;s:@@bootstrap_num@@:$bootstrap_num:g;s:@@run_los@@:0:g" < bootstrap_pbs_rerun.sh > $_data_loc/$_write_loc/script_rerun.sh

		# Submit Script to PBS via qsub
		qsub $_data_loc/$_write_loc/script_rerun.sh 

		m=$((m+1))

	done	
done












