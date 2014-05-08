#!/bin/bash
# This program locates all flux_stack_pbs.sh runs that failed and re-runs them.

######################
#### Begin Script ####
######################

## Initialize Configuration Arrays and Other Constants
self_stack=True                                                 # Self Stacking or Bin Stacking?
mass_mix=False                                                  # If Bin Stacking, Mass Mixing or not?
center_guess=None						# Vary halo center? 'r', 'v', or None

cell_num=($(seq 1 49))						# Number of Cells
line_num=(2 5 10 15 25 50 100)					# Line of Sight Number 
gal_num=(5 10 15 25 50 100 150)					# Ngal number
clus_num=(75 30 15 10 6 3 1)					# Number of Ens Clusters done per instance
job_num=(14 14 14 14 14 14 21)					# Number of Jobs Submitted
halo_num=2100							# Number of Halos in Sample
method_num=1							# Ensemble Build Method
table_num=1							# Version of entire run table
mass_scat=None							# If mass mixing, induced fractional scatter
data_loc="selfstack/ss_run_table"$table_num			# Highest Directory for Data
write_stem="ss_m1_run"						# Stem of write_loc directory
data_loc="selfstack/ss_run_table$table_num"			# Highest Directory for Data
base_dir="/glusterfs/users/caustics1/nkern/OSDCStacking"	# Base Directory


## Go To Stacking Directory
cd $base_dir 

## Check Directory ##
echo "Loaded Directory is: $data_loc"
echo -n "Is this the desired directory? (y/n): "
read accept
if [ $accept != 'y' ];
	then echo 'Quitting...'
	exit
fi 


## Find failed jobs ##
DIRS=()
NUMS=()
# i loops over Ngal, j loops over LOS
for i in {0..6}
do
	for j in {0..6}
	do
		let "k=$i*7+($j+1)"
		dir=$write_stem$k
		echo "Working on Directory: $dir..."
		m=0
		if [ $self_stack == 'True' ]
		then 
			iter=$(seq 1 $(($halo_num/${job_num[$j]})) $(($halo_num)))
		else 
			iter=$(seq 1 $(($halo_num/${line_num[$j]}/${job_num[$j]})) $(($halo_num/${line_num[$j]})))
		fi
		for n in $iter
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
echo "Beginning FLUX job submission..."
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
cd $base_dir

m=0
for k in ${DIRS[*]}
do
	k=$((k-1))
	run_num=${NUMS[$m]}

	echo '----------------------------------------------------------'
	echo -e "cell_num=${cell_num[$k]}\tgal_num=${gal_num[$k]}\tline_num=${line_num[$k]}"
	# Submit Job Array To PBS by feeding "table_rerun_pbs.sh" job parameters
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

	# Define constants to be replaced in pbs script
	_run_num="$run_num"
	_clus_num="$clus_num"
	_gal_num="${gal_num[$k]}"
	_line_num="${line_num[$k]}"
	_method_num="$method_num"
	_cell_num="${cell_num[$k]}"
	_data_loc="$data_loc"
	_write_loc="$write_stem${cell_num[$k]}"

	# Create script_rerun.sh file
	sed -e "s:@@write_loc@@:$_write_loc:g;s:@@data_loc@@:$_data_loc:g;s:@@run_num@@:$_run_num:g" < table_rerun_pbs.sh > $_data_loc/$_write_loc/script_rerun.sh

	# Submit Script to PBS via qsub
	qsub $_data_loc/$_write_loc/script_rerun.sh 
	echo ""

	m=$((m+1))
done


