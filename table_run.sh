#!/bin/bash
# This program iteratively submits qsub table_run_pbs.sh to PBS.

echo -n "Are you sure you want to run the script table_run.sh? (y/n):"
read accept
if [ $accept != 'y' ];
	then echo 'Quitting...'
	exit
fi 

##################
## Begin Script ##
##################

## Initialize Configuration Arrays and Other Constants
### FLAGS ###
self_stack=False				# Run self-stack or bin-stack
scale_data=True					# Scale data by r200 if True
write_data=True					# Write Data to Result directories if True
clean_ens=False					# Do an extra shiftgapper on ensemble before the lines of sight get stacked.
small_set=False					# 100 Halo Set or 2000 Halo Set
mass_mix=False					# Incorporate Mass Mixing Models?
bootstrap=False					# Perform a bootstrapping technique to estimate error in mass estimation?
new_halo_cent=False				# Use Updated Halo Centers instead of BCG Values
true_mems=False					# Run with only gals within r200?
run_los=False					# Run line of sight mass estimation or not
cent_offset=None				# Either 'r', 'v', 'full', or None.

### CONSTANTS ###
# Run Dependent
gal_num=(5 10 15 25 50 100 150)			# Ngal number
line_num=(2 5 10 15 25 50 100)			# Line of Sight Number 
method_num=0					# Ensemble Build Method Number
cell_num=($(seq 1 49))				# Number of Cells
table_num=1					# Table Re-Run Version  
data_loc="binstack/bs_run_table$table_num"	# Highest Directory for Data
base_dir="/glusterfs/users/caustics1/nkern"	# Base Directory
job_name="BIN-STACK"				# PBS Job Name Stem
write_stem="bs_m0_run"				# Stem of write_loc directory

# Other Techniques
mass_scat=None					# If mass_mix = True, fractional scatter induced into table mass, feed as string, ex. "'0.25'"
center_scat=None				# If guessing halo center, fractional induced scatter into known center
avg_meth="'median'"				# If bin stacking, by which method do you average bin properties? (ex. median, mean)
bootstrap_num=None				# Highest directory marker for bootstrap data, ex. bootstrap1
bootstrap_rep=None				# Bootstrap repetition directory marker, ex. bootstrap1/rep1


## Go To Stacking Directory ##
cd $base_dir/OSDCStacking

## Double Check if Directories Exist
# Parent Directory Check
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

# Sub Directory Check (including LOS)
echo ""
echo "...Checking subdirectories"
for i in ${cell_num[*]}
do
	dir=$write_stem$i
	if [ -d $data_loc\/$dir ]
		then
		echo -n
		else
		mkdir -v $data_loc/$dir	
	fi
done


## Begin Nested Loops Through Table
for i in $(seq 0 6)
do
	for j in $(seq 0 6)
	do
		let "k=($i*7)+$j"
		echo '----------------------------------------------------------'
		echo -e "cell_num=${cell_num[$k]}\tgal_num=${gal_num[$i]}\tline_num=${line_num[$j]}"
		# Submit Job Array To PBS by feeding "table_run_pbs.sh" job parameters
		# $1 : ens_num			(first caustic_mass_stack2D.py positional parameter)
		# $2 : gal_num			(second positional parameter)
		# $3 : line_num			(third positional parameter)
		# $4 : method_num		(fourth positional parameter)
		# $5 : cell_num			(fifth positional parameter)
		# $6 : job_array		(number of runs per job)
		# $7 : Submission Number	(appended to data_loc, ex: binstack_run_table3)
		# $8 : write_loc		(2nd level directory to write .pkl files in)

		# Define Constants
		job_array=(13 13 13 13 13 13 20)
		if [ $self_stack == 'True' ]
			then
			job_array=(13 13 13 13 13 13 20) #349) Not doing LOS runs for now
			let "ens_num=$halo_num/(${job_array[$j]}+1)"
			else
			let "ens_num=$halo_num/(${job_array[$j]}+1)/${line_num[$j]}"
		fi
		write_loc=$write_stem${cell_num[$k]}

		# Define constants to be replaced in pbs script 
		_gal_num="${gal_num[$i]}"
		_line_num="${line_num[$j]}"
		_cell_num="${cell_num[$k]}"
		_job_array="${job_array[$j]}"

		# Create caustic_params.py file in working directory
		sed -e "s:@@self_stack@@:$self_stack:g;s:@@scale_data@@:$scale_data:g;s:@@write_data@@:$write_data:g;s:@@clean_ens@@:$clean_ens:g;s:@@small_set@@:$small_set:g;s:@@mass_mix@@:$mass_mix:g;s:@@bootstrap@@:$bootstrap:g;s:@@new_halo_cent@@:$new_halo_cent:g;s:@@true_mems@@:$true_mems:g;s:@@run_los@@:$run_los:g;s:@@cent_offset@@:$cent_offset:g;s:@@ens_num@@:$ens_num:g;s:@@gal_num@@:$_gal_num:g;s:@@line_num@@:$_line_num:g;s:@@method_num@@:$method_num:g;s:@@cell_num@@:$_cell_num:g;s:@@table_num@@:$table_num:g;s:@@data_loc@@:$data_loc:g;s:@@write_loc@@:$write_loc:g;s:@@mass_scat@@:$mass_scat:g;s:@@center_scat@@:$center_scat:g;s:@@avg_meth@@:$avg_meth:g;s:@@bootstrap_num@@:$bootstrap_num:g;s:bootstrap_rep@@:$bootstrap_rep:g" < caustic_params_pbs.py > $data_loc/$write_loc/caustic_params.py

		# Create script.sh file
		sed -e "s:@@job_name@@:$job_name:g;s:@@write_loc@@:$write_loc:g;s:@@data_loc@@:$data_loc:g;s:@@job_array@@:$_job_array:g" < table_run_pbs.sh > $data_loc/$write_loc/script.sh

		# Change run_los = True if line_num == 100 in caustic_params.py
		let "a=${cell_num[$k]}%7"
		if [ $a == 0 ]
			then
			echo -n
		#	sed -e "s:@@self_stack@@:$self_stack:g;s:@@scale_data@@:$scale_data:g;s:@@write_data@@:$write_data:g;s:@@clean_ens@@:$clean_ens:g;s:@@small_set@@:$small_set:g;s:@@mass_mix@@:$mass_mix:g;s:@@bootstrap@@:$bootstrap:g;s:@@new_halo_cent@@:$new_halo_cent:g;s:@@true_mems@@:$true_mems:g;s:@@run_los@@:True:g;s:@@cent_offset@@:$cent_offset:g;s:@@ens_num@@:$ens_num:g;s:@@gal_num@@:$_gal_num:g;s:@@line_num@@:$_line_num:g;s:@@method_num@@:$method_num:g;s:@@cell_num@@:$_cell_num:g;s:@@table_num@@:$table_num:g;s:@@data_loc@@:$data_loc:g;s:@@write_loc@@:$write_loc:g;s:@@mass_scat@@:$mass_scat:g;s:@@center_scat@@:$center_scat:g;s:@@avg_meth@@:$avg_meth:g;s:@@bootstrap_num@@:$bootstrap_num:g;s:bootstrap_rep@@:$bootstrap_rep:g" < caustic_params_pbs.py > $data_loc/$write_loc/caustic_params.py
		fi

		# Submit Script to PBS via qsub
#		qsub $data_loc/$write_loc/script.sh 
		echo ""

		echo '----------------------------------------------------------'
	done
done




