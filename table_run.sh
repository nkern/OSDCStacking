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
filename=millennium_stack			# Primary file to be run

### FLAGS ###
run_qsub=False					# If True perform qsub of scripts, if False create scripts but don't qsub

self_stack=False				# Run self-stack or bin-stack
scale_data=True					# Scale data by r200 if True
lightcone=False					# Run over Henriques Lightcone or Guo Data Cube?
write_data=True					# Write Data to Result directories if True
init_shiftgap=False				# Run a shiftgapper on individual LOS phase space before it gets stacked
shiftgapper=True				# Run a Shiftgapper Technique on Phase Space Before Caustic Technique
edge_int_remove=True				# Run Inflection Interloper Removal Technique
small_set=False					# 100 Halo Set or 2000 Halo Set
mass_mix=True					# Incorporate Mass Mixing Models?
bootstrap=False					# Perform a bootstrapping technique to estimate error in mass estimation?
new_halo_cent=True				# Use Updated Halo Centers instead of BCG Values
true_mems=False					# Run with only gals within r200?
run_los=False					# Run line of sight mass estimation or not
mirror=False					# Mirror Phase Space in Caustic Surface Estimation?
cent_offset=None				# Either 'r', 'v', 'full', or None.

### CONSTANTS ###
# Run Dependent
gal_num=(5 10 15 25 50 100 150)			# Ngal number
line_num=(2 5 10 15 25 50 100)			# Line of Sight Number 
method_num=0					# Ensemble Build Method Number
cell_num=($(seq 1 49))				# Number of Cells
table_num=1					# Table Re-Run Version  
job_name="MASS-MIX"				# PBS Job Name Stem
#halo_num=6000					# Total number of halos to work with
#halo_num=2100                                   # Total number of halos to work with
root="/glusterfs/users/caustics1/nkern"       # Base Directory

# Other Techniques
mm_est='richness'				# If mass mixing, observable to bin on, 'richness', 'veldisp'
edge_perc=0.1					# Percent of Top galaxies used in edge detection technique
center_scat=None				# If guessing halo center, fractional induced scatter into known center
avg_meth="'median'"				# If bin stacking, by which method do you average bin properties? (ex. median, mean)
bootstrap_num=None				# Highest directory marker for bootstrap data, ex. bootstrap1
bootstrap_rep=None				# Bootstrap repetition directory marker, ex. bootstrap1/rep1

# Location
write_stem="mm_m0_run"				# Stem of write_loc directory
data_loc="mass_mix/mm_run_table$table_num"	# Highest Directory for Data


## Go To Stacking Directory ##
cd $root/OSDCStacking

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

# Set Constants and Write README
echo ""
echo "...Setting Constants and Writing README file"
# Define Constants
if [ $lightcone == 'True' ]
then
        job_array=(19 19 19 9 9 9 9)
        halo_num=6000
else
        job_array=(9 9 9 9 6 6 6)
        if [ $small_set == 'True' ]
        then
                halo_num=100
        else
                halo_num=2100
        fi
fi

# Write README
sed -e "s:@@self_stack@@:$self_stack:g;s:@@scale_data@@:$scale_data:g;s:@@write_data@@:$write_data:g;s:@@init_shiftgap@@:$init_shiftgap:g;s:@@shiftgapper@@:$shiftgapper:g;s:@@edge_int_remove@@:$edge_int_remove:g;s:@@small_set@@:$small_set:g;s:@@mass_mix@@:$mass_mix:g;s:@@bootstrap@@:$bootstrap:g;s:@@new_halo_cent@@:$new_halo_cent:g;s:@@true_mems@@:$true_mems:g;s:@@run_los@@:$run_los:g;s:@@mirror@@:$mirror:g;s:@@cent_offset@@:$cent_offset:g;s:@@table_num@@:$table_num:g;s:@@data_loc@@:$data_loc:g;s:@@write_loc@@:$write_loc:g;s:@@mm_est@@:$mm_est:g;s:@@edge_perc@@:$edge_perc:g;s:@@center_scat@@:$center_scat:g;s:@@avg_meth@@:$avg_meth:g;s:@@root@@:$root:g;s:@@lightcone@@:$lightcone:g;" < readme_pbs.txt > $data_loc/readme.txt

# If Mass Mixing, create Halo_Arrays.pkl file
if [ -e $data_loc\/halo_arrays.pkl ]
then
	echo -n
else
	python millennium_massmix.py 0 mm	
fi

## Begin Nested Loops Through Table
for i in $(seq 0 6)
do
	for j in $(seq 0 6)
	do
		let "k=($i*7)+$j"


		echo '----------------------------------------------------------'
		echo -e "cell_num=${cell_num[$k]}\tgal_num=${gal_num[$i]}\tline_num=${line_num[$j]}"
		# Submit Job Array To PBS by feeding "table_run_pbs.sh" job parameters

		if [ $self_stack == 'True' ]
		then
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
		sed -e "s:@@self_stack@@:$self_stack:g;s:@@scale_data@@:$scale_data:g;s:@@write_data@@:$write_data:g;s:@@init_shiftgap@@:$init_shiftgap:g;s:@@shiftgapper@@:$shiftgapper:g;s:@@edge_int_remove@@:$edge_int_remove:g;s:@@small_set@@:$small_set:g;s:@@mass_mix@@:$mass_mix:g;s:@@bootstrap@@:$bootstrap:g;s:@@new_halo_cent@@:$new_halo_cent:g;s:@@true_mems@@:$true_mems:g;s:@@run_los@@:$run_los:g;s:@@mirror@@:$mirror:g;s:@@cent_offset@@:$cent_offset:g;s:@@ens_num@@:$ens_num:g;s:@@gal_num@@:$_gal_num:g;s:@@line_num@@:$_line_num:g;s:@@method_num@@:$method_num:g;s:@@cell_num@@:$_cell_num:g;s:@@table_num@@:$table_num:g;s:@@data_loc@@:$data_loc:g;s:@@write_loc@@:$write_loc:g;s:@@mm_est@@:$mm_est:g;s:@@edge_perc@@:$edge_perc:g;s:@@center_scat@@:$center_scat:g;s:@@avg_meth@@:$avg_meth:g;s:@@bootstrap_num@@:$bootstrap_num:g;s:@@bootstrap_rep@@:$bootstrap_rep:g;s:@@root@@:$root:g;s:@@lightcone@@:$lightcone:g;s:@@halo_num@@:$halo_num:g;" < caustic_params_pbs.py > $data_loc/$write_loc/caustic_params.py

		# Create run_pbs.sh file
		sed -e "s:@@filename@@:$filename:g;s:@@job_name@@:$job_name:g;s:@@write_loc@@:$write_loc:g;s:@@data_loc@@:$data_loc:g;s:@@job_array@@:$_job_array:g;s:@@root@@:$root:g" < table_run_pbs.sh > $data_loc/$write_loc/run_pbs.sh

		# Change run_los = True if line_num == 100 in caustic_params.py
#		let "a=${cell_num[$k]}%7"
#		if [ $a == 0 ]
#			then
#			echo -n
#			sed -e "s:@@self_stack@@:$self_stack:g;s:@@scale_data@@:$scale_data:g;s:@@write_data@@:$write_data:g;s:@@clean_ens@@:$clean_ens:g;s:@@small_set@@:$small_set:g;s:@@mass_mix@@:$mass_mix:g;s:@@bootstrap@@:$bootstrap:g;s:@@new_halo_cent@@:$new_halo_cent:g;s:@@true_mems@@:$true_mems:g;s:@@run_los@@:True:g;s:@@cent_offset@@:$cent_offset:g;s:@@ens_num@@:$ens_num:g;s:@@gal_num@@:$_gal_num:g;s:@@line_num@@:$_line_num:g;s:@@method_num@@:$method_num:g;s:@@cell_num@@:$_cell_num:g;s:@@table_num@@:$table_num:g;s:@@data_loc@@:$data_loc:g;s:@@write_loc@@:$write_loc:g;s:@@mass_scat@@:$mass_scat:g;s:@@center_scat@@:$center_scat:g;s:@@avg_meth@@:$avg_meth:g;s:@@bootstrap_num@@:$bootstrap_num:g;s:bootstrap_rep@@:$bootstrap_rep:g" < caustic_params_pbs.py > $data_loc/$write_loc/caustic_params.py
#		fi

		if [ $run_qsub == 'True' ]
			then
			echo "Submitting PBS Job"
			qsub $data_loc/$write_loc/run_pbs.sh 
		fi
	
		echo ""

		echo '----------------------------------------------------------'
	done
done




