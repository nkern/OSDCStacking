#!/bin/bash
# This program locates all flux_stack_pbs.sh runs that failed and re-runs them.

######################
#### Begin Script ####
######################

## Initialize Configuration Arrays and Other Constants
filename=millennium_bootstrap			# Primary file to be run

### FLAGS ###
self_stack=False				# Run self-stack or bin-stack
scale_data=True					# Scale data by r200 if True
write_data=True					# Write Data to Result directories if True
init_clean=False				# Do an extra shiftgapper on ensemble before the lines of sight get stacked.
small_set=False					# 100 Halo Set or 2000 Halo Set
mass_mix=False					# Incorporate Mass Mixing Models?
bootstrap=True					# Perform a bootstrapping technique to estimate error in mass estimation?
new_halo_cent=True				# Use Updated Halo Centers instead of BCG Values
true_mems=False					# Run with only gals within r200?
run_los=False					# Run line of sight mass estimation or not
mirror=True                                     # Mirror Phase Space in Caustic Surface Estimation?
cent_offset=None				# Either 'r', 'v', 'full', or None.

### CONSTANTS ###
# Run Dependent
gal_num=(5 10 15 25 50 100 150)			# Ngal number
line_num=(2 5 10 15 25 50 100)			# Line of Sight Number 
method_num=0					# Ensemble Build Method Number
cell_num=($(seq 1 49))				# Number of Cells
table_num=1					# Table Re-Run Version  
job_name="BOOTSTRAP"				# PBS Job Name Stem
job_num=(14 14 14 14 14 14 21)			# Number of Jobs Submitted
halo_num=2100					

# Other Techniques
edge_perc=0.1					# Percent of Top galaxies used in edge detection technique
mass_scat=None					# If mass_mix = True, fractional scatter induced into table mass, feed as string, ex. "'0.25'"
center_scat=None				# If guessing halo center, fractional induced scatter into known center
avg_meth="'median'"				# If bin stacking, by which method do you average bin properties? (ex. median, mean)
bootstrap_num=1					# Highest directory marker for bootstrap data, ex. bootstrap1
bootstrap_rep=$1				# Bootstrap repetition directory marker, ex. bootstrap1/rep1

# Location
write_stem="bo_m0_run"  		                       	    # Stem of write_loc directory
data_loc="binstack/bootstrap$bootstrap_num/rep$bootstrap_rep"     # Highest Directory for Data
base_dir="/glusterfs/users/caustics1/nkern" 			    # Base Directory

## For Bootstrapping Only
if [ $bootstrap == 'True' ]
then
        cell_select=(9 13 23 27)        # Cells to run bootstrap over
fi

## Go To Stacking Directory ##
cd $base_dir/OSDCStacking

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

		### Bootstrap ###
		pass=0  
		for z in ${cell_select[*]}
		do      
			if [ $z == $k ]
			then    
				pass=$((pass+1))
			fi      
		done    
		if [ $pass == 0 ]
		then    
			continue
		fi      
		#################

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
cd $base_dir/OSDCStacking

m=0
for k in ${DIRS[*]}
do
	k=$((k-1))
	run_num=${NUMS[$m]}

	echo '----------------------------------------------------------'
	echo -e "cell_num=${cell_num[$k]}\tgal_num=${gal_num[$k]}\tline_num=${line_num[$k]}"
	# Submit Job Array To PBS by feeding "table_rerun_pbs.sh" job parameters
	# $1 : ens_num			(first caustic_mass_stack2D.py positional parameter)
	# $2 : gal_num			(second positional parameter)
	# $3 : line_num			(third positional parameter)
	# $4 : method_num		(fourth positional parameter)
	# $5 : cell_num			(fifth positional parameter)
	# $6 : job_array		(number of runs per job)
	# $7 : Submission Number	(appended to data_loc, ex: binstack_run_table3)
	# $8 : write_loc		(2nd level directory to write .pkl files in)

	# Define Constants
	let "ens_num=$halo_num/(${job_array[$k]}+1)/${line_num[$k]}"

	# Define constants to be replaced in pbs script
	_run_num="$run_num"
	_gal_num="${gal_num[$k]}"
	_line_num="${line_num[$k]}"
	_cell_num="${cell_num[$k]}"
	write_loc="$write_stem${cell_num[$k]}"

	# Create script_rerun.sh file
	sed -e "s:@@filename@@:$filename:g;s:@@job_name@@:$job_name:g;s:@@write_loc@@:$write_loc:g;s:@@data_loc@@:$data_loc:g;s:@@run_num@@:$_run_num:g" < table_rerun_pbs.sh > $data_loc/$write_loc/script_rerun.sh

	# Submit Script to PBS via qsub	
	echo "Submitting Job to PBS"
	qsub $data_loc/$write_loc/script_rerun.sh 
	echo ""

	m=$((m+1))
done


