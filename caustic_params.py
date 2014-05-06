'''
This file contains run-dependent parameters used by caustic_mass_stack2D.py

file location : /glusterfs/users/caustics1/nkern/OSDCStacking 
'''

## Run Dependent Constants ##

clus_num	= 3				# Number of Ensembles to build and solve for in this run
gal_num		= 25				# Number of galaxies taken per line of sight
line_num	= 5				# Number of lines of sight to stack over
method_num	= 1				# Ensemble Build Method Number
cell_num	= 0				# Cell Number ID corresponding to given gal_num & line_num geometry in a Run Table
table_num	= 1				# Table Re-Run Version	
run_los 	= False				# If fed 8th arg value as True, run_los
mass_scat	= '0.05'			# If mass_mix = True, fractional scatter induced into table mass
avg_meth	= 'median'			# If bin stacking, by which method do you average bin properties? (ex. median, mean)
bootstrap_num	= 1				# Highest directory marker for bootstrap data, ex. bootstrap1
bootstrap_rep	= 1				# Bootstrap repetition directory marker, ex. bootstrap1/rep1


