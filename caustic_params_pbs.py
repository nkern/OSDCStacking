'''
This file contains run-dependent parameters used by caustic_mass_stack2D.py

file location : /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@
'''

## Run Dependent Constants ##

clus_num	= @@clus_num@@			# Number of Ensembles to build and solve for in this run
gal_num		= @@gal_num@@			# Number of galaxies taken per line of sight
line_num	= @@line_num@@			# Number of lines of sight to stack over
method_num	= @@method_num@@		# Ensemble Build Method Number
cell_num	= @@cell_num@@			# Cell Number ID corresponding to given gal_num & line_num geometry in a Run Table
table_num	= @@table_num@@			# Table Re-Run Version	
run_los 	= @@run_los@@			# If fed 8th arg value as True, run_los
mass_scat	= @@mass_scat@@			# If mass_mix = True, fractional scatter induced into table mass
center_scat	= @@center_scat@@		# If guessing halo center, fractional induced scatter into known center
avg_meth	= @@avg_meth@@			# If bin stacking, by which method do you average bin properties? (ex. median, mean)
cent_offset	= @@cent_offset@@		# Either 'r', 'v', 'full', or None.
bootstrap_num	= @@bootstrap_num@@		# Highest directory marker for bootstrap data, ex. bootstrap1
bootstrap_rep	= @@bootstrap_rep@@		# Bootstrap repetition directory marker, ex. bootstrap1/rep1

