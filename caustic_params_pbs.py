'''
This file contains run-dependent parameters used by millennium_stack.py

file location : /glusterfs/users/caustics1/nkern/OSDCStacking/@@data_loc@@/@@write_loc@@
'''
import time

### FLAGS ###
self_stack	= @@self_stack@@		# Run self-stack or bin-stack
scale_data	= @@scale_data@@		# Scale data by r200 if True
write_data 	= @@write_data@@		# Write Data to Result directories if True
init_clean	= @@init_clean@@			# Do an extra shiftgapper on ensemble before the lines of sight get stacked.
small_set	= @@small_set@@			# 100 Halo Set or 2000 Halo Set
mass_mix	= @@mass_mix@@			# Incorporate Mass Mixing Models?
bootstrap	= @@bootstrap@@			# Perform a bootstrapping technique to estimate error in mass estimation?
new_halo_cent	= @@new_halo_cent@@		# Use Updated Halo Centers instead of BCG Values
mirror		= @@mirror@@			# Mirror Phase Space in Caustic Surface Estimation?
true_mems	= @@true_mems@@			# Run with only gals within r200?
run_los         = @@run_los@@			# Run line of sight mass estimation or not
cent_offset     = @@cent_offset@@		# Either 'r', 'v', 'full', or None.

### CONSTANTS ###
# Run Dependent
ens_num		= @@ens_num@@			# Number of Ensembles to build and solve for IN THIS RUN
gal_num		= @@gal_num@@			# Number of galaxies taken per line of sight
line_num	= @@line_num@@			# Number of lines of sight to stack over
method_num	= @@method_num@@		# Ensemble Build Method Number
cell_num	= @@cell_num@@			# Cell Number ID corresponding to given gal_num & line_num geometry in a Run Table
table_num	= @@table_num@@			# Table Re-Run Version  
data_loc	= '@@data_loc@@'		# Alternative data_loc, either None or String
write_loc	= '@@write_loc@@'		# Alternative write_loc, either None or String

# Other Techniques
edge_perc       = @@edge_perc@@			# Percent of Top galaxies used in edge detection technique
mass_scat	= @@mass_scat@@			# If mass_mix = True, fractional scatter induced into table mass, feed as string, ex. "'0.25'"
center_scat	= @@center_scat@@		# If guessing halo center, fractional induced scatter into known center
avg_meth	= @@avg_meth@@			# If bin stacking, by which method do you average bin properties? (ex. median, mean)
bootstrap_num	= @@bootstrap_num@@		# Highest directory marker for bootstrap data, ex. bootstrap1
bootstrap_rep	= @@bootstrap_rep@@		# Bootstrap repetition directory marker, ex. bootstrap1/rep1

# Caustic Technique Dependent
q		= 50.0				# Scale of Gaussian Kernel Density Estimator
beta		= 0.2				# Velocity Anisotropy Beta parameter, if constant profile
fbeta		= 0.65				# fbeta value, see 'Diaferio 1999'
r_limit 	= 1.5				# Phase space radius Cut Scaled by R200
v_limit		= 3500.0			# Phase space velocity Cut in km/s

# Data Set
data_set	= 'Guo30_2'			# Data set to draw semi analytic data from
halo_num	= 2100				# Total number of halos to work with

# Physical
c               = 2.99792e5                     # speed of light in km/s
h               = 1.0                           # Hubble Constant, unitless
H0              = h*100.0                       # Hubble Constant, km s-1 Mpc-1

# Other
run_time	= time.asctime()                # Time when program was started
root 		= '/glusterfs/users/caustics1/nkern'  # Root for OSDC

