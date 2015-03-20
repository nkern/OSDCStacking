'''
Readme file detailing parameters of run for each full run_table
'''

### FLAGS ###
self_stack=@@self_stack@@			# Run self-stack or bin-stack
scale_data=@@scale_data@@			# Scale data by r200 if True
lightcone=@@lightcone@@				# Run over Henriques Lightcone or Guo Data Cube?
write_data=@@write_data@@			# Write Data to Result directories if True
init_shiftgap=@@init_shiftgap@@			# Run a shiftgapper on individual LOS phase space before it gets stacked
shiftgapper=@@shiftgapper@@			# Run a Shiftgapper Technique on Phase Space Before Caustic Technique
edge_int_remove=@@edge_int_remove@@		# Run Inflection Interloper Removal Technique
small_set=@@small_set@@				# 100 Halo Set or 2000 Halo Set
mass_mix=@@mass_mix@@				# Incorporate Mass Mixing Models?
bootstrap=@@bootstrap@@				# Perform a bootstrapping technique to estimate error in mass estimation?
new_halo_cent=@@new_halo_cent@@			# Use Updated Halo Centers instead of BCG Values
true_mems=@@true_mems@@				# Run with only gals within r200?
run_los=@@run_los@@				# Run line of sight mass estimation or not
mirror=@@mirror@@				# Mirror Phase Space in Caustic Surface Estimation?
cent_offset=@@cent_offset@@			# Either 'r', 'v', 'full', or None.

### CONSTANTS ###
# Run Dependent
gal_num=@@gal_num@@				# Ngal number
line_num=@@line_num@@				# Line of Sight Number 
jobs_array=@@job_array@@			# Number of PBS runs to feed for every run_table row (distinct line_of_sight values)
method_num=@@method_num@@			# Ensemble Build Method Number
cell_num=@@cell_num@@				# Number of Cells
table_num=@@table_num@@				# Table Re-Run Version  
job_name=@@job_name@@				# PBS Job Name Stem
halo_num=@@halo_num@@				# Total number of halos to work with
root='@@root@@'					# Base Directory

# Other Techniques
mm_est=@@mm_est@@				# If mass mxixing, observable to bin on, 'richness', 'veldisp'
edge_perc=@@edge_perc@@				# Percent of Top galaxies used in edge detection technique
center_scat=@@center_scat@@               # If guessing halo center, fractional induced scatter into known center
avg_meth=@@avg_meth@@				# If bin stacking, by which method do you average bin properties? (ex. median, mean)

# Location
data_loc='@@data_loc@@'				# Highest Directory for Data



