## henriques_stack.py
"""
This program uses data from the henriques light cone from the millennium simulation to stack galaxy clusters and produce an 'ensemble cluster'.
We estimate the mass of this ensemble using the Caustic Technique.
"""

###################### BEGIN STANDARD PREAMBLE #####################
print '...importing modules'
from caustic_stack import *
from stack_class import *

sys.path.insert(0,os.getcwd())
__import__('caustic_params')
from caustic_params import *
print "Loaded caustic_params from",sys.modules['caustic_params']

## Fed Positional Parameter
run_num			= int(sys.argv[1])		# run_num th iteration of the whole job array in PBS script
	
if self_stack == True:
	""" Self Stacking """
        stack_range = np.arange(run_num*ens_num,run_num*ens_num+ens_num)	# Range of halos, each to be stacked individually
	if data_loc == None:
		data_loc = 'selfstack/ss_run_table'+str(table_num)		# Parent Directory where write_loc directories live
	if write_loc == None:
		write_loc = 'ss_m'+str(method_num)+'_run'+str(cell_num)		# Self Stack data-write location
else:
	""" Bin Stacking """
        stack_range = np.arange(run_num*ens_num*line_num,run_num*ens_num*line_num+ens_num*line_num)
	if data_loc == None:
		data_loc = 'binstack/bs_run_table'+str(table_num)
	if write_loc == None:
		write_loc = 'bs_m'+str(method_num)+'_run'+str(cell_num)			# Bin Stack data-write location
	if cent_offset == 'r':
		data_loc = 'binstack/roff_run_table'+str(table_num)
	elif cent_offset == 'v':
		data_loc = 'binstack/voff_run_table'+str(table_num)
	if mass_mix == True:							# Change write_loc if mass mixing
		data_loc = 'mass_mix/mm_'+str(mass_scat)+'_run_table'+str(table_num)
		write_loc = 'mm_m'+str(method_num)+'_run'+str(cell_num)
		if cent_offset == 'r':
			data_loc = 'mass_mix/roff_run_table'+str(table_num)
		elif cent_offset == 'v':
			data_loc = 'mass_mix/voff_run_table'+str(table_num)
if bootstrap == True:
	write_loc = 'bo_m'+str(method_num)+'_run'+str(cell_num)
	data_loc = 'binstack/bootstrap'+str(bootstrap_num)+'/rep'+str(bootstrap_rep)

## Make dictionary for loaded constants, doesn't matter if some don't exist
keys = ['c','h','H0','q','beta','fbeta','r_limit','v_limit','data_set','halo_num','gal_num','line_num','method_num','write_loc','data_loc','root','self_stack','scale_data','write_data','run_time','init_clean','small_set','run_los','bootstrap','run_num','ens_num','cell_num','stack_range','mass_mix','mass_scat','bootstrap_num','bootstrap_rep','avg_meth','cent_offset','center_scat','new_halo_cent','true_mems','mirror','edge_perc','lightcone']
varib = ez.create(keys,locals())
varib.update({'_name_':'varib'})

## INITIALIZATION ##
S = Stack(varib)
U = Universal(varib)
M = Millennium(varib)

U.print_separation('## Running millennium_stack.py')
names = ['run_time','','run_num','gal_num','line_num','cell_num','ens_num','halo_num','method_num','avg_meth','','self_stack','lightcone','mass_mix','write_data','scale_data','run_los','mirror','new_halo_cent','cent_offset','true_mems','init_clean','bootstrap','','mass_scat','center_scat','bootstrap_num','bootstrap_rep','','data_loc','write_loc','root']
U.print_varibs(names,varib)

####################### END STANDARD PREAMBLE ##############################

## Load Halo Data
U.print_separation('# ...Loading Halos',type=1)
HaloID,HaloData = M.load_halos()

# Sort Halos by A Priori Known Descending Mass (Mass Critical 200)
HaloID,HaloData = M.sort_halos(HaloID,HaloData)





