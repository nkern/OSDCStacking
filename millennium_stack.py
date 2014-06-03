## millennium_stack.py
"""
This program uses the Caustic Technique to estimate the masses of galaxy clusters,
after having applied a stacking technique to create ensemble clusters.
This script does strictly self stacking or bin stacking.
"""

###################### BEGIN STANDARD PREAMBLE #####################
print '...importing modules'
from caustic_stack import *
from analysis import *
from stack_class import *

sys.path.insert(0,os.getcwd())
__import__('caustic_params')
from caustic_params import *
print "Loaded caustic_params from",sys.modules['caustic_params']

## Fed Positional Parameter ##
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
keys = ['c','h','H0','q','beta','fbeta','r_limit','v_limit','data_set','halo_num','gal_num','line_num','method_num','write_loc','data_loc','root','self_stack','scale_data','write_data','run_time','init_clean','small_set','run_los','bootstrap','run_num','ens_num','cell_num','stack_range','mass_mix','mass_scat','bootstrap_num','bootstrap_rep','avg_meth','cent_offset','center_scat','new_halo_cent','true_mems','mirror']
varib = ez.create(keys,locals())
varib.update({'_name_':'varib'})

## INITIALIZATION ##
S = Stack(varib)
U = Universal(varib)
A = Analysis(varib)
M = Millennium(varib)

U.print_separation('## Running caustic_mass_stack.py')
U.print_varibs(varib)

## Load Halo Data
U.print_separation('# ...Loading Halos',type=1)
HaloID,HaloData = M.load_halos()

# Sort Halos by A Priori Known Descending Mass (Mass Critical 200)
HaloID,HaloData = M.sort_halos(HaloID,HaloData)

####################### END STANDARD PREAMBLE ##############################

# Unpack HaloData array into local namespace
M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData
HaloData = M_crit200,R_crit200,HVD
Halo_P,Halo_V = np.vstack([HPX,HPY,HPZ]).T,np.vstack([HVX,HVY,HVZ]).T

# Initialize Multi-Ensemble Array to hold resultant data
STACK_DATA = []

U.print_separation('# ...Starting Ensemble Loop',type=1)
#  j: index for each ensemble, w.r.t. the FOR loop
#  k: index for each final ensemble, w.r.t. total number of final ensembles
#  l: index for line of sight, w.r.t. total lines of sight of run

for j in range(ens_num):

	# Update S.j
	S.j = j

	# Total Halo Index
	k = run_num*ens_num + j

	# Define Container to Hold Phase Space Data
	PS = Data()

	# Iterate through lines of sight
	U.print_separation('## Working on Ensemble '+str(j),type=1)
	ens_gal_count = 0
	for l in range(line_num):
		print '...Loading galaxy data and projecting line of sight #'+str(l)

		# Load galaxy data, project it, then append to PS
		if self_stack == True:
			M.load_project_append(HaloID[stack_range][j],M_crit200[stack_range][j],R_crit200[stack_range][j],HVD[stack_range][j],Z[stack_range][j],Halo_P[stack_range][j],Halo_V[stack_range][j],PS)

		else:
			M.load_project_append(HaloID[stack_range][j*line_num:(j+1)*line_num][l],M_crit200[stack_range][j*line_num:(j+1)*line_num][l],R_crit200[stack_range][j*line_num:(j+1)*line_num][l],HVD[stack_range][j*line_num:(j+1)*line_num][l],Z[stack_range][j*line_num:(j+1)*line_num][l],Halo_P[stack_range][j*line_num:(j+1)*line_num][l],Halo_V[stack_range][j*line_num:(j+1)*line_num][l],PS)


	PS.to_array(['Rdata','Vdata','HaloID','M200','R200','HVD','G_mags','R_Mags','I_Mags'])

	# Build Ensemble and Run Caustic Technique
	stack_data = S.caustic_stack(PS.Rdata,PS.Vdata,PS.HaloID,np.vstack([PS.M200,PS.R200,PS.HVD]),line_num,feed_mags=True,G_Mags=PS.G_Mags,R_Mags=PS.R_Mags,I_Mags=PS.I_Mags)

	# Append other Arrays
	extra = {'pro_pos':PS.pro_pos,'_name_':'stack_data'}
	stack_data.update(extra)

	# Append to STACK_DATA
	STACK_DATA.append(stack_data)


# Finished Loop
U.print_separation('#...Finished Ensemble Loop',type=2)

# Create run_dict
names = ['HaloID']
run_dict = ez.create(names,locals())
run_dict.update({'_name_':'run_dict'})

### Save Data into .pkl Files ###
if write_data == True:
	U.print_separation('##...Starting Data Write',type=2)
	for m in range(ens_num):
		n = run_num*ens_num + m
		U.print_separation("Writing Data for Ensemble #"+str(n),type=2)
		print 'Writing File: '+root+'/OSDCStacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(n)+'_Data.pkl'
		pkl_file = open(root+'/OSDCStacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(n)+'_Data.pkl','wb')
		output = pkl.Pickler(pkl_file)
		output.dump(STACK_DATA[m])
		output.dump(varib)
		output.dump(run_dict)
		pkl_file.close()

	U.print_separation('#...Finished Data Write',type=2)

## Wait until at least 60 seconds is up
duration = (float(time.asctime()[11:13])*3600+float(time.asctime()[14:16])*60+float(time.asctime()[17:19])) - (float(run_time[11:13])*3600+float(run_time[14:16])*60+float(run_time[17:19]))
if duration < 60:
	time.sleep(60-duration)
	

U.print_separation('## Finished millennium_stack.py'+'\n'+'Start:'+'\t'+run_time+'\n'+'End:'+'\t'+time.asctime(),type=1)


