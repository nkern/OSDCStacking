## millennium_massmix.py
"""
This program uses the Caustic Technique to estimate the masses of galaxy clusters,
after having applied a stacking technique to create ensemble clusters.
Before it does so, however, it introduces a scatter into the "assumed" mass of
each individual cluster, producing a "mass mixing" effect into our final analysis.
"""

###################### BEGIN STANDARD PREAMBLE #####################
print '...importing modules'
from caustic_stack import *
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
		data_loc = 'mass_mix/mm_'+"%.2f" % mass_scat+'_run_table'+str(table_num)
		write_loc = 'mm_m'+str(method_num)+'_run'+str(cell_num)
		if cent_offset == 'r':
			data_loc = 'mass_mix/roff_run_table'+str(table_num)
		elif cent_offset == 'v':
			data_loc = 'mass_mix/voff_run_table'+str(table_num)
if bootstrap == True:
	write_loc = 'bo_m'+str(method_num)+'_run'+str(cell_num)
	data_loc = 'binstack/bootstrap'+str(bootstrap_num)+'/rep'+str(bootstrap_rep)

## Make dictionary for loaded constants, doesn't matter if some don't exist
keys = ['c','h','H0','q','beta','fbeta','r_limit','v_limit','data_set','halo_num','gal_num','line_num','method_num','write_loc','data_loc','root','self_stack','scale_data','write_data','run_time','init_clean','small_set','run_los','bootstrap','run_num','ens_num','cell_num','stack_range','mass_mix','mass_scat','bootstrap_num','bootstrap_rep','avg_meth','cent_offset','center_scat','new_halo_cent','true_mems','mirror','edge_perc','lightcone','mm_est']
varib = ez.create(keys,locals())
varib.update({'_name_':'varib'})

## INITIALIZATION ##
S = Stack(varib)
U = Universal(varib)
M = Millennium(varib)

U.print_separation('## Running millennium_massmix.py')
names = ['run_time','','run_num','gal_num','line_num','cell_num','ens_num','halo_num','method_num','avg_meth','','self_stack','lightcone','mass_mix','write_data','scale_data','run_los','mirror','new_halo_cent','cent_offset','true_mems','init_clean','bootstrap','','mass_scat','center_scat','bootstrap_num','bootstrap_rep','','data_loc','write_loc','root']
U.print_varibs(names,varib)

## Load Halo Data
U.print_separation('# ...Loading Halos',type=1)
if lightcone == True:
	HaloID,RA,DEC,HaloData = M.load_halos()
else:
	HaloID,HaloData = M.load_halos()

# Sort Halos by A Priori Known Descending Mass (Mass Critical 200)
if lightcone == True:
	HaloID,RA,DEC,HaloData = M.sort_halos(HaloID,HaloData)
else:
	HaloID,HaloData = M.sort_halos(HaloID,HaloData)

############### END STANDARD PREAMBLE ######################

# Load in Halo Data from previously written file, or write file if doesn't exist
try:
	f = open(root+'/OSDCStacking/'+data_loc+'/halo_arrays.pkl','rb')
	input = pkl.Unpickler(f)	
	run_dict = input.load()
	globals().update(run_dict)
	f.close()
except IOError:
	try:
		f = open(root+'/OSDCStacking/'+data_loc+'/halo_arrays.pkl','wb')	
		output = pkl.Pickler(f)	
		if lightcone == True:
			run_dict = {'HaloID':HaloID,'RA':RA,'DEC':DEC,'HaloData':HaloData}
		else:
			run_dict = {'HaloID':HaloID,'HaloData':HaloData}
		output.dump(run_dict)
		f.close()
	except IOError:
		pass

# Unpack HaloData array into local namespace
M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData
Halo_P,Halo_V = np.vstack([HPX,HPY,HPZ]).T,np.vstack([HVX,HVY,HVZ]).T

###############Outdated Mass Mixing Technique##################
#if mass_mix == True:
#	HaloID_init,HaloData_init = np.copy(HaloID),np.copy(HaloData)
#	try:
#		""" Need to create 1 mass mixed array that every job loads in to preserve consistency"""
#		# Mass Mix arrays
#		HaloID,HaloData,massmix_sort = M.mass_mixing(HaloID,HaloData,float(mass_scat))
#		# Create Assumed and True BinData, former sorted by mixed mass, latter sorted by table mass
#		BinData = U.Bin_Calc(HaloData[0][:halo_num],HaloData[1][:halo_num],HaloData[3][:halo_num])
#		TrueBinData = U.Bin_Calc(HaloData_init[0][:halo_num],HaloData_init[1][:halo_num],HaloData_init[3][:halo_num])
#		# Dump Data in to halo_arrays.pkl
#		f = open(root+'/OSDCStacking/'+data_loc+'/halo_arrays.pkl','wb')
#		output = pkl.Pickler(f)
#		data = {'massmix_sort':massmix_sort,'HaloID':HaloID,'HaloData':HaloData,'HaloID_init':HaloID_init,'HaloData_init':HaloData_init,'BinData':BinData,'TrueBinData':TrueBinData}
#		output.dump(data)
#		f.close()
#	except IOError:
#		f = open(root+'/OSDCStacking/'+data_loc+'/halo_arrays.pkl','rb')
#		input = pkl.Unpickler(f)
#		globals().update(input.load())
################################################################

# Initialize Multi-Ensemble Array to hold resultant data
STACK_DATA = []

#  j: index for each ensemble, w.r.t. the FOR loop
#  k: index for each final ensemble, w.r.t. total number of final ensembles
#  l: index for line of sight, w.r.t. total lines of sight of run
U.print_separation('# ...Starting Ensemble Loop',type=2)

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
		pkl_file.close()

	U.print_separation('#...Finished Data Write',type=2)

## Wait until at least 60 seconds is up
duration = (float(time.asctime()[11:13])*3600+float(time.asctime()[14:16])*60+float(time.asctime()[17:19])) - (float(run_time[11:13])*3600+float(run_time[14:16])*60+float(run_time[17:19]))
if duration < 60:
	time.sleep(60-duration)
	

U.print_separation('## Finished millennium_massmix.py'+'\n'+'Start:'+'\t'+run_time+'\n'+'End:'+'\t'+time.asctime(),type=1)
