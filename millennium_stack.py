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
keys = ['c','h','H0','q','beta','fbeta','r_limit','v_limit','data_set','halo_num','gal_num','line_num','method_num','write_loc','data_loc','root','self_stack','scale_data','write_data','run_time','clean_ens','small_set','run_los','bootstrap','run_num','ens_num','cell_num','stack_range','mass_mix','mass_scat','bootstrap_num','bootstrap_rep','avg_meth','cent_offset','center_scat','new_halo_cent','true_mems']
varib = ez.create(keys,locals())

## INITIALIZATION ##
S = Stack(varib)
U = Universal(varib)
A = Analysis(varib)
M = Millennium(varib)

U.print_separation('## Running caustic_mass_stack.py')
U.print_varibs(varib)

## Load Halo Data
U.print_separation('# ...Loading Halos',type=2)
HaloID,HaloData = M.load_halos()

# Sort Halos by A Priori Known Descending Mass (Mass Critical 200)
HaloID,HaloData = M.sort_halos(HaloID,HaloData)

####################### END STANDARD PREAMBLE ##############################

# Unpack HaloData array into local namespace
M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData

## Load Galaxy Data
U.print_separation('# ...Loading Galaxies',type=2)
Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,HaloData = M.configure_galaxies(HaloID,HaloData)
M_crit200,R_crit200,Z,HVD = HaloData

# Get Gal_P in physical coordinates
Gal_P2 = []
if self_stack == True:
	for [i,k] in zip(np.arange(ens_num),stack_range):
		Gal_P2.append((Gal_P[i].T-Halo_P[k]).T)
else:
	for [i,l] in zip(np.arange(ens_num*line_num),stack_range):
		Gal_P2.append((Gal_P[i].T-Halo_P[l]).T)


# If Bin Stacking:
#	- Create Ensemble R200 and HVD arrays
#	- Change order of halos to bin upon
#	- Create any other arrays needed
if self_stack == False:
	U.print_separation("Average Method for Construction of Bin Properties is "+avg_meth,type=2)
	BinData = U.Bin_Calc(HaloData,varib,avg_meth=avg_meth)
	BIN_M200,BIN_R200,BIN_HVD = BinData

# Initialize Multi-Ensemble Array to hold resultant data
STACK_DATA = []

U.print_separation('# ...Starting Ensemble Loop',type=2)
#  j: index for each ensemble, w.r.t. the FOR loop
#  k: index for each final ensemble, w.r.t. total number of final ensembles
#  l: index for line of sight, w.r.t. total lines of sight of run

# Define Container to Hold Phase Space Data 
for j in range(ens_num):

	# Total Halo Index
	k = run_num*ens_num + j

	# Define Container to Hold Phase Space Data
	PS = Data()

	# Iterate through lines of sight
	ens_gal_count = 0
	for l in range(line_num):
		if self_stack == True:
			# Do Projection
			r, v, projected_pos = U.line_of_sight(Gal_P[j],Gal_V[j],Halo_P[stack_range][j],Halo_V[stack_range][j],project=False)	

			#  


		else:
			# Do Projection


	# Build Ensemble and Run Caustic Technique
	stack_data = S.caustic_stack()









	if self_stack:
		if bootstrap == True:
			stack_data = SS.self_stack_bootstrap(HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j)
		else:
			stack_data = SS.self_stack_clusters(HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j)
	else:
		if bootstrap == True:
			stack_data = BS.bin_stack_bootstrap(HaloID,HaloData,BinData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j)
		else:
			stack_data = BS.bin_stack_clusters(HaloID,HaloData,BinData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j)

	# Unpack data
	ens_r,ens_v,ens_gal_id,ens_clus_id,ens_gmags,ens_rmags,ens_imags,ens_hvd,ens_caumass,ens_caumass_est,ens_causurf,ens_nfwsurf,los_r,los_v,los_gal_id,los_gmags,los_rmags,los_imags,los_hvd,los_caumass,los_caumass_est,los_causurf,los_nfwsurf,x_range,sample_size,pro_pos,ens_r200_est,vel_avg = stack_data

	# Operate on Center Offset Data if Available
	if cent_offset != None:
		offset_data = BS.bin_stack_clusters(HaloID,HaloData,BinData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j,PRO_POS=pro_pos,VEL_AVG=vel_avg)
		off_ens_r,off_ens_v,off_ens_gal_id,off_ens_clus_id,off_ens_gmags,off_ens_rmags,off_ens_imags,off_ens_hvd,off_ens_caumass,off_ens_caumass_est,off_ens_causurf,off_ens_nfwsurf,off_los_r,off_los_v,off_los_gal_id,off_los_gmags,off_los_rmags,off_los_imags,off_los_hvd,off_los_caumass,off_los_caumass_est,off_los_causurf,off_los_nfwsurf,off_x_range,off_sample_size,off_pro_pos,off_ens_r200_est,off_vel_avg = offset_data

	# Get 3D data
	if self_stack == True:
		if bootstrap == True:
			pass # Don't yet know how to do this for bootstrap
		else:	
			ens_gp3d,ens_gv3d,los_gp3d,los_gv3d = U.get_3d(np.array(Gal_P2),np.array(Gal_V),ens_gal_id,los_gal_id,stack_range,ens_num,self_stack,j)	
	else:
		if bootstrap == True:
			pass # Don't yet know how to do this for bootstrap
		else:
			ens_gp3d,ens_gv3d,los_gp3d,los_gv3d = U.get_3d(np.array(Gal_P2),np.array(Gal_V),ens_gal_id,los_gal_id,stack_range,ens_num,self_stack,j)

	# Combine into stack_data
	if run_los == False:
                keys = ['ens_r','ens_v','ens_gal_id','ens_clus_id','ens_gmags','ens_rmags','ens_imags','ens_hvd','ens_caumass','ens_caumass_est','ens_causurf','ens_nfwsurf','x_range','sample_size','pro_pos','ens_gp3d','ens_gv3d','BS.bootstrap_select','ens_r200_est','vel_avg']
	
	else:
		keys = ['ens_r','ens_v','ens_gal_id','ens_clus_id','ens_gmags','ens_rmags','ens_imags','ens_hvd','ens_caumass','ens_caumass_est','ens_causurf','ens_nfwsurf','los_r','los_v','los_gal_id','los_gmags','los_rmags','los_imags','los_hvd','los_caumass','los_caumass_est','los_causurf','los_nfwsurf','x_range','sample_size','pro_pos','ens_gp3d','ens_gv3d','los_gp3d','los_gv3d','BS.bootstrap_select','ens_r200_est','vel_avg']
	stack_data = ez.create(keys,locals())
	
	if cent_offset != None:
		keys = ['off_ens_r','off_ens_v','off_ens_hvd','off_ens_caumass','off_ens_caumass_est','off_ens_causurf','off_los_r','off_los_v','off_los_hvd','off_los_caumass','off_los_caumass_est','off_los_causurf']
		off_dict = ez.create(keys,locals())
		stack_data.update(off_dict)

	# Append to STACK_DATA
	STACK_DATA.append(stack_data)


# Finished Loop
U.print_separation('#...Finished Ensemble Loop',type=2)

# Create run_dict
keys = ['HaloID','HaloData','HaloID_init','HaloData_init','mass_mix_match','M_crit200_match']
run_dict = ez.create(keys,locals())

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
	

U.print_separation('## Finished caustic_mass_stack2D.py'+'\n'+'Start:'+'\t'+run_time+'\n'+'End:'+'\t'+time.asctime(),type=1)


