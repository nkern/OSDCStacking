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
		data_loc = 'mass_mix/mm_run_table'+str(table_num)
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
names = ['run_time','','run_num','gal_num','line_num','cell_num','ens_num','halo_num','method_num','avg_meth','','self_stack','lightcone','mass_mix','mm_est','write_data','scale_data','run_los','mirror','new_halo_cent','cent_offset','true_mems','init_clean','bootstrap','','mass_scat','center_scat','bootstrap_num','bootstrap_rep','','data_loc','write_loc','root']
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

# Unpack HaloData array into local namespace
M_crit200,R_crit200,HVD,Z,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData
Halo_P,Halo_V = np.vstack([HPX,HPY,HPZ]).T,np.vstack([HVX,HVY,HVZ]).T

############### END STANDARD PREAMBLE ######################

# Load in Halo Data from previously written file that has phase space data that has been mass mixed
try:
	print '...Loading from halo_arrays.pkl'
	f = open(root+'/OSDCStacking/'+data_loc+'/halo_arrays.pkl','rb')
	input = pkl.Unpickler(f)	
	run_dict = input.load()
	globals().update(run_dict)
	f.close()

# Do Mass Mixing if Required
except IOError:
	if mass_mix == True and sys.argv[-1] == 'mm':
		# Check to see if correct path
		if os.path.exists(root+'/OSDCStacking/'+data_loc): pass
		else: print data_loc+" Doesn't exist!"; raise NameError
		# Start Mass Mixing
		print "...Do Mass Mixing for "+str(halo_num)+" clusters using "+mm_est+" estimator"
		richness = []
		Rdata,Vdata,Gmags,Rmags,Imags,Pro_Pos = [],[],[],[],[],[]
		sys.stdout.write('...Mass Mixing\n')
		n = -1
		for i in range(len(HaloID)):
			if i%(len(HaloID)/100) == 0: sys.stdout.write(str(n)+'%\r'); sys.stdout.flush(); n += 1

			galdata = M.configure_galaxies(HaloID[i],HaloData.T[i])
			if lightcone == True:
				# Get R and V data
				gal_ra,gal_dec,gal_z,gmags,rmags,imags,gal_p,gal_v = galdata
				clus_ra = RA[i]
				clus_dec = DEC[i]
				clus_z = Z[i]
				ang_d,lum_d = S.C.zdistance(clus_z,H0)
				angles = S.C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
				rdata = angles * ang_d
				vdata = c * (gal_z - clus_z) / (1 + clus_z)
				# Make a rough cut
				cut = np.where( (rdata < 7) & (rmags < -18.0) )[0]
				rdata,vdata,gmags,rmags,imags = rdata[cut],vdata[cut],gmags[cut],rmags[cut],imags[cut]
				# Append and get richness
				Rdata.append(rdata)
				Vdata.append(vdata)
				Gmags.append(gmags)
				Rmags.append(rmags)
				Imags.append(imags)
				richness.append(M.richness_est(rdata,vdata,gmags,rmags,imags))
			else:
				# Get R and V data
				gal_p,gal_v,gmags,rmags,imags = galdata
				rdata, vdata, pro_pos = U.line_of_sight(gal_p,gal_v,Halo_P[i],Halo_V[i])
				rdata = np.array(rdata)
				vdata = np.array(vdata)
				pro_pos = np.array(pro_pos)
				# Make rough cut
				cut = np.where( (rdata < 7) & (rmags < -18.0) )[0]
				rdata,vdata,gmags,rmags,imags = rdata[cut],vdata[cut],gmags[cut],rmags[cut],imags[cut]
				# Append and get richness
				Rdata.append(rdata)
				Vdata.append(vdata)
				Gmags.append(gmags)
				Rmags.append(rmags)
				Imags.append(imags)
				Pro_Pos.append(pro_pos)
				richness.append(M.richness_est(rdata,vdata,gmags,rmags,imags))

		Rdata,Vdata,Gmags,Rmags,Imags = np.array(Rdata),np.array(Vdata),np.array(Gmags),np.array(Rmags),np.array(Imags)
		Pro_Pos = np.array(Pro_Pos)
		richness = np.array(richness)

		# Sort by Richness
		if mm_est == 'richness':
			sort = np.argsort(richness)[::-1]
			richness = richness[sort]
			HaloID,HaloData = HaloID[sort],HaloData.T[sort].T
			Rdata,Vdata,Gmags,Rmags,Imags = Rdata[sort],Vdata[sort],Gmags[sort],Rmags[sort],Imags[sort]

		# Write Data file halo_arrays.pkl and output data for each halo
		try: os.system('mkdir '+root+'/OSDCStacking/'+data_loc+'/cluster_data')
		except: pass

		f = open(root+'/OSDCStacking/'+data_loc+'/halo_arrays.pkl','wb')
		output = pkl.Pickler(f)
		run_dict = {'HaloID':HaloID,'HaloData':HaloData}

		if lightcone == True: run_dict['RA'] = RA; run_dict['DEC'] = DEC
		else: run_dict['Pro_Pos'] = Pros_Pos

		if mm_est == 'richness': run_dict['richness'] = richness

		output.dump(run_dict)
		f.close()

		sys.stdout.write('\n...Writing Files\n')
		n = -1
		for i in range(len(HaloID)):
			if i%(len(HaloID)/100) == 0: sys.stdout.write(str(n)+'%\r'); sys.stdout.flush(); n += 1

			f = open(root+'/OSDCStacking/'+data_loc+'/cluster_data/'+str(HaloID[i])+'_RVdata.tab','w')
			f.write('#Rdata, Vdata, Gmags, Rmags, Imags\n')
			for j in range(len(Rdata[i])):
				f.write(str(Rdata[i][j])+'\t'+str(Vdata[i][j])+'\t'+str(Gmags[i][j])+'\t'+str(Rmags[i][j])+'\t'+str(Imags[i][j])+'\n')
			f.close()	

		U.print_separation('## Finished mass mixing!!'+'\n'+'Start:'+'\t'+run_time+'\n'+'End:'+'\t'+time.asctime(),type=1)
		raise NameError	

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

## Create Mass Mixing Functions and assign to "Stack" class

# Initialize Multi-Ensemble Array to hold resultant data
STACK_DATA = []

#  j: index for each ensemble, w.r.t. the FOR loop
#  k: index for each final ensemble, w.r.t. total number of final ensembles
#  l: index for line of sight, w.r.t. total lines of sight of run
U.print_separation('# ...Starting Ensemble Loop',type=2)

for j in range(ens_num):
	U.print_separation('## Working on Ensemble '+str(j),type=1)

	# Update S.j
	S.j = j

	# Total Halo Index
	k = run_num*ens_num + j

	# Load galaxy data from pre-written files
	if self_stack == True:
		Rdata,Vdata,G_Mags,R_Mags,I_Mags = [],[],[],[],[]
		for i in stack_range[j]:
			rdata,vdata,gmags,rmags,imags = np.loadtxt(root+'/OSDCStacking/'+data_loc+'/cluster_data/'+str(HaloID[i])+'_RVdata.tab',unpack=True)
			Rdata.append(rdata)
			Vdata.append(vdata)
			G_Mags.append(gmags)
			R_Mags.append(rmags)
			I_Mags.append(imags)
		Rdata,Vdata,G_Mags,R_Mags,I_Mags = np.array(Rdata),np.array(Vdata),np.array(G_Mags),np.array(R_Mags),np.array(I_Mags)	

	else:
		Rdata,Vdata,G_Mags,R_Mags,I_Mags = [],[],[],[],[]
		for i in stack_range[j*line_num:(j+1)*line_num]:
			rdata,vdata,gmags,rmags,imags = np.loadtxt(root+'/OSDCStacking/'+data_loc+'/cluster_data/'+str(HaloID[i])+'_RVdata.tab',unpack=True)
			Rdata.append(rdata)
			Vdata.append(vdata)
			G_Mags.append(gmags)
			R_Mags.append(rmags)
			I_Mags.append(imags)
		Rdata,Vdata,G_Mags,R_Mags,I_Mags = np.array(Rdata),np.array(Vdata),np.array(G_Mags),np.array(R_Mags),np.array(I_Mags)	

	# Build Ensemble and Run Caustic Technique
	stack_data = S.caustic_stack(Rdata,Vdata,HaloID[stack_range[j*line_num:(j+1)*line_num]],np.vstack([M_crit200[stack_range[j*line_num:(j+1)*line_num]],R_crit200[stack_range[j*line_num:(j+1)*line_num]],HVD[stack_range[j*line_num:(j+1)*line_num]]]),line_num,feed_mags=True,G_Mags=G_Mags,R_Mags=R_Mags,I_Mags=I_Mags)

	# Append other Arrays
	#if lightcone == False:
	#	extra = {'pro_pos':Pro_Pos,'_name_':'stack_data'}
	#	stack_data.update(extra)

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
