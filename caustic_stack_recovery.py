"""
This program loads previously written data files from PBS batch job running on the stacking method.
The data is serialized using cPickle.
"""

## Import Modules ##
import cPickle as pkl
import numpy as np
import matplotlib.pyplot as mp
import numpy.ma as ma
import astStats
from mpl_toolkits.mplot3d import Axes3D
import sys
from AttrDict import AttrDict
import os.path
from caustic_universal_stack2D import universal
from caustic_class_stack2D import SelfStack,BinStack
import warnings
import scipy as sc
from numpy import random
import DictEZ as ez

## Flags ##
use_flux	= True				# Running on flux or sophie?
get_los		= True				# Upload Line of Sight Data as well?
data_loc	= 'binstack_run_table2'		# Parent directory where write_loc lives
write_loc	= 'ss_m1_run1'			# Which directory within data_loc to load ensembles from?
root		= '/glusterfs/users/caustics1/'

## Constants ##
warnings.filterwarnings("module",message="Warning: converting a masked element to nan.")


## Functions ##
class Recover(universal):
	'''This class contains functions that recover serialized data'''

	def __init__(self):
		pass


	def recover(self,write_loc=write_loc,write_stem=None,raw_data=False,ss=True,mm=False,go_global=True,ens_only=True,data_loc=None,bootstrap=False):

		"""
		This function uploads the pickle files from directory stack_data and configures them into multi dimensional arrays.
		It is meant to work with the self-stacked ensembles.
		go_global = True makes variables uploaded to global dictionary, False makes it returned to a dictionary
		write_loc = place where data lives
		raw_data: if True, output just mass estimates and caustic surfaces, no statistical calculations
		"""
		# For now, the "varib" dictionary and "HALODATA" array only need to be uploaded from one halo, say the first.
		# However, we must loop over the ensemble data to combine the stack_data arrays.

		## Load Data ##
		# Make wanted variables global
		self.ss = ss
		self.mm = mm

                # Create Final Dictionary w/ Data to Return
                final_data = {}

		# Create them as lists
		ENS_CAUMASS,ENS_CAUMASS_EST,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUMASS_EST,LOS_CAUSURF,LOS_NFWSURF = [],[],[],[],[],[],[],[]
		ENS_R,ENS_V,ENS_GMAGS,ENS_RMAGS,ENS_IMAGS,LOS_R,LOS_V,LOS_GMAGS,LOS_RMAGS,LOS_IMAGS,ENS_HVD,LOS_HVD = [],[],[],[],[],[],[],[],[],[],[],[]
		SAMS,PRO_POS,ENS_GP3D,ENS_GV3D,LOS_GP3D,LOS_GV3D = [],[],[],[],[],[]
		ENS_GAL_ID,LOS_GAL_ID,ENS_CLUS_ID,ENS_R3D,ENS_V3D = [],[],[],[],[]
		ENS_MEM = []
		BOOTSTRAP_SELECT = []

		# Initialization step 
		if data_loc==None:	
			if self.ss:	data_loc = 'old_selfstack_run_table'
			else:	data_loc = 'binstack_run_table2'

		pkl_file = open(root+'/nkern/OSDCStacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(0)+'_Data.pkl','rb')
		input = pkl.Unpickler(pkl_file)
		stack_data 	= input.load()
		varib		= input.load()
		run_dict	= input.load()

                # Add varib and run_dict to final_data
                final_data.update(varib)
                final_data.update(run_dict)

	        # Add varib to Classes  
		self.__dict__.update(varib)
		self.U = universal(varib)
		if self.ss:	self.halo_range = range(varib['halo_num'])
		else:	self.halo_range = range(varib['halo_num']/varib['line_num'])

		# Assign stack_data vars to object 'd'
		d = AttrDict(stack_data)

                ## Get Halo Data
                # Load and Sort Halos by Mass
                HaloID,HaloData = self.U.load_halos()
                HaloID,HaloData = self.U.sort_halos(HaloID,HaloData)
                HaloID,M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = np.vstack((HaloID,HaloData))
                HaloID = np.array(HaloID,int)
                # Build Halo_P, Halo_V
                Halo_P = np.vstack([HPX,HPY,HPZ])
                Halo_V = np.vstack([HVX,HVY,HVZ])

		## Make Bin Arrays if ss == False
		if self.ss == False:
			BinData = final_data['HaloData']
			if mm == True:
				mass_mix_match = final_data['mass_mix_match']
				HaloID,HaloID_init = final_data['HaloID'],final_data['HaloID_init']
				HaloData,HaloData_init = final_data['HaloData'],final_data['HaloData_init'][0:6]
				M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData
				BinData = HaloData_init.T[mass_mix_match].T
			BinData = self.U.Bin_Calc(BinData,varib,avg_meth=varib['avg_meth'])
			BIN_M200,BIN_R200,BIN_HVD = BinData


		# Append stack_data to major lists
		ENS_R.append(d.ens_r)
		ENS_V.append(d.ens_v)
		ENS_GMAGS.append(d.ens_gmags)
		ENS_RMAGS.append(d.ens_rmags)
		ENS_IMAGS.append(d.ens_imags)
		ENS_HVD.append(float(d.ens_hvd))
		ENS_CAUMASS.append(float(d.ens_caumass))
		ENS_CAUMASS_EST.append(float(d.ens_caumass_est))
		ENS_CAUSURF.append(d.ens_causurf)
		ENS_NFWSURF.append(d.ens_nfwsurf)
		LOS_R.append(d.los_r)
		LOS_V.append(d.los_v)
		LOS_GMAGS.append(d.los_gmags)
		LOS_RMAGS.append(d.los_rmags)
		LOS_IMAGS.append(d.los_imags)
		LOS_HVD.append(d.los_hvd)
		LOS_CAUMASS.append(d.los_caumass)
		LOS_CAUMASS_EST.append(d.los_caumass_est)
		LOS_CAUSURF.append(d.los_causurf)
		LOS_NFWSURF.append(d.los_nfwsurf)
		SAMS.append(d.sample_size)
		PRO_POS.append(d.pro_pos)
		ENS_GP3D.append(d.ens_gp3d)
		ENS_GV3D.append(d.ens_gv3d)
		LOS_GP3D.append(d.los_gp3d)
		LOS_GV3D.append(d.los_gv3d)
		ENS_GAL_ID.append(d.ens_gal_id)
		LOS_GAL_ID.append(d.los_gal_id)
		ENS_CLUS_ID.append(d.ens_clus_id)

		# Load 3D Data
		if bootstrap == False:
			ens_r3d,ens_v3d,ens_clus_index,ens_member = self.get_RV3D(ENS_GP3D,ENS_GV3D,ENS_CLUS_ID,HaloID,R_crit200,HVD,k=0)
			ENS_R3D.append(ens_r3d)
			ENS_V3D.append(ens_v3d)	
			ENS_MEM.append(ens_member)

		# Loop over ensembles
		j = 0
		for i in self.halo_range[1:]:
			# Progress Bar
			sys.stdout.write("Progress... "+str(j)+"\r")
			sys.stdout.flush()
			j += 1
			pkl_file = open(root+'/nkern/OSDCStacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(i)+'_Data.pkl','rb')
			input = pkl.Unpickler(pkl_file)
			stack_data = input.load()
				

			# Update locals() dictionary with stack_data dictionary			
			d = AttrDict(stack_data)

			# Append stack_data to major lists
			ENS_R.append(d.ens_r)
			ENS_V.append(d.ens_v)
			ENS_GMAGS.append(d.ens_gmags)
			ENS_RMAGS.append(d.ens_rmags)
			ENS_IMAGS.append(d.ens_imags)
			ENS_HVD.append(float(d.ens_hvd))
			ENS_CAUMASS.append(float(d.ens_caumass))
			ENS_CAUMASS_EST.append(float(d.ens_caumass_est))
			ENS_CAUSURF.append(d.ens_causurf)
			ENS_NFWSURF.append(d.ens_nfwsurf)
			LOS_R.append(d.los_r)
			LOS_V.append(d.los_v)
			LOS_GMAGS.append(d.los_gmags)
			LOS_RMAGS.append(d.los_rmags)
			LOS_IMAGS.append(d.los_imags)
			LOS_HVD.append(d.los_hvd)
			LOS_CAUMASS.append(d.los_caumass)
			LOS_CAUMASS_EST.append(d.los_caumass_est)
			LOS_CAUSURF.append(d.los_causurf)
			LOS_NFWSURF.append(d.los_nfwsurf)
			SAMS.append(d.sample_size)
			PRO_POS.append(d.pro_pos)
			ENS_GP3D.append(d.ens_gp3d)
			ENS_GV3D.append(d.ens_gv3d)
			LOS_GP3D.append(d.los_gp3d)
			LOS_GV3D.append(d.los_gv3d)
			ENS_GAL_ID.append(d.ens_gal_id)
			LOS_GAL_ID.append(d.los_gal_id)
			ENS_CLUS_ID.append(d.ens_clus_id)

			# Load 3D Data
			#if bootstrap == False:
			#	ens_r3d,ens_v3d,ens_clus_index,ens_member = self.get_RV3D(ENS_GP3D,ENS_GV3D,ENS_CLUS_ID,HaloID,R_crit200,HVD,k=i)
			#	ENS_R3D.append(ens_r3d)
			#	ENS_V3D.append(ens_v3d)	
			#	ENS_MEM.append(ens_member)

		print ''

		# Convert to arrays
		ENS_R = np.array(ENS_R)
		ENS_V = np.array(ENS_V)
		ENS_GMAGS = np.array(ENS_GMAGS)
		ENS_RMAGS = np.array(ENS_RMAGS)
		ENS_IMAGS = np.array(ENS_IMAGS)
		ENS_HVD = np.array(ENS_HVD)
		ENS_CAUMASS = np.array(ENS_CAUMASS)
		ENS_CAUMASS_EST = np.array(ENS_CAUMASS_EST)
		ENS_CAUSURF = np.array(ENS_CAUSURF)
		ENS_NFWSURF = np.array(ENS_NFWSURF)
		LOS_R = np.array(LOS_R)
		LOS_V = np.array(LOS_V)
		LOS_GMAGS = np.array(LOS_GMAGS)
		LOS_RMAGS = np.array(LOS_RMAGS)
		LOS_IMAGS = np.array(LOS_IMAGS)
		LOS_HVD = np.array(LOS_HVD,object)
		LOS_CAUMASS = np.array(LOS_CAUMASS,object)
		LOS_CAUMASS_EST = np.array(LOS_CAUMASS_EST,object)
		LOS_CAUSURF = np.array(LOS_CAUSURF,object)
		LOS_NFWSURF = np.array(LOS_NFWSURF,object)
		SAMS = np.array(SAMS)
		PRO_POS = np.array(PRO_POS)
		ENS_GP3D = np.array(ENS_GP3D)
		ENS_GV3D = np.array(ENS_GV3D)
		LOS_GP3D = np.array(LOS_GP3D)
		LOS_GV3D = np.array(LOS_GV3D)
		ENS_GAL_ID = np.array(ENS_GAL_ID)
		LOS_GAL_ID = np.array(LOS_GAL_ID)
		ENS_CLUS_ID = np.array(ENS_CLUS_ID)
		ENS_R3D = np.array(ENS_R3D)
		ENS_V3D = np.array(ENS_V3D)
		ENS_MEM = np.array(ENS_MEM)
		if bootstrap == True:
			BOOTSTRAP_SELECT = np.array(BOOTSTRAP_SELECT)

		if raw_data == True:
			# Return to Namespaces depending on go_global
			names = ['varib','HaloID','Halo_P','Halo_V','ENS_GP3D','ENS_GV3D','LOS_GP3D','LOS_GV3D','M_crit200','R_crit200','Z','SRAD','ESRAD','HVD','x_range','ENS_CAUMASS','ENS_CAUMASS_EST','ENS_CAUSURF','ENS_NFWSURF','LOS_CAUMASS','LOS_CAUMASS_EST','LOS_CAUSURF','LOS_NFWSURF','ENS_R','ENS_V','ENS_GMAGS','ENS_RMAGS','ENS_IMAGS','LOS_R','LOS_V','LOS_GMAGS','LOS_RMAGS','LOS_IMAGS','ENS_HVD','LOS_HVD','SAMS','PRO_POS','ENS_GAL_ID','LOS_GAL_ID','ENS_CLUS_ID','ENS_R3D','ENS_V3D','ENS_MEM','BIN_M200','BIN_R200','BIN_HVD','BOOSTRAP_SELECT']
			mydict = ez.create(names,locals())

                        # Add to final_data
                        final_data.update(mydict)

			# Return either dictionary of dump to globals()
			if go_global == True:
				globals().update(final_data)
				return
			elif go_global == False:
				return  final_data

		################################
		### Statistical Calculations ###
		################################
		if self.ss:	# If self stacking is True
			M200 = M_crit200[0:self.halo_num]
			HVD200 = HVD[0:self.halo_num]	
			ENS_MFRAC,ens_mbias,ens_mscat,ENS_VFRAC,ens_vbias,ens_vscat = self.stat_calc(ENS_CAUMASS,M200,ENS_HVD,HVD200)
			if ens_only == True:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = None,None,None,None,None,None
			else:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = self.stat_calc(LOS_CAUMASS.ravel(),M_crit200[0:self.halo_num],LOS_HVD.ravel(),HVD[0:self.halo_num],ens=False)
		else:		# If self stacking is False
			if mm == True:
				ENS_MFRAC,ens_mbias,ens_mscat,ENS_VFRAC,ens_vbias,ens_vscat = self.stat_calc(ENS_CAUMASS,BIN_M200,ENS_HVD,BIN_HVD,data_set=None)#'cut_low_mass')
			else:	
				ENS_MFRAC,ens_mbias,ens_mscat,ENS_VFRAC,ens_vbias,ens_vscat = self.stat_calc(ENS_CAUMASS,BIN_M200,ENS_HVD,BIN_HVD)	
			if ens_only == True:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = None,None,None,None,None,None
			else:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = self.stat_calc(LOS_CAUMASS.ravel(),M_crit200[0:self.halo_num],LOS_HVD.ravel(),HVD[0:self.halo_num],ens=False)


		# Create a dictionary
		names = ['varib','HaloID','Halo_P','Halo_V','ENS_GP3D','ENS_GV3D','LOS_GP3D','LOS_GV3D','M_crit200','R_crit200','Z','SRAD','ESRAD','HVD','x_range','ENS_CAUMASS','ENS_CAUMASS_EST','ENS_CAUSURF','ENS_NFWSURF','LOS_CAUMASS','LOS_CAUMASS_EST','LOS_CAUSURF','LOS_NFWSURF','ENS_R','ENS_V','ENS_GMAGS','ENS_RMAGS','ENS_IMAGS','LOS_R','LOS_V','LOS_GMAGS','LOS_RMAGS','LOS_IMAGS','ENS_HVD','LOS_HVD','SAMS','PRO_POS','ENS_GAL_ID','LOS_GAL_ID','ENS_CLUS_ID','ens_mbias','ens_mscat','los_mbias','los_mscat','ens_vbias','ens_vscat','los_vbias','los_vscat','ENS_MFRAC','ENS_VFRAC','LOS_MFRAC','LOS_VFRAC','ENS_R3D','ENS_V3D','ENS_MEM','BIN_M200','BIN_R200','BIN_HVD','BOOSTRAP_SELECT']
		mydict = ez.create(names,locals())

		# Append to final_data
		final_data.update(mydict)

		# Return either a dictionary or dump to globals()
		if go_global == True:
			globals().update(final_data)
			return
		elif go_global == False:
			return  final_data



	def stat_calc(self,MASS_EST,MASS_TRUE,HVD_EST,HVD_TRUE,data_set='full',ens=True):
		''' Does bias and scatter calculations '''
                # Cut data set if necessary
                if data_set == 'cut_low_mass':
                        '''Cutting all 'true' mass estimates below 1e14 off'''
                        cut = np.where(MASS_TRUE>1e14)[0]
                        MASS_EST = MASS_EST[cut]
                        MASS_TRUE = MASS_TRUE[cut]
                        HVD_EST = HVD_EST[cut]
                        HVD_TRUE = HVD_TRUE[cut]

		# Define a Masked array for sometimes zero terms
		epsilon = 10.0
		use_est = False				# Use MassCalc estimated r200 mass values if true 
		maMASS_EST	= ma.masked_array(MASS_EST,mask=MASS_EST<epsilon)		# Mask essentially zero values
		maHVD_EST	= ma.masked_array(HVD_EST,mask=HVD_EST<epsilon)


		# Mass / HVD Fractions
		if ens == True:
			# Ensemble Arrays
			MFRAC = np.log(maMASS_EST/MASS_TRUE)
			VFRAC = np.log(maHVD_EST/HVD_TRUE)
		else:
			# LOS Mass Fraction Arrays: 0th axis is halo number, 1st axis is line of sight number
			MFRAC,VFRAC = [],[]
			for a in range(len(MASS_EST)):
				MFRAC.append( ma.log( maMASS_EST[a]/MASS_TRUE[a] ) )
				VFRAC.append( ma.log( maHVD_EST[a]/HVD_TRUE[a] ) )
			MFRAC,VFRAC = np.array(MFRAC),np.array(VFRAC)

		if ens == True:
			mbias,mscat = astStats.biweightLocation(MFRAC,6.0),astStats.biweightScale(MFRAC,9.0)
			vbias,vscat = astStats.biweightLocation(VFRAC,6.0),astStats.biweightScale(VFRAC,9.0)
			return MFRAC,mbias,mscat,VFRAC,vbias,vscat
		else:
			if self.ss:
				# Create vertically averaged (by halo averaged) arrays, with line_num elements
				# biweightLocation takes only arrays with 4 or more elements
				HORZ_MFRAC,HORZ_VFRAC = [],[]
				VERT_MFRAC,VERT_VFRAC = [],[]
				for a in range(self.line_num):
					if len(ma.compressed(MFRAC[:,a])) > 4:
						VERT_MFRAC.append( astStats.biweightLocation( ma.compressed( MFRAC[:,a] ), 6.0 ) )
						VERT_VFRAC.append( astStats.biweightLocation( ma.compressed( VFRAC[:,a] ), 6.0 ) )
					else:
						VERT_MFRAC.append( np.median( ma.compressed( MFRAC[:,a] ) ) )
						VERT_VFRAC.append( np.median( ma.compressed( VFRAC[:,a] ) ) )
				VERT_MFRAC,VERT_VFRAC = np.array(VERT_MFRAC),np.array(VERT_VFRAC)
				# Create horizontally averaged (by line of sight) arrays, with halo_num elements
				for a in self.halo_range:
					if len(ma.compressed(MFRAC[a])) > 4:
						HORZ_MFRAC.append( astStats.biweightLocation( ma.compressed( MFRAC[a] ), 6.0 ) )
						HORZ_VFRAC.append( astStats.biweightLocation( ma.compressed( VFRAC[a] ), 6.0 ) )
					else:
						HORZ_MFRAC.append( np.median( ma.compressed( MFRAC[a] ) ) )
						HORZ_VFRAC.append( np.median( ma.compressed( VFRAC[a] ) ) )
				HORZ_MFRAC,HORZ_VFRAC = np.array(HORZ_MFRAC),np.array(HORZ_VFRAC)
				# Bias and Scatter Calculations
				mbias,mscat = astStats.biweightLocation(VERT_MFRAC,6.0),astStats.biweightScale(VERT_MFRAC,9.0)
				vbias,vscat = astStats.biweightLocation(VERT_VFRAC,6.0),astStats.biweightScale(VERT_VFRAC,9.0)
			else:
				# Bin stack LOS systems need only one average
				mbias,mscat = astStats.biweightLocation(MFRAC,6.0),astStats.biweightScale(MFRAC,9.0)
				vbias,vscat = astStats.biweightLocation(VFRAC,6.0),astStats.biweightScale(VFRAC,9.0)

			return MFRAC,mbias,mscat,VFRAC,vbias,vscat
	

	def get_RV3D(self,ENS_GP3D,ENS_GV3D,ENS_CLUS_ID,HaloID,R_crit200,HVD,k=0):
		''' 
		This function takes 3D position and velocities of ensemble galaxies to compute 3D R and V vectors
		Additionally, it takes the R200 of each individual cluster, and assigns a numeric to each galaxy,
		either a 0 or 1, meaning it is not a member (0) or is a member (1) of its host cluster.
		This is done by the galaxy being within 3D R200 radius.
		k : This is the halo index currently being worked on, w.r.t. final HaloID array.
		'''
		# First Create an Indexing Array that assigns each galaxy to its host cluster
		ens_clus_id = ENS_CLUS_ID[k]
		ens_gp3d = ENS_GP3D[k]
		ens_gv3d = ENS_GV3D[k]
		ens_clus_index = np.array([],int)
		for i in range(len(ens_clus_id)):
			ens_clus_index = np.append(ens_clus_index,np.where(HaloID==ens_clus_id[i])[0][0])

		# Create 3D R and V arrays	
		ens_r3d = np.sqrt(ens_gp3d[0]**2 + ens_gp3d[1]**2 + ens_gp3d[2]**2)
		ens_v3d = np.sqrt(ens_gv3d[0]**2 + ens_gv3d[1]**2 + ens_gv3d[2]**2) / np.sqrt(3)

		# Make Membership Array
		ens_member = np.array([],int)
		for i in range(len(ens_clus_index)):
			radius_ratio = ens_r3d[i]/R_crit200[ens_clus_index[i]]
			if radius_ratio <= 1.0:
				ens_member = np.append(ens_member,[1])
			else:
				ens_member = np.append(ens_member,[0])

		return ens_r3d,ens_v3d,ens_clus_index,ens_member	




class Work(Recover):
	'''This class contains functions that works with the data previously loaded'''

	def __init__(self,Recover):
		pass	
	
	def append_data(self,kwargs,i):
		''' This function was created so as to reclaim the mydict dictionary memory after exiting the function.'''
		
		# Load in Data from Run Table and append
		mydict = self.recover(**kwargs)
		d = AttrDict(mydict)
		self.RUN_NUM.append(i)
		self.GAL_NUM.append(d.gal_num)
		self.LINE_NUM.append(d.line_num)
		self.ENS_MBIAS.append(d.ens_mbias)
		self.ENS_MSCAT.append(d.ens_mscat)
		self.ENS_VBIAS.append(d.ens_vbias)
		self.ENS_VSCAT.append(d.ens_vscat)
		self.LOS_MBIAS.append(d.los_mbias)
		self.LOS_MSCAT.append(d.los_mscat)
		self.LOS_VBIAS.append(d.los_vbias)
		self.LOS_VSCAT.append(d.los_vscat)	


	def load_all(self,iter_array=None,tab_shape=None,ens_only=True,kwargs=None,write_stem=None):
		'''
		This iterates over different richness geometry configurations and runs statistics on data.
		It is recommended to do any calculations (statistics, plots etc.) within the for loop, and
		then feed results back out via a global variable.
		'''
		# Feed Local Variables
		self.ens_only = ens_only	
		if kwargs == None:
			kwargs = {'write_loc':'mm_m0_run1','raw_data':False,'ss':False,'go_global':False,'ens_only':True,'data_loc':'mass_mix/mm_0.05_run_table1'}
		if write_stem == None:
			write_stem = 'mm_m0_run'
	
		# Configure Variables
		if iter_array == None:
			iter_array = np.arange(1,50)
			tab_shape = (7,7)

		
		## Calculate Ensemble Only Statistics
		if self.ens_only == True:

			# Create arrays
			self.ENS_MBIAS,self.ENS_MSCAT,self.ENS_VBIAS,self.ENS_VSCAT = [],[],[],[]
			self.LOS_MBIAS,self.LOS_MSCAT,self.LOS_VBIAS,self.LOS_VSCAT = [],[],[],[]
			self.RUN_NUM,self.GAL_NUM,self.LINE_NUM = [],[],[]

			# Iterate over runs for ensemble data
			print '...Loading Data from '+str(len(iter_array))+' runs'
			for i in iter_array:
				print ''
				print 'Working on Run #'+str(i)
				print '-'*25
				## Define Recover Keyword Arguments!  ##
				kwargs['write_loc'] = write_stem+str(i)
				print 'Recover Keyword Arguments:'
				print '-'*40
				print kwargs
				
				## Load and Append Data		
				self.append_data(kwargs,i)

			# Make into arrays that resemble table
			print 'Table Shape =',tab_shape
			ENS_MBIAS,ENS_MSCAT,ENS_VBIAS,ENS_VSCAT = np.array(self.ENS_MBIAS).reshape(tab_shape),np.array(self.ENS_MSCAT).reshape(tab_shape),np.array(self.ENS_VBIAS).reshape(tab_shape),np.array(self.ENS_VSCAT).reshape(tab_shape)
			LOS_MBIAS,LOS_MSCAT,LOS_VBIAS,LOS_VSCAT = np.array(self.LOS_MBIAS).reshape(tab_shape),np.array(self.LOS_MSCAT).reshape(tab_shape),np.array(self.LOS_VBIAS).reshape(tab_shape),np.array(self.LOS_VSCAT).reshape(tab_shape)
			RUN_NUM,GAL_NUM,LINE_NUM = np.array(self.RUN_NUM).reshape(tab_shape),np.array(self.GAL_NUM).reshape(tab_shape),np.array(self.LINE_NUM).reshape(tab_shape) 

			# Other Data Arrays
			RICH_NUM = GAL_NUM*LINE_NUM

			return (ENS_MBIAS,ENS_MSCAT,ENS_VBIAS,ENS_VSCAT,LOS_MBIAS,LOS_MSCAT,LOS_VBIAS,LOS_VSCAT,RUN_NUM,GAL_NUM,LINE_NUM,RICH_NUM)


		else:
			pass


	def bootstrap_load_write(self):
		''' 
		Performs a load_all() on bootstrap run tables over all repetitions, then writes out pertinent data
		'''
		pass


	def oto(self,xarray,yarray,style='ko',alpha=None):
		'''Simple log log one to one plot setup'''
		p1, = mp.plot(xarray,yarray,style,alpha=alpha)
		mp.plot([xarray[0],xarray[-1]],[xarray[0],xarray[-1]],'b')
		mp.xscale('log')
		mp.yscale('log')
		return p1


	def get_3d(self):
		pass
	
	def sample_histogram(self,caumass,truemass,bins=20,ax=None):
		if ax == None:
			mp.hist(caumass,bins=bins,color='b',alpha=.6)
			p1 = mp.axvline(truemass,ymin=.8,color='k',lw=1.5)
			p2 = mp.axvline(np.median(caumass*1),ymin=.8,color='g',lw=1.5)
			p3 = mp.axvline(np.mean(caumass*1),ymin=.8,color='c',lw=1.5)
			p4 = mp.axvline(astStats.biweightLocation(caumass*1,6.0),ymin=.8,color='r',lw=1.5)
			mp.legend([p1,p2,p3,p4],["Table Mass","Median","Mean","biweightLocation"])
		else:
			ax.hist(caumass,bins=bins,color='b',alpha=.6)
			p1 = ax.axvline(truemass,ymin=.8,color='k',lw=1.5)
			p2 = ax.axvline(np.median(caumass*1),ymin=.8,color='g',lw=1.5)
			p3 = ax.axvline(np.mean(caumass*1),ymin=.8,color='c',lw=1.5)
			p4 = ax.axvline(astStats.biweightLocation(caumass*1,6.0),ymin=.8,color='r',lw=1.5)
			ax.legend([p1,p2,p3,p4],["Table Mass","Median","Mean","biweightLocation"],fontsize=8)






###############################################################
############## END CLASSES, BEGIN PROGRAM #####################
###############################################################


## Initialize Classes
R = Recover()
W = Work(Recover)

work = True
if work == True:
	
	table_num	= str(sys.argv[1])
	iter_array	= np.arange(1,35)
	tab_shape	= (4,7)
	write_stem	= 'ss_m1_run'
	write_loc	= 'ss_m1_run1'
	data_loc	= 'selfstack/ss_run_table'+str(table_num)
	ss		= True
	mm		= False
	

	kwargs = {'write_loc':write_loc,'raw_data':False,'ss':ss,'mm':mm,'go_global':False,'ens_only':True,'data_loc':data_loc}

	data = W.load_all(kwargs=kwargs,write_stem=write_stem,iter_array=iter_array,tab_shape=tab_shape)

	names = ('ENS_MBIAS','ENS_MSCAT','ENS_VBIAS','ENS_VSCAT','LOS_MBIAS','LOS_MSCAT','LOS_VBIAS','LOS_VSCAT','RUN_NUM','GAL_NUM','LINE_NUM','RICH_NUM')

	values = np.copy(data)

	dictionary = dict(zip(names,values))

	file = open(root+'/nkern/OSDCStacking/'+data_loc+'/table_analysis.pkl','wb')

	output = pkl.Pickler(file)

	output.dump(dictionary)









