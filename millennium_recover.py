## millennium_recover.py
"""
This program loads previously written data files from stacking code.
The data is serialized with CPickle.
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
import warnings
import scipy as sc
from numpy import random as npr
import DictEZ as ez
from stack_class import *

## Flags ##
root		= '/glusterfs/users/caustics1/nkern'


## Constants ##
warnings.filterwarnings("module",message="Warning: converting a masked element to nan.")


## Functions ##
class Recover(Universal):
	'''This class contains functions that recovers serialized data'''

	def __init__(self):
		pass


	def recover(self,write_loc=None,raw_data=False,ss=True,mm=False,go_global=True,ens_only=True,data_loc=None,bootstrap=False,avg_meth=None,cent_offset=None):

		"""
		This function uploads the pickle files from directory stack_data and configures them into multi dimensional arrays.
		It is meant to work with the self-stacked ensembles.
		go_global = True makes variables uploaded to global dictionary, False makes it returned to a dictionary
		write_loc = place where data lives
		raw_data : if True, output just mass estimates and caustic surfaces, no statistical calculations
		"""

		## Load Data ##
		# Attach certain variables to class
		self.ss = ss
		self.mm = mm

		# Create Open Container for Data
		D = Data()

		# Create list of variables you want to load in and stack into large array
		# Not all names in list need to be in stack_data, it will take those that exist
		load_names = [
		'ens_r','ens_v','ens_gmags','ens_rmags','ens_imags','ens_caumass','ens_caumass_est','ens_edgemass','ens_edgemass_est','ens_hvd',
		'ens_causurf','ens_edgesurf','ens_nfwsurf','ens_gal_id','ens_clus_id',
		'los_r','los_v','los_gmags','los_rmags','los_imags','los_caumass','los_caumass_est','los_edgemass','los_edgemass_est','los_hvd',
		'los_causurf','los_edgesurf','los_nfwsurf','los_gal_id',
		'pro_pos','bootstrap_select','ens_r200_est',
		'BinM200','BinR200','BinHVD'
		]

		# Open First Data file and load in varib dictionary
		pkl_file = open(root+'/OSDCStacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(0)+'_Data.pkl','rb')
		input = pkl.Unpickler(pkl_file)
		stack_data 	= input.load()
		varib		= input.load()

		# Add varib and first Ensemble and other single arrays to Data
		D.__dict__.update(varib)
		D.append(stack_data,keys=load_names)
		D.add({'x_range':stack_data['x_range']})

	        # Add varib to Classes  
		self.__dict__.update(varib)
		self.U = Universal(varib)
		self.M = Millennium(varib)

		# Define halo_range: number of ensembles to load
		if self.ss:	self.halo_range = range(varib['halo_num'])
		else:		self.halo_range = range(varib['halo_num']/varib['line_num'])

                ## Get Halo Data
		# Try and get halos from data_loc dictionary
		try:
			file = open(root+'/OSDCStacking/'+data_loc+'/halo_arrays.pkl','rb')
			input = pkl.Unpickler(file)
			data = input.load()
			D.add(data,keys=load_names)
			D.M_crit200,D.R_crit200,D.Z,D.HVD,D.HPX,D.HPY,D.HPZ,D.HVX,D.HVY,D.HVZ = D.HaloData
		except:			
                	# Load and Sort Halos by Mass
                	D.HaloID,D.HaloData = self.M.load_halos()
                	D.HaloID,D.HaloData = self.M.sort_halos(D.HaloID,D.HaloData)
                	D.HaloID,D.M_crit200,D.R_crit200,D.Z,D.HVD,D.HPX,D.HPY,D.HPZ,D.HVX,D.HVY,D.HVZ = np.vstack((D.HaloID,D.HaloData))
                	D.HaloID = np.array(D.HaloID,int)

                # Build Halo_P, Halo_V
                Halo_P = np.vstack([D.HPX,D.HPY,D.HPZ])
                Halo_V = np.vstack([D.HVX,D.HVY,D.HVZ])

		# Loop over ensembles
		j = 2
		for i in self.halo_range[1:]:
			# Progress Bar
			sys.stdout.write("Progress... "+str(j)+" out of "+str(len(self.halo_range))+"\r")
			sys.stdout.flush()
			j += 1
			pkl_file = open(root+'/OSDCStacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(i)+'_Data.pkl','rb')
			input = pkl.Unpickler(pkl_file)
			stack_data = input.load()
		
			# Append variables to data container D	
			D.append(stack_data,keys=load_names)

		print ''

		# Convert to arrays
		D.to_array(load_names)
		D.upper(names=load_names)

		# Return Data
		if raw_data == True:
			# Return either dictionary of dump to globals()
			if go_global == True:
				globals().update(D.__dict__)
				return
			elif go_global == False:
				return  D.__dict__

		################################
		### Statistical Calculations ###
		################################
		if self.ss:	# If self stacking is True
			M200 = D.M_crit200[0:self.halo_num]
			HVD200 = D.HVD[0:self.halo_num]	
			ENS_MFRAC,ens_mbias,ens_mscat,ENS_VFRAC,ens_vbias,ens_vscat = self.stat_calc(D.ENS_CAUMASS,D.M200,D.ENS_HVD,D.HVD200)
			if ens_only == True:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = None,None,None,None,None,None
			else:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = self.stat_calc(D.LOS_CAUMASS.ravel(),D.M_crit200[0:self.halo_num],D.LOS_HVD.ravel(),D.HVD[0:self.halo_num],ens=False)
		else:		# If self stacking is False
			if mm == True:
				ENS_MFRAC,ens_mbias,ens_mscat,ENS_VFRAC,ens_vbias,ens_vscat = self.stat_calc(D.ENS_CAUMASS,D.BINM200,D.ENS_HVD,D.BINHVD,data_set=None)#'cut_low_mass')
			else:	
				ENS_MFRAC,ens_mbias,ens_mscat,ENS_VFRAC,ens_vbias,ens_vscat = self.stat_calc(D.ENS_CAUMASS,D.BINM200,D.ENS_HVD,D.BINHVD)	
			if ens_only == True:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = None,None,None,None,None,None
			else:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = self.stat_calc(D.LOS_CAUMASS.ravel(),D.M_crit200[0:self.halo_num],D.LOS_HVD.ravel(),D.HVD[0:self.halo_num],ens=False)
	
		if cent_offset != None:
			OFF_MFRAC,off_ens_mbias,off_ens_mscat,OFF_ENS_VFRAC,off_ens_vbias,off_ens_vscat = self.stat_calc(D.OFF_ENS_CAUMASS,D.BIN_M200,D.OFF_ENS_HVD,D.BIN_HVD,data_set=None)


		# Create a dictionary
		names = ['ens_mbias','ens_mscat','los_mbias','los_mscat','ens_vbias','ens_vscat','los_vbias','los_vscat','ENS_MFRAC','ENS_VFRAC','LOS_MFRAC','LOS_VFRAC']
		mydict = ez.create(names,locals())
		D.add(mydict)

		# Return either a dictionary or dump to globals()
		if go_global == True:
			globals().update(D.__dict__)
			return
		elif go_global == False:
			return  D.__dict__



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


	def app(self,key,dictionary,LIST):
		'''
		This function sees if the string key is one of the keys of dictionary
		and if so, appends the dictionary[key] to LIST
		'''
		if key in dictionary.keys():
			LIST.append(dictionary[key])


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
		if kwargs['cent_offset'] != None:
			self.OFF_ENS_MBIAS.append(d.off_ens_mbias)
			self.OFF_ENS_MSCAT.append(d.off_ens_mscat)
			self.OFF_ENS_VBIAS.append(d.off_ens_vbias)
			self.OFF_ENS_VSCAT.append(d.off_ens_vscat)


	def load_all(self,iter_array=None,tab_shape=None,ens_only=True,kwargs=None,write_stem=None):
		'''
		This iterates over different richness geometry configurations and runs statistics on data.
		It is recommended to do any calculations (statistics, plots etc.) within the for loop, and
		then feed results back out via a global variable.
		'''
		# Feed Local Variables
		self.ens_only = ens_only	
		if kwargs == None:
			kwargs = {'write_loc':'mm_m0_run1','raw_data':False,'ss':False,'go_global':False,'ens_only':True,'data_loc':'mass_mix/mm_0.05_run_table1','cent_offset':None}
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
			if kwargs['cent_offset'] != None:
				self.OFF_ENS_MBIAS,self.OFF_ENS_MSCAT,self.OFF_ENS_VBIAS,self.OFF_ENS_VSCAT = [],[],[],[]

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
			self.ENS_MBIAS,self.ENS_MSCAT,self.ENS_VBIAS,self.ENS_VSCAT = np.array(self.ENS_MBIAS).reshape(tab_shape),np.array(self.ENS_MSCAT).reshape(tab_shape),np.array(self.ENS_VBIAS).reshape(tab_shape),np.array(self.ENS_VSCAT).reshape(tab_shape)
			self.LOS_MBIAS,self.LOS_MSCAT,self.LOS_VBIAS,self.LOS_VSCAT = np.array(self.LOS_MBIAS).reshape(tab_shape),np.array(self.LOS_MSCAT).reshape(tab_shape),np.array(self.LOS_VBIAS).reshape(tab_shape),np.array(self.LOS_VSCAT).reshape(tab_shape)
			self.RUN_NUM,self.GAL_NUM,self.LINE_NUM = np.array(self.RUN_NUM).reshape(tab_shape),np.array(self.GAL_NUM).reshape(tab_shape),np.array(self.LINE_NUM).reshape(tab_shape) 

			if kwargs['cent_offset'] != None:
				self.OFF_ENS_MBIAS,self.OFF_ENS_MSCAT,self.OFF_ENS_VBIAS,self.OFF_ENS_VSCAT = np.array(self.OFF_ENS_MBIAS).reshape(tab_shape),np.array(self.OFF_ENS_MSCAT).reshape(tab_shape),np.array(self.OFF_ENS_VBIAS).reshape(tab_shape),np.array(self.OFF_ENS_VSCAT).reshape(tab_shape)

			# Other Data Arrays
			self.RICH_NUM = self.GAL_NUM*self.LINE_NUM

			# Create Dictionary
			names = ['ENS_MBIAS','ENS_MSCAT','ENS_VBIAS','ENS_VSCAT','LOS_MBIAS','LOS_VBIAS','LOS_VSCAT','RUN_NUM','GAL_NUM','LINE_NUM','RICH_NUM','OFF_ENS_MBIAS','OFF_ENS_MSCAT','OFF_ENS_VBIAS','OFF_ENS_VSCAT']
			dictionary = ez.create(names,self.__dict__)
	
			return dictionary


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
	iter_array	= np.arange(1,50)
	tab_shape	= (7,7)
	write_stem	= 'bs_m0_run'
	write_loc	= 'bs_m0_run1'
	data_loc	= 'binstack/bs_run_table'+str(table_num)
	ss		= False
	mm		= False
	cent_offset	= None

	kwargs = {'write_loc':write_loc,'raw_data':False,'ss':ss,'mm':mm,'go_global':False,'ens_only':True,'data_loc':data_loc,'cent_offset':cent_offset}

	data = W.load_all(kwargs=kwargs,write_stem=write_stem,iter_array=iter_array,tab_shape=tab_shape)

	names = ['ENS_MBIAS','ENS_MSCAT','ENS_VBIAS','ENS_VSCAT','LOS_MBIAS','LOS_MSCAT','LOS_VBIAS','LOS_VSCAT','RUN_NUM','GAL_NUM','LINE_NUM','RICH_NUM','OFF_ENS_MBIAS','OFF_ENS_MSCAT','OFF_ENS_VBIAS','OFF_ENS_VSCAT']

	dictionary = ez.create(names,data)

	file = open(root+'/OSDCStacking/'+data_loc+'/table_analysis.pkl','wb')

	output = pkl.Pickler(file)

	output.dump(dictionary)

















