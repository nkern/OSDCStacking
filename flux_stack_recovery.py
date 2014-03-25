"""
This program loads previously written data files from flux batch job running on the stacking method.
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

## Flags ##
use_flux	= True				# Running on flux or sophie?
get_los		= True				# Upload Line of Sight Data as well?
data_loc	= 'binstack_run_table2'		# Parent directory where write_loc lives
write_loc	= 'ss_m1_run1'			# Which directory within data_loc to load ensembles from?

if use_flux == True:
	root = str('/nfs/christoq_ls')
else:
	root = str('/n/Christoq1')

## Constants ##
warnings.filterwarnings("module",message="Warning: converting a masked element to nan.")


## Functions ##
class Recover(universal):
	'''This class contains functions that recover serialized data'''

	def __init__(self):
		pass
	
	def recover(self,write_loc=write_loc,raw_data=False,ss=True,go_global=True,ens_only=True,data_loc=None):

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

		# Create them as lists
		ENS_CAUMASS,ENS_CAUMASS_EST,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUMASS_EST,LOS_CAUSURF,LOS_NFWSURF = [],[],[],[],[],[],[],[]
		ENS_R,ENS_V,ENS_GMAGS,ENS_RMAGS,ENS_IMAGS,LOS_R,LOS_V,LOS_GMAGS,LOS_RMAGS,LOS_IMAGS,ENS_HVD,LOS_HVD = [],[],[],[],[],[],[],[],[],[],[],[]
		SAMS,PRO_POS,GPX3D,GPY3D,GPZ3D,GVX3D,GVY3D,GVZ3D = [],[],[],[],[],[],[],[]

		# Initialization step 
		if data_loc==None:	
			if self.ss:	data_loc = 'selfstack_run_table'
			else:	data_loc = 'binstack_run_table2'

		pkl_file = open(root+'/nkern/Stacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(0)+'_Data.pkl','rb')
		input = pkl.Unpickler(pkl_file)
		stack_data 			= input.load()
		varib				= input.load()

		self.__dict__.update(varib)
		self.U = universal(varib)
		if self.ss:	self.halo_range = range(2124)
		else:	self.halo_range = range(varib['halo_num']/varib['line_num'])

		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags,ens_hvd,ens_caumass,ens_caumass_est,ens_causurf,ens_nfwsurf,los_r,los_v,los_gmags,los_rmags,los_imags,los_hvd,los_caumass,los_caumass_est,los_causurf,los_nfwsurf,x_range,sample_size,pro_pos,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = stack_data

		# Append stack_data to major lists
		ENS_R.append(ens_r)
		ENS_V.append(ens_v)
		ENS_GMAGS.append(ens_gmags)
		ENS_RMAGS.append(ens_rmags)
		ENS_IMAGS.append(ens_imags)
		ENS_HVD.append(float(ens_hvd))
		ENS_CAUMASS.append(float(ens_caumass))
		ENS_CAUMASS_EST.append(float(ens_caumass_est))
		ENS_CAUSURF.append(ens_causurf)
		ENS_NFWSURF.append(ens_nfwsurf)
		LOS_R.append(los_r)
		LOS_V.append(los_v)
		LOS_GMAGS.append(los_gmags)
		LOS_RMAGS.append(los_rmags)
		LOS_IMAGS.append(los_imags)
		LOS_HVD.append(los_hvd)
		LOS_CAUMASS.append(los_caumass)
		LOS_CAUMASS_EST.append(los_caumass_est)
		LOS_CAUSURF.append(los_causurf)
		LOS_NFWSURF.append(los_nfwsurf)
		SAMS.append(sample_size)
		PRO_POS.append(pro_pos)
		GPX3D.append(gpx3d)
		GPY3D.append(gpy3d)
		GPZ3D.append(gpz3d)
		GVX3D.append(gvx3d)
		GVY3D.append(gvy3d)
		GVZ3D.append(gvz3d)

		# Loop over ensembles
		j = 0
		for i in self.halo_range[1:]:
			# Progress Bar
			sys.stdout.write("Progress... "+str(j)+"\r")
			sys.stdout.flush()
			j += 1
			pkl_file = open(root+'/nkern/Stacking/'+data_loc+'/'+write_loc+'/Ensemble_'+str(i)+'_Data.pkl','rb')
			input = pkl.Unpickler(pkl_file)
			stack_data = input.load()

			# Unpack Variables
			ens_r,ens_v,ens_gmags,ens_rmags,ens_imags,ens_hvd,ens_caumass,ens_caumass_est,ens_causurf,ens_nfwsurf,los_r,los_v,los_gmags,los_rmags,los_imags,los_hvd,los_caumass,los_caumass_est,los_causurf,los_nfwsurf,x_range,sample_size,pro_pos,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = stack_data

			# Append stack_data to major lists
			ENS_R.append(ens_r)
			ENS_V.append(ens_v)
			ENS_GMAGS.append(ens_gmags)
			ENS_RMAGS.append(ens_rmags)
			ENS_IMAGS.append(ens_imags)
			ENS_HVD.append(float(ens_hvd))
			ENS_CAUMASS.append(float(ens_caumass))
			ENS_CAUMASS_EST.append(float(ens_caumass_est))
			ENS_CAUSURF.append(ens_causurf)
			ENS_NFWSURF.append(ens_nfwsurf)
			LOS_R.append(los_r)
			LOS_V.append(los_v)
			LOS_GMAGS.append(los_gmags)
			LOS_RMAGS.append(los_rmags)
			LOS_IMAGS.append(los_imags)
			LOS_HVD.append(los_hvd)
			LOS_CAUMASS.append(los_caumass)
			LOS_CAUMASS_EST.append(los_caumass_est)
			LOS_CAUSURF.append(los_causurf)
			LOS_NFWSURF.append(los_nfwsurf)
			SAMS.append(sample_size)
			PRO_POS.append(pro_pos)
			GPX3D.append(gpx3d)
			GPY3D.append(gpy3d)
			GPZ3D.append(gpz3d)
			GVX3D.append(gvx3d)
			GVY3D.append(gvy3d)
			GVZ3D.append(gvz3d)

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
		LOS_HVD = np.array(LOS_HVD)
		LOS_CAUMASS = np.array(LOS_CAUMASS)
		LOS_CAUMASS_EST = np.array(LOS_CAUMASS_EST)
		LOS_CAUSURF = np.array(LOS_CAUSURF)
		LOS_NFWSURF = np.array(LOS_NFWSURF)
		SAMS = np.array(SAMS)
		PRO_POS = np.array(PRO_POS)
		GPX3D = np.array(GPX3D)
		GPY3D = np.array(GPY3D)
		GPZ3D = np.array(GPZ3D)
		GVX3D = np.array(GVX3D)
		GVY3D = np.array(GVY3D)
		GVZ3D = np.array(GVZ3D)

		## Get Halo Data
		# Load and Sort Halos by Mass
		HaloID,HaloData = self.U.load_halos()
		HaloID,HaloData = self.U.sort_halos(HaloID,HaloData)
		HaloID,M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = np.vstack((HaloID,HaloData))
		HaloID = np.array(HaloID,int)
		# Build Halo_P, Halo_V
		Halo_P = np.vstack([HPX,HPY,HPZ])
		Halo_V = np.vstack([HVX,HVY,HVZ])

		# Make Bin Arrays if ss == False
		if self.ss == False:
			BinData = self.U.Bin_Calc(HaloData[0:6],varib)
			BIN_M200,BIN_R200,BIN_HVD = BinData

		if raw_data == True:
			# Return to Namespaces depending on go_global
			names = ['varib','HaloID','Halo_P','Halo_V','GPX3D','GPY3D','GPZ3D','GVX3D','GVY3D','GVZ3D','M_crit200','R_crit200','Z','SRAD','ESRAD','HVD','x_range','ENS_CAUMASS','ENS_CAUMASS_EST','ENS_CAUSURF','ENS_NFWSURF','LOS_CAUMASS','LOS_CAUMASS_EST','LOS_CAUSURF','LOS_NFWSURF','ENS_R','ENS_V','ENS_GMAGS','ENS_RMAGS','ENS_IMAGS','LOS_R','LOS_V','LOS_GMAGS','LOS_RMAGS','LOS_IMAGS','ENS_HVD','LOS_HVD','SAMS','PRO_POS']
			data = [varib,HaloID,Halo_P,Halo_V,GPX3D,GPY3D,GPZ3D,GVX3D,GVY3D,GVZ3D,M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,x_range,ENS_CAUMASS,ENS_CAUMASS_EST,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUMASS_EST,LOS_CAUSURF,LOS_NFWSURF,ENS_R,ENS_V,ENS_GMAGS,ENS_RMAGS,ENS_IMAGS,LOS_R,LOS_V,LOS_GMAGS,LOS_RMAGS,LOS_IMAGS,ENS_HVD,LOS_HVD,SAMS,PRO_POS]
			mydict = dict( zip(names,data) )

			# Append Bin Arrays if bin stack
			if self.ss == False:
				mydict.update({'BIN_M200':BIN_M200,'BIN_R200':BIN_R200,'BIN_HVD':BIN_HVD})	

			if go_global == True:
				globals().update(mydict)
				globals().update(varib)
				return
			elif go_global == False:
				return  mydict	

		################################
		### Statistical Calculations ###
		################################
		if self.ss:
			ENS_MFRAC,ens_mbias,ens_mscat,ENS_VFRAC,ens_vbias,ens_vscat = self.stat_calc(ENS_CAUMASS,M_crit200,ENS_HVD,HVD)
			if ens_only == True:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = None,None,None,None,None,None
			else:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = self.stat_calc(LOS_CAUMASS.ravel(),M_crit200[0:self.halo_num],LOS_HVD.ravel(),HVD[0:self.halo_num],ens=False)
		else:
			ENS_MFRAC,ens_mbias,ens_mscat,ENS_VFRAC,ens_vbias,ens_vscat = self.stat_calc(ENS_CAUMASS,BIN_M200,ENS_HVD,BIN_HVD)
			if ens_only == True:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = None,None,None,None,None,None
			else:
				LOS_MFRAC,los_mbias,los_mscat,LOS_VFRAC,los_vbias,los_vscat = self.stat_calc(LOS_CAUMASS.ravel(),M_crit200[0:self.halo_num],LOS_HVD.ravel(),HVD[0:self.halo_num],ens=False)


		# Return to Namespaces depending on go_global
		names = ['varib','HaloID','Halo_P','Halo_V','GPX3D','GPY3D','GPZ3D','GVX3D','GVY3D','GVZ3D','M_crit200','R_crit200','Z','SRAD','ESRAD','HVD','x_range','ENS_CAUMASS','ENS_CAUMASS_EST','ENS_CAUSURF','ENS_NFWSURF','LOS_CAUMASS','LOS_CAUMASS_EST','LOS_CAUSURF','LOS_NFWSURF','ENS_R','ENS_V','ENS_GMAGS','ENS_RMAGS','ENS_IMAGS','LOS_R','LOS_V','LOS_GMAGS','LOS_RMAGS','LOS_IMAGS','ENS_HVD','LOS_HVD','SAMS','PRO_POS','ens_mbias','ens_mscat','los_mbias','los_mscat','ens_vbias','ens_vscat','los_vbias','los_vscat','ENS_MFRAC','ENS_VFRAC','LOS_MFRAC','LOS_VFRAC']
		data = [varib,HaloID,Halo_P,Halo_V,GPX3D,GPY3D,GPZ3D,GVX3D,GVY3D,GVZ3D,M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,x_range,ENS_CAUMASS,ENS_CAUMASS_EST,ENS_CAUSURF,ENS_NFWSURF,LOS_CAUMASS,LOS_CAUMASS_EST,LOS_CAUSURF,LOS_NFWSURF,ENS_R,ENS_V,ENS_GMAGS,ENS_RMAGS,ENS_IMAGS,LOS_R,LOS_V,LOS_GMAGS,LOS_RMAGS,LOS_IMAGS,ENS_HVD,LOS_HVD,SAMS,PRO_POS,ens_mbias,ens_mscat,los_mbias,los_mscat,ens_vbias,ens_vscat,los_vbias,los_vscat,ENS_MFRAC,ENS_VFRAC,LOS_MFRAC,LOS_VFRAC]
		mydict = dict( zip(names,data) )

		# Append Bin Arrays if bin stack
		if self.ss == False:
			mydict.update({'BIN_M200':BIN_M200,'BIN_R200':BIN_R200,'BIN_HVD':BIN_HVD})	

		if go_global == True:
			globals().update(mydict)
			globals().update(varib)
			return
		elif go_global == False:
			return  mydict	


	def stat_calc(self,MASS_EST,MASS_TRUE,HVD_EST,HVD_TRUE,data_set='full',ens=True):
		''' Does bias and scatter calculations '''
		# Define a Masked array for sometimes zero terms
		epsilon = 10.0
		use_est = False				# Use MassCalc estimated r200 mass values if true 
		maMASS_EST	= ma.masked_array(MASS_EST,mask=MASS_EST<epsilon)		# Mask essentially zero values
		maHVD_EST	= ma.masked_array(HVD_EST,mask=HVD_EST<epsilon)

		# 

		# Mass / HVD Fractions
		if ens == True:
			# Ensemble Arrays
			MFRAC = np.log(maMASS_EST/MASS_TRUE)
			VFRAC = np.log(maHVD_EST/HVD_TRUE)
		else:
			# LOS Mass Fraction Arrays: 0th axis is halo number, 1st axis is line of sight number
			MFRAC,VFRAC = [],[]
			for a in range(self.halo_num):
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
		self.GAL_NUM.append(d.varib['gal_num'])
		self.LINE_NUM.append(d.varib['line_num'])
		self.ENS_MBIAS.append(d.ens_mbias)
		self.ENS_MSCAT.append(d.ens_mscat)
		self.ENS_VBIAS.append(d.ens_vbias)
		self.ENS_VSCAT.append(d.ens_vscat)
		self.LOS_MBIAS.append(d.los_mbias)
		self.LOS_MSCAT.append(d.los_mscat)
		self.LOS_VBIAS.append(d.los_vbias)
		self.LOS_VSCAT.append(d.los_vscat)	


	def load_all(self,iter_array=None,tab_shape=None,ens_only=True):
		'''
		This iterates over different richness geometry configurations and runs statistics on data.
		It is recommended to do any calculations (statistics, plots etc.) within the for loop, and
		then feed results back out via a global variable.
		'''
		# Feed Local Variables
		self.ens_only = ens_only	
	
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
				kwargs = {'write_loc':'bs_m0_run'+str(i),'raw_data':False,'ss':False,'go_global':False,'ens_only':True,'data_loc':'binstack_run_table2'}
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

work = False
if work == True:
	data = W.load_all()

	names = ('ENS_MBIAS','ENS_MSCAT','ENS_VBIAS','ENS_VSCAT','LOS_MBIAS','LOS_MSCAT','LOS_VBIAS','LOS_VSCAT','RUN_NUM','GAL_NUM','LINE_NUM','RICH_NUM')

	values = np.copy(data)

	dictionary = dict(zip(names,values))

	file = open('binstack_run_table2_analysis.pkl','wb')

	output = pkl.Pickler(dictionary)

	output.dump(data)









