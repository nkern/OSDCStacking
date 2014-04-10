#caustic_universal_stack2D.py
'''
This script contains the class universal, which is used by caustic_class_stack2D.py and caustic_mass_stack2D.py

universal:
	functions:
	attributes:
'''

## Import Modules ##
import CausticMass as cm
import numpy as np
import astropy.io.fits as fits
import random
from scipy import weave
from scipy.weave import converters
import cosmolopy.distance as cd
import time
from numpy import random
from scipy.linalg import norm
import cPickle as pkl

## Program ##

class universal:

	def __init__(self,varib):
		''' Adding permanent program variables to class namespace'''
		self.__dict__.update(varib)

	def load_halos(self):
		'''This function loads halo data and makes cosmology corrections'''
		if self.small_set == True:
			HaloID = np.loadtxt(self.root+'/nkern/Caustic/biglosclusters.csv', delimiter=',', dtype='string', usecols=(0,), unpack=True)
			HPX,HPY,HPZ,HVX,HVY,HVZ = np.loadtxt(self.root+'/nkern/Caustic/biglosclusters.csv', delimiter=',', dtype='float', usecols=(9,10,11,12,13,14), unpack=True)
			SRAD,ESRAD,R_crit200,M_crit200,HVD,Z = np.loadtxt(self.root+'/nkern/Caustic/Millbig_concentrations.phys_phys.csv', delimiter=',', dtype='float', usecols=(1,2,5,7,9,12), unpack=True)
		else:
			HaloID,HPX,HPY,HPZ,HVX,HVY,HVZ,R_crit200,M_crit200,HVD,Z = np.loadtxt(self.root+'/nkern/Millennium/Large_Halo_Set/halos.csv',usecols=(0,8,9,10,11,12,13,16,5,7,4),delimiter=',',unpack=True)
			M_crit200 *= 1e10
			HaloID = np.array(HaloID,int)
			SRAD,ESRAD=np.ones(HaloID.size),np.ones(HaloID.size)

		# Hubble Constant Coefficient
		R_crit200,M_crit200,HPX,HPY,HPZ = R_crit200,M_crit200,HPX,HPY,HPZ

		# Cosmological Correction
		for l in xrange(len(HaloID)):	
			HPX[l],HPY[l],HPZ[l] = HPX[l]/(1+Z[l]),HPY[l]/(1+Z[l]),HPZ[l]/(1+Z[l])	
			# Fix weird SRAD values, if R_crit200/SRAD = Conc > 2, set SRAD=R_crit200/2
			if R_crit200[l]/SRAD[l] < 2.0:
				SRAD[l] = R_crit200[l] / 2.0
		#Ordering of halos in arrays is identical to biglosclusters' inherent ordering.
			
		return HaloID, np.vstack([M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ])


	def sort_halos(self,HaloID,HaloData):
		''' Sort Halo Data by some Criteria '''
		# Unpack Array HaloData into local namespace for easier use and clarity
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData	
		# Sort Arrays by descending M_crit200	
		sort = np.argsort(M_crit200)[::-1]	
		HaloID = HaloID[sort]
		M_crit200 = M_crit200[sort]
		R_crit200 = R_crit200[sort]
		Z = Z[sort]

		SRAD = SRAD[sort]
		ESRAD = ESRAD[sort]
		HVD = HVD[sort]
		HPX = HPX[sort]
		HPY = HPY[sort]
		HPZ = HPZ[sort]
		HVX = HVX[sort]
		HVY = HVY[sort]
		HVZ = HVZ[sort]
		# Return packed array
		return HaloID, np.vstack([M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ]) 

	def configure_galaxies(self,HaloID,HaloData):
		''' Loads galaxy data from halo list, and converts to physical coordinates and corrects cosmological factors '''
		# Unpack Array HaloData into local namespace for easier use and clarity
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData	

		G_Mags = []
		R_Mags = []
		I_Mags = []
		Gal_V = []
		Gal_P = []
		for k in self.stack_range:
			galdata = self.load_galaxies(HaloID[k],HaloData.T[k],R_crit200[k])
			# unpack array galdata into namespace
			gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags = galdata	
			gal_p	= np.array([gpx,gpy,gpz],float)
			gal_v	= np.array([gvx,gvy,gvz],float)
		
			G_Mags.append(gmags)
			R_Mags.append(rmags)	
			I_Mags.append(imags)
			Gal_P.append(gal_p)
			Gal_V.append(gal_v)			

		Halo_P,Halo_V = np.vstack([HPX,HPY,HPZ]).T,np.vstack([HVX,HVY,HVZ]).T
		
		return Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,HaloData[0:6]

	def load_galaxies(self,haloid,halodata,r_crit200):
		''' Loads haloid galaxies from a local directory '''
		# Unpack array halodata into local namespace
		m_crit200,r_crit200,z,srad,esrad,hvd,hpx,hpy,hpz,hvx,hvy,hvz = halodata
		# load galaxy data
		if self.small_set == True:
			# 100 Halo Sample	
			f = fits.open(self.root+'/giffordw/Millenium/30Mpchalos/'+haloid+'.'+self.data_set+'.fits')
			data = f[1].data
			gal_z,gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags = data.field(13),data.field(17),data.field(18),data.field(19),data.field(20),data.field(21),data.field(22),data.field(62),data.field(63),data.field(64)
		else:
			# 2,000 Halo Sample
			data = fits.getdata(self.root+'/nkern/Millennium/Large_Halo_Set/Halo_'+str(haloid)+'.Guo2010.fits')

			gal_z,gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags = data.field(3),data.field(6),data.field(7),data.field(8),data.field(9),data.field(10),data.field(11),data.field(14),data.field(15),data.field(16)

		# Cosmology corrections
		gpx,gpy,gpz = (gpx/(1+z)),(gpy/(1+z)),(gpz/(1+z))

		# Turn into array
		gmags,rmags,imags = np.array(gmags,float),np.array(rmags,float),np.array(imags,float)
		# remove BCG from sample
		BCG = np.where(gpx != hpx)

		return np.vstack([ gpx[BCG], gpy[BCG], gpz[BCG], gvx[BCG], gvy[BCG], gvz[BCG], gmags[BCG], rmags[BCG], imags[BCG] ])


	def scale_gals(self,r,r_crit200):
		''' Scales galaxy projected radius by a given radius'''
		r /= r_crit200
		return r

	def rand_pos(self,distance):
	        '''Picks a random position for the observer a given distance away from the center'''
		theta = random.normal(np.pi/2,np.pi/4)
		phi = random.uniform(0,2*np.pi)
		x = np.sin(theta)*np.cos(phi)
		y = np.sin(theta)*np.sin(phi)
		z = np.cos(theta)
	
		unit = np.array([x,y,z])/(x**2+y**2+z**2)**(.5)
		# move the position a random 'distance' Mpc away
	        return distance*unit


	def def_pos(self,distance,k):
		major_axis = True		# Project along major or minor axes?

		# Load Ellipticity Data for 100 Halos
		pkl_file = open(self.root+'/nkern/Stacking/Halo_Shape/100_halo_ellipticities.pkl','rb')
		input = pkl.Unpickler(pkl_file)
		d = input.load()
		eig_vec = d['eig_vec']
		eig_val = d['eig_val']

		# Define theta and phi for unit vector and construt unit vector in cartesian coordinates
		if major_axis == True:
			eig_vec = eig_vec[k][0]
			r = norm(eig_vec)
			theta = np.arccos(eig_vec[2]/r)
			phi = np.arctan(eig_vec[1]/eig_vec[0])

		if major_axis == False:
			eig_vec = eig_vec[k][1]
			r = norm(eig_vec)
			theta = np.arccos(eig_vec[2]/r)
			phi = np.arctan(eig_vec[1]/eig_vec[0])

	#	if random.random()<.5:
	#		eig_vec *= -1.

		theta = random.normal(theta,.075)
		phi = random.normal(phi,.075)
		x = np.sin(theta)*np.cos(phi)
		y = np.sin(theta)*np.sin(phi)
		z = np.cos(theta)	

		unit = np.array([x,y,z])/norm(np.array([x,y,z]))
		return unit*distance	




	def line_of_sight(self,gal_p,gal_v,halo_p,halo_v,k):
		'''Line of Sight Calculations to mock projected data, if given 3D data'''
		# Pick Position
		new_pos = self.rand_pos(30)
		new_pos += halo_p 

		# New Halo Information
		halo_dist = ((halo_p[0]-new_pos[0])**2 + (halo_p[1]-new_pos[1])**2 + (halo_p[2]-new_pos[2])**2)**0.5
		halo_pos_unit = np.array([halo_p[0]-new_pos[0],halo_p[1]-new_pos[1],halo_p[2]-new_pos[2]]) / halo_dist
		halo_vlos = np.dot(halo_pos_unit, halo_v)

		# New Galaxy Information
		gal_p = np.array(gal_p)
		gal_v = np.array(gal_v)
		gal_dist = ((gal_p[0]-new_pos[0])**2 + (gal_p[1]-new_pos[1])**2 + (gal_p[2]-new_pos[2])**2)**0.5
		gal_vlos = np.zeros(gal_dist.size)
		gal_pos_unit = np.zeros((3,gal_dist.size))	#vector from new_p to gal	
		n = gal_dist.size
		# Line of sight
		code = """
		int u,w;
		for (u=0;u<n;++u){
		for(w=0;w<3;++w){
		gal_pos_unit(w,u) = (gal_p(w,u)-new_pos(w))/gal_dist(u);
		}
		gal_vlos(u) = gal_pos_unit(0,u)*gal_v(0,u)+gal_pos_unit(1,u)*gal_v(1,u)+gal_pos_unit(2,u)*gal_v(2,u);
		}
		"""
		fast = weave.inline(code,['gal_pos_unit','n','gal_dist','gal_vlos','gal_v','new_pos','gal_p'],type_converters=converters.blitz,compiler='gcc')
		angles = np.arccos(np.dot(halo_pos_unit,gal_pos_unit))
		r = angles*halo_dist
		#v_pec = gal_vlos-halo_vlos*np.dot(halo_pos_unit,gal_pos_unit)
		z_clus_cos = self.H0*halo_dist/self.c
		z_clus_pec = halo_vlos/self.c
		z_clus_obs = (1+z_clus_pec)*(1+z_clus_cos)-1
		z_gal_cos = self.H0*gal_dist/self.c
		z_gal_pec = gal_vlos/self.c
		z_gal_obs = (1+z_gal_pec)*(1+z_gal_cos)-1
		v = self.c*(z_gal_obs-z_clus_obs)/(1+z_clus_obs)
		#gal_vdisp3d[i] = np.sqrt(astStats.biweightScale(gal_v[0][np.where(gal_radius<=HaloR200[i])]-Halo_V[0],9.0)**2+astStats.biweightScale(gal_v[1][np.where(gal_radius<=HaloR200[i])]-Halo_V[1],9.0)**2+astStats.biweightScale(gal_v[2][np.where(gal_radius<=HaloR200[i])]-Halo_V[2],9.0)**2)/np.sqrt(3)
		#print 'MY VELOCITY OF GALAXIES', gal_vdisp3d[i]
#		particle_vdisp3d[i] = HVD*np.sqrt(3)
#		gal_rmag_new = gal_abs_rmag# + 5*np.log10(gal_dist*1e6/10.0)

		return r, v, new_pos

	def limit_gals(self,r,v,en_gal_id,en_clus_id,ln_gal_id,gmags,rmags,imags,r200,hvd):
		''' Sort data by magnitude, and elimite values outside phase space limits '''
		# Sort by ascending r magnitude (bright to dim)
		sorts = np.argsort(rmags)
		r,v,en_gal_id,en_clus_id,ln_gal_id,gmags,rmags,imags = r[sorts],v[sorts],en_gal_id[sorts],en_clus_id[sorts],ln_gal_id[sorts],gmags[sorts],rmags[sorts],imags[sorts]

		# Limit Phase Space
		sample = np.where( (r < r200*self.r_limit) & (v > -self.v_limit) & (v < self.v_limit) )[0] 
		r,v,en_gal_id,en_clus_id,ln_gal_id,gmags,rmags,imags = r[sample],v[sample],en_gal_id[sample],en_clus_id[sample],ln_gal_id[sample],gmags[sample],rmags[sample],imags[sample]
		samp_size = len(sample)


		# Eliminate galaxies w/ mag = 99.
		cut = np.where((gmags!=99)&(rmags!=99)&(imags!=99))[0]
		r,v,en_gal_id,en_clus_id,ln_gal_id,gmags,rmags,imags = r[cut],v[cut],en_gal_id[cut],en_clus_id[cut],ln_gal_id[cut],gmags[cut],rmags[cut],imags[cut]
		samp_size = len(cut)
	
		return r,v,en_gal_id,en_clus_id,ln_gal_id,gmags,rmags,imags,samp_size


	def Bin_Calc(self,HaloData,varib):
		'''
		This function does pre-technique binning analysis
		'''
		# Unpack Arrays
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData

		# Sort Arrays based on specific Binning Method

		# Calculate Bin R200 and Bin HVD, use median
		BIN_M200,BIN_R200,BIN_HVD = [],[],[]
		for i in range(varib['halo_num']/varib['line_num']):
			BIN_M200.append( np.median( M_crit200[i*varib['line_num']:(i+1)*varib['line_num']] ) )
			BIN_R200.append( np.median( R_crit200[i*varib['line_num']:(i+1)*varib['line_num']] ) )
			BIN_HVD.append( np.median( HVD[i*varib['line_num']:(i+1)*varib['line_num']] ) )
	
		BIN_M200,BIN_R200,BIN_HVD = np.array(BIN_M200),np.array(BIN_R200),np.array(BIN_HVD)

		# Re-pack arrays
		BinData = np.vstack([BIN_M200,BIN_R200,BIN_HVD])

		return BinData


	def mass_mixing(self,HaloID,HaloData,mass_scat):
		'''
		This function performs a mass mixing procedure with a given fractional scatter in assumed mass
		'''

		# Unpack Array
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData
		
		# Create Gaussian distribution about 1 with width mass_scat, length HaloID.size
		mass_mix = random.lognormal(0,mass_scat,len(HaloID))

		# Apply Mass Scatter
		M_crit200 *= mass_mix

		# Sort by Descending Mass
		sort = np.argsort(M_crit200)[::-1]
		M_crit200 = M_crit200[sort]
		R_crit200 = R_crit200[sort]
		Z = Z[sort]
		SRAD = SRAD[sort]
		ESRAD = ESRAD[sort]
		HVD = HVD[sort]
		HPX = HPX[sort]
		HPY = HPY[sort]
		HPZ = HPZ[sort]
		HVX = HVX[sort]
		HVY = HVY[sort]
		HVZ = HVZ[sort]
		HaloID = HaloID[sort]

		# Re-pack
		HaloData = np.array([M_crit200,R_crit200,Z,SRAD,ESRAD,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ])

		return HaloID,HaloData
		
 
	def id_match(self,ID1,ID2,data):
		'''
		This function takes two IDs (one-dim arrays) and matches them up by identical value, the conversion being matching ID2 to ID1, then applying the match to data.
		This means that ID2 and data match in primary key
		'''
		sort = np.zeros(len(ID2),int)
		for i in np.arange(len(sort)):
			sort[i] = np.where(ID2[i] == ID1)[0][0]
		
		ID2 = ID2[sort]
		data = data.T[sort].T

		return sort,ID2,data			


	def mag_get_3d(self,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,gmags,rmags,imags):
		'''
		This function retreives the 3D position and velocity data for finalized galaxies in phase spaces for a given Halo.
		It matches the galaxies using the magnitude as a key, b/c rarely are magnitudes degenerate, however, sometimes they are.
		Therefore, three different magnitudes are provided, given the fact that degeneracy on all 3 levels for 2 or more galaxies is extremely low.
		'''
		gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = [],[],[],[],[],[]
		# Match Galaxies by Magnitude
		select = []
		for i in xrange(len(gmags)):
			pick = np.where(G_Mags==gmags[i])[0]
			size = len(pick)
			if size > 1:
				pick = np.where(R_Mags==rmags[i])[0]
				size = len(pick)
				if size > 1:
					pick = np.where(I_Mags==imags[i])[0]
					size = len(pick)
					if size > 1:
						print 'degeneracy on all magnitudes!'
						break

			select.append(int(pick[0]))
		select = np.array(select)
		gpx3d,gpy3d,gpz3d = Gal_P[0][select],Gal_P[1][select],Gal_P[2][select]
		gvx3d,gvy3d,gvz3d = Gal_V[0][select],Gal_V[1][select],Gal_V[2][select]

		return gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d

	def get_3d(self,Gal_P,Gal_V,ens_gal_id,los_gal_id,stack_range,clus_num,self_stack,j):
		'''
		This function recovers the 3D positions and velocities of the galaxies in the ensemble and los phase space.
		'''

		if self_stack == True:
			pass

		else:	
			# Create fully concatenated arrays to draw ensemble data from
			[BIN_GPX3D,BIN_GPY3D,BIN_GPZ3D] = map(np.concatenate,Gal_P[j*self.line_num:(j+1)*self.line_num].T)
			[BIN_GVX3D,BIN_GVY3D,BIN_GVZ3D] = map(np.concatenate,Gal_V[j*self.line_num:(j+1)*self.line_num].T)

			# Recover ensemble 3D data	
			ens_gpx3d,ens_gpy3d,ens_gpz3d = BIN_GPX3D[ens_gal_id],BIN_GPY3D[ens_gal_id],BIN_GPZ3D[ens_gal_id]
			ens_gvx3d,ens_gvy3d,ens_gvz3d = BIN_GVX3D[ens_gal_id],BIN_GVY3D[ens_gal_id],BIN_GVZ3D[ens_gal_id]

			# Recover line_of_sight 3D data
			los_gpx3d,los_gpy3d,los_gpz3d = [],[],[]
			los_gvx3d,los_gvy3d,los_gvz3d = [],[],[]
			for i,k in zip( np.arange(j*self.line_num,(j+1)*self.line_num), np.arange(self.line_num) ):
				los_gpx3d.append(Gal_P[i][0][los_gal_id[k]])
				los_gpy3d.append(Gal_P[i][1][los_gal_id[k]])
				los_gpz3d.append(Gal_P[i][2][los_gal_id[k]])
				los_gvx3d.append(Gal_V[i][0][los_gal_id[k]])
				los_gvy3d.append(Gal_V[i][1][los_gal_id[k]])
				los_gvz3d.append(Gal_P[i][2][los_gal_id[k]])

			los_gpx3d,los_gpy3d,los_gpz3d = np.array(los_gpx3d),np.array(los_gpy3d),np.array(los_gpz3d)
			los_gvx3d,los_gvy3d,los_gvz3d = np.array(los_gvx3d),np.array(los_gvy3d),np.array(los_gvz3d)

		return np.array([ens_gpx3d,ens_gpy3d,ens_gpz3d]),np.array([ens_gvx3d,ens_gvy3d,ens_gvz3d]),np.array([los_gpx3d,los_gpy3d,los_gpz3d]),np.array([los_gvx3d,los_gvy3d,los_gvz3d])


	def print_varibs(self,varib):
		print "Start Time		=",time.asctime()
		print "run_num			=",varib['run_num']
		print "clus_num		=",varib['clus_num']
		print "gal_num			=",varib['gal_num']
		print "line_num		=",varib['line_num']
		print "halo_num		=",varib['halo_num']
		print "method_num		=",varib['method_num']
		print "data_loc		=",varib['data_loc']
		print "write_loc		=",varib['write_loc']
		print "data_set		=",varib['data_set']
		print "small_set		=",varib['small_set']
		print "self_stack		=",varib['self_stack']
		print "write_data		=",varib['write_data']
		print "run_los			=",varib['run_los']
		print "mass_mix		=",varib['mass_mix']
		print "mass_scat		=",varib['mass_scat']
		return	

	def print_separation(self,text,type=1):
		if type==1:
			print ''
			print '-'*60
			print str(text)
			print '-'*60
		elif type==2:
			print ''
			print str(text)
			print '-'*30
		return




