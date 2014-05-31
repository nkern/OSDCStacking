## stack_class.py
"""
This file contains classes and functions necessary in stacking Millennium Simulation
halos and running the caustic technique over them, utilizing both 'causticpy' and 'caustic_stack'
programs.
"""

from caustic_stack import *
from analysis import *



class Millennium(object):
	"""
	A class with functions designed to work with Millennium Simulation data
	"""

	def __init__(self,varib):
		self.__dict__.update(varib)


	def load_halos(self,new_coords='median'):
		'''This function loads halo data and makes cosmology corrections'''
		if self.small_set == True:
			HaloID = np.loadtxt(self.root+'/Caustic/biglosclusters.csv', delimiter=',', dtype='string', usecols=(0,), unpack=True)
			HPX,HPY,HPZ,HVX,HVY,HVZ = np.loadtxt(self.root+'/Caustic/biglosclusters.csv', delimiter=',', dtype='float', usecols=(9,10,11,12,13,14), unpack=True)
			SRAD,ESRAD,R_crit200,M_crit200,HVD,Z = np.loadtxt(self.root+'/Caustic/Millbig_concentrations.phys_phys.csv', delimiter=',', dtype='float', usecols=(1,2,5,7,9,12), unpack=True)
		else:
			HaloID,HPX,HPY,HPZ,HVX,HVY,HVZ,R_crit200,M_crit200,HVD,Z = np.loadtxt(self.root+'/Millennium/Large_Halo_Set/halos.csv',usecols=(0,8,9,10,11,12,13,16,5,7,4),delimiter=',',unpack=True)
			M_crit200 *= 1e10
			HaloID = np.array(HaloID,int)
			SRAD,ESRAD=np.ones(HaloID.size),np.ones(HaloID.size)
			if self.new_halo_cent == True:
				if new_coords == 'mean':
					HALOID,HPX,HPY,HPZ,HVX,HVY,HVZ = np.loadtxt(self.root+'/OSDCStacking/new_halo_coords.csv',delimiter=',',usecols=(0,1,2,3,4,5,6),unpack=True)
				elif new_coords == 'median':
					HALOID,HPX,HPY,HPZ,HVX,HVY,HVZ = np.loadtxt(self.root+'/OSDCStacking/new_halo_coords.csv',delimiter=',',usecols=(0,7,8,9,10,11,12),unpack=True)
				else:
					print 'Old Coordinates Used...'
				HALOID = np.array(HALOID,int)
				to_sort = np.array(map(lambda x: np.where(HALOID==x), HaloID)).ravel()
				HALOID,HPX,HPY,HPZ,HVX,HVY,HVZ = HALOID[to_sort],HPX[to_sort],HPY[to_sort],HPZ[to_sort],HVX[to_sort],HVY[to_sort],HVZ[to_sort]

		# Hubble Constant Coefficient
		R_crit200,M_crit200,HPX,HPY,HPZ = R_crit200,M_crit200,HPX,HPY,HPZ

		# Cosmological Correction
		for l in xrange(len(HaloID)):	
			HPX[l],HPY[l],HPZ[l] = HPX[l]/(1+Z[l]),HPY[l]/(1+Z[l]),HPZ[l]/(1+Z[l])	
			# Fix weird SRAD values, if R_crit200/SRAD = Conc > 2, set SRAD=R_crit200/2
			if R_crit200[l]/SRAD[l] < 2.0:
				SRAD[l] = R_crit200[l] / 2.0
		#Ordering of halos in arrays is identical to biglosclusters' inherent ordering.
			
		return HaloID, np.vstack([M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ])


	def sort_halos(self,HaloID,HaloData):
		''' Sort Halo Data by some Criteria '''
		# Unpack Array HaloData into local namespace for easier use and clarity
		M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData	
		# Sort Arrays by descending M_crit200	
		sort = np.argsort(M_crit200)[::-1]	
		HaloID = HaloID[sort]
		M_crit200 = M_crit200[sort]
		R_crit200 = R_crit200[sort]
		Z = Z[sort]
		HVD = HVD[sort]
		HPX = HPX[sort]
		HPY = HPY[sort]
		HPZ = HPZ[sort]
		HVX = HVX[sort]
		HVY = HVY[sort]
		HVZ = HVZ[sort]
		# Return packed array
		return HaloID, np.vstack([M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ]) 


	def configure_galaxies(self,HaloID,HaloData):
		''' Loads galaxy data from halo list, and converts to physical coordinates and corrects cosmological factors '''
		# Unpack Array HaloData into local namespace for easier use and clarity
		M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData	

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
		
		return Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,HaloData[0:4]


	def load_galaxies(self,haloid,halodata,r_crit200):
		''' Loads haloid galaxies from a local directory '''
		# Unpack array halodata into local namespace
		m_crit200,r_crit200,z,hvd,hpx,hpy,hpz,hvx,hvy,hvz = halodata
		# load galaxy data
		if self.small_set == True:
			# 100 Halo Sample	
			f = fits.open(self.root+'/giffordw/Millenium/30Mpchalos/'+haloid+'.'+self.data_set+'.fits')
			data = f[1].data
			gal_z,gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags = data.field(13),data.field(17),data.field(18),data.field(19),data.field(20),data.field(21),data.field(22),data.field(62),data.field(63),data.field(64)
		else:
			# 2,000 Halo Sample
			data = fits.getdata(self.root+'/Millennium/Large_Halo_Set/Halo_'+str(haloid)+'.Guo2010.fits')

			gal_z,gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags = data.field(3),data.field(6),data.field(7),data.field(8),data.field(9),data.field(10),data.field(11),data.field(14),data.field(15),data.field(16)

		# Cosmology corrections
		gpx,gpy,gpz = (gpx/(1+z)),(gpy/(1+z)),(gpz/(1+z))

		# Turn into array
		gmags,rmags,imags = np.array(gmags,float),np.array(rmags,float),np.array(imags,float)

		# remove BCG from sample
		BCG = np.where((gpx != hpx)&(gpy != hpy)&(gpz != hpz))
		gpx, gpy, gpz, gvx, gvy, gvz, gmags, rmags, imags = gpx[BCG], gpy[BCG], gpz[BCG], gvx[BCG], gvy[BCG], gvz[BCG], gmags[BCG], rmags[BCG], imags[BCG]

		# Cut down to only members if desired
		if self.true_mems == True:
			gpr = np.sqrt((gpx-hpx)**2 + (gpy-hpy)**2 + (gpz-hpz)**2 )
			cut = np.where(gpr < r_crit200)
			gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags = gpx[cut],gpy[cut],gpz[cut],gvx[cut],gvy[cut],gvz[cut],gmags[cut],rmags[cut],imags[cut]

		return np.vstack([ gpx, gpy, gpz, gvx, gvy, gvz, gmags, rmags, imags ])









