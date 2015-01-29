## stack_class.py
"""
This file contains classes and functions necessary in stacking Millennium Simulation
halos and running the caustic technique over them, utilizing both 'causticpy' and 'caustic_stack'
programs.
"""

from caustic_stack import *

class Millennium(object):
	"""
	A class with functions designed to work with Millennium Simulation data
	"""

	def __init__(self,varib):
		self.__dict__.update(varib)
		self.U = Universal(varib)


	def load_halos(self,new_coords='median'):
		'''This function loads halo data and makes cosmology corrections'''

		if self.lightcone == True:
			HaloID = np.loadtxt(self.root+'/OSDCStacking/centralsClusters2.csv',delimiter=',',usecols=(0,),unpack=True,dtype='str')
			RA,DEC,HPX,HPY,HPZ,HVX,HVY,HVZ,R_crit200,M_crit200,HVD,Z,Clus_rmag = np.loadtxt(self.root+'/OSDCStacking/centralsClusters2.csv',delimiter=',',usecols=(20,21,14,15,16,17,18,19,8,2,4,22,12),unpack=True)
			M_crit200 *= 1e10
			cut = []
			duplicate = []
			for i in range(len(HaloID)):
				x = np.where(HaloID == HaloID[i])[0]
				if len(x) == 1 and HaloID[i] not in duplicate:
					cut.append(i)
					continue
				elif len(x) > 1 and HaloID[i] not in duplicate:
					cut.append(i)
					duplicate.append(HaloID[i])
			cut = np.array(cut)
			HaloID = HaloID[cut]	
			RA,DEC,HPX,HPY,HPZ,HVX,HVY,HVZ,R_crit200,M_crit200,HVD,Z,Clus_rmag = RA[cut],DEC[cut],HPX[cut],HPY[cut],HPZ[cut],HVX[cut],HVY[cut],HVZ[cut],R_crit200[cut],M_crit200[cut],HVD[cut],Z[cut],Clus_rmag[cut]	

		else:
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
						raise NameError
					HALOID = np.array(HALOID,int)
					to_sort = np.array(map(lambda x: np.where(HALOID==x), HaloID)).ravel()
					HALOID,HPX,HPY,HPZ,HVX,HVY,HVZ = HALOID[to_sort],HPX[to_sort],HPY[to_sort],HPZ[to_sort],HVX[to_sort],HVY[to_sort],HVZ[to_sort]

		# Hubble Constant Coefficient
		# R_crit200,M_crit200,HPX,HPY,HPZ = R_crit200,M_crit200,HPX,HPY,HPZ

		# Cosmological Correction
		for l in xrange(len(HaloID)):	
			HPX[l],HPY[l],HPZ[l] = HPX[l]/(1+Z[l]),HPY[l]/(1+Z[l]),HPZ[l]/(1+Z[l])
			if self.lightcone == False:
				# Fix weird SRAD values, if R_crit200/SRAD = Conc > 2, set SRAD=R_crit200/2
				if R_crit200[l]/SRAD[l] < 2.0:
					SRAD[l] = R_crit200[l] / 2.0

		# Ordering of halos in arrays is identical to biglosclusters' inherent ordering.
		if self.lightcone == True:
			return HaloID, RA, DEC, np.vstack([M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ])
		else:
			return HaloID, np.vstack([M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ])


	def sort_halos(self,HaloID,HaloData,RA=None,DEC=None):
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
		if self.lightcone == True:
			RA = RA[sort]
			DEC = DEC[sort]

		# Return packed array
		if self.lightcone == True:
			return HaloID, RA, DEC, np.vstack([M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ]) 
		else:
			return HaloID, np.vstack([M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ]) 


	def configure_galaxies(self,haloid,halodata):
		''' Loads galaxy data from halo list, and converts to physical coordinates and corrects cosmological factors '''
		# Unpack Array HaloData into local namespace for easier use and clarity
		m_crit200,r_crit200,hvd,z,hpx,hpy,hpz,hvx,hvy,hvz = halodata

		galdata = self.load_galaxies(haloid,halodata)

		# unpack array galdata into namespace
		if self.lightcone == True:
			gal_ra,gal_dec,gal_z,gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags,abs_gmags,abs_rmags,abs_imags = galdata
		else:
			gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags = galdata	

		gal_p	= np.array([gpx,gpy,gpz],float)
		gal_v	= np.array([gvx,gvy,gvz],float)

		if self.lightcone == True:
			return gal_ra,gal_dec,gal_z,abs_gmags,abs_rmags,abs_imags,gal_p,gal_v
		else:	
			return gal_p,gal_v,gmags,rmags,imags


	def load_galaxies(self,haloid,halodata):
		''' Loads haloid galaxies from a local directory '''
		# Unpack array halodata into local namespace
		m_crit200,r_crit200,hvd,z,hpx,hpy,hpz,hvx,hvy,hvz = halodata

		# load galaxy data
		if self.lightcone == True:
			gal_id = np.loadtxt(self.root+'/Caustic/lowz_data2_2/'+str(haloid)+'.galaxies.tab',delimiter='\t',unpack=True,usecols=(0,),dtype='str')
			gal_ra,gal_dec,gal_z,gmags,rmags,imags,gpx,gpy,gpz,gvx,gvy,gvz,mem = np.loadtxt(self.root+'/Caustic/lowz_data2_2/'+str(haloid)+'.galaxies.tab',delimiter='\t',unpack=True,usecols=(1,2,3,5,6,7,9,10,11,12,13,14,15))
			abs_gmags,abs_rmags,abs_imags = np.loadtxt(self.root+'/Caustic/lowz_data2_2/'+str(haloid)+'.abs_rmags.tab',delimiter='\t',usecols=(2,3,4),unpack=True)
		
		else:
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
		BCG = np.where((gpx != hpx)&(gpy != hpy)&(gpz != hpz))[0]
		gpx, gpy, gpz, gvx, gvy, gvz, gmags, rmags, imags = gpx[BCG], gpy[BCG], gpz[BCG], gvx[BCG], gvy[BCG], gvz[BCG], gmags[BCG], rmags[BCG], imags[BCG]
		if self.lightcone == True:
			gal_ra, gal_dec, gal_z, abs_gmags, abs_rmags, abs_imags = gal_ra[BCG], gal_dec[BCG], gal_z[BCG], abs_gmags[BCG], abs_rmags[BCG], abs_imags[BCG]

		# Cut down to only members if desired
		if self.true_mems == True:
			gpr = np.sqrt((gpx-hpx)**2 + (gpy-hpy)**2 + (gpz-hpz)**2 )
			cut = np.where(gpr < r_crit200)
			gpx,gpy,gpz,gvx,gvy,gvz,gmags,rmags,imags = gpx[cut],gpy[cut],gpz[cut],gvx[cut],gvy[cut],gvz[cut],gmags[cut],rmags[cut],imags[cut]

		if self.lightcone == True:
			return np.vstack([ gal_ra, gal_dec, gal_z, gpx, gpy, gpz, gvx, gvy, gvz, gmags, rmags, imags, abs_gmags, abs_rmags, abs_imags])
		else:
			return np.vstack([ gpx, gpy, gpz, gvx, gvy, gvz, gmags, rmags, imags ])


	def load_project_append(self,haloid,m200,r200,hvd,clus_z,halo_p,halo_v,PS,clus_ra=None,clus_dec=None):
		"""
		This function loads galaxy data then projects it, discarding the galaxy data after using it to save memory
		"""

		if self.lightcone == True:
                        # Load Galaxies
			gal_ra,gal_dec,gal_z,gmags,rmags,imags,gal_p,gal_v = self.configure_galaxies(haloid,np.array([m200,r200,hvd,z,halo_p[0],halo_p[1],halo_p[2],halo_v[0],halo_v[1],halo_v[2]]))

			# Get angles and phase spaces
			ang_d,lum_d = S.C.zdistance(clus_z,self.H0)
			angles = S.C.findangle(gal_ra,gal_dec,clus_ra,clus_dec)
			rdata = angles * ang_d
			vdata = self.c * (gal_z - clus_z) / (1 + clus_z)
			pro_pos = [None]

			# Do Mass Mixing if Applicable	
			if self.mass_mix == True:
				if self.mm_est == 'richness':
					pass

				elif self.mm_est == 'vel_disp':
					pass

				elif self.mm_est == 'luminosity':
					pass

				else:
					print 'No match for mm_est'
					raise NameError

		else:
                        # Load Galaxies
			gal_p,gal_v,gmags,rmags,imags = self.configure_galaxies(haloid,np.array([m200,r200,hvd,z,halo_p[0],halo_p[1],halo_p[2],halo_v[0],halo_v[1],halo_v[2]]))

			# Do Projection, get phase spaces
			rdata, vdata, pro_pos = self.U.line_of_sight(gal_p,gal_v,halo_p,halo_v)
			rdata = np.array(rdata)
			vdata = np.array(vdata)
			pro_pos = np.array(pro_pos)

			# Do Mass Mixing if Applicable	
			if self.mass_mix == True:
				if self.mm_est == 'richness':
					pass

				elif self.mm_est == 'vel_disp':
					pass

				elif self.mm_est == 'luminosity':
					pass

				else:
					print 'No match for mm_est'
					raise NameError
		

		# Append to PS
		PS.append( {'Rdata':rdata,'Vdata':vdata,'pro_pos':np.array(pro_pos),'G_Mags':gmags,'R_Mags':rmags,'I_Mags':imags,'HaloID':haloid,'M200':m200,'R200':r200,'HVD':hvd} )


        def load_project_append_bootstrap(self,HaloID,M200,R200,HVD,Z,Halo_P,Halo_V,PS,weight):
                """
                This function loads galaxy data then projects it, discarding the galaxy data after using it to save memory
                """

                # Load Galaxies
		if self.lightcone == True:
			Gal_RA,Gal_DEC,Gal_Z,G_Mags,R_Mags,I_Mags,Gal_P,Gal_V = self.configure_galaxies(HaloID,np.array([M200,R200,HVD,Z,Halo_P[0],Halo_P[1],Halo_P[2],Halo_V[0],Halo_V[1],Halo_V[2]]))
		else:
                	Gal_P,Gal_V,G_Mags,R_Mags,I_Mags = self.configure_galaxies(HaloID,np.array([M200,R200,HVD,Z,Halo_P[0],Halo_P[1],Halo_P[2],Halo_V[0],Halo_V[1],Halo_V[2]]))

                # Do Projection
                r, v, pro_pos = self.U.line_of_sight(Gal_P,Gal_V,Halo_P,Halo_V,project=False)

                r = np.array(r)
                v = np.array(v)
                pro_pos = np.array(pro_pos)
                G_Mags = np.array(G_Mags)
                R_Mags = np.array(R_Mags)
                I_Mags = np.array(I_Mags)

		for m in range(weight):
                	# Append to PS
                	PS.append( {'Rdata':r,'Vdata':v,'pro_pos':np.array(pro_pos),'G_Mags':G_Mags,'R_Mags':R_Mags,'I_Mags':I_Mags,'HaloID':HaloID,'M200':M200,'R200':R200,'HVD':HVD} )


#	Outdated Mass Mixing Routine
#	def mass_mixing(self,HaloID,HaloData,mass_scat):
#		'''
#		This function performs a mass mixing procedure with a given fractional scatter in assumed mass
#		'''
#		# Unpack Array
#		M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ = HaloData
#		# Create lognormal distribution about 1 with width mass_scat, length HaloID.size
#		mass_mix = npr.lognormal(0,mass_scat,len(HaloID))
#		# Apply Mass Scatter
#		M_crit200 *= mass_mix
#		# Create M200_match array
#		M_crit200_match = np.copy(M_crit200)
#		# Sort by Descending Mass
#		sort = np.argsort(M_crit200)[::-1]
#		M_crit200 = M_crit200[sort]
#		R_crit200 = R_crit200[sort]
#		Z = Z[sort]
#		HVD = HVD[sort]
#		HPX = HPX[sort]
#		HPY = HPY[sort]
#		HPZ = HPZ[sort]
#		HVX = HVX[sort]
#		HVY = HVY[sort]
#		HVZ = HVZ[sort]
#		HaloID = HaloID[sort]
#		# Re-pack
#		HaloData = np.array([M_crit200,R_crit200,Z,HVD,HPX,HPY,HPZ,HVX,HVY,HVZ])
#		return HaloID,HaloData,sort


	def richness_mass_mixing(self,rdata,vdata,gmags,rmags,imags):
		''' Richness estimator combined with mass scatter calculation for millennium galaxy clusters'''
		# Rough Cut at > 3 Mpc and +/- 5000 km/s
		cut = np.where((rdata < 7) & (np.abs(vdata) < 5000))[0]
		rdata = rdata[cut]
		vdata = vdata[cut]
		gmags = gmags[cut]
		rmags = rmags[cut]
		imags = imags[cut]

		# Shiftgapper for Interlopers
		#clus_data = np.vstack([rdata,vdata,gmags,rmags,imags])
		#clus_data = self.C.shiftgapper(clus_data.T).T
		#rdata,vdata,gmags,rmags,imags = clus_data

		# Measure Velocity Dispersion of all galaxies within 3 Mpc
		vel_disp = astats.biweight_midvariance(vdata[np.where(rdata<3)])

		# Take rough virial radius measurement
		r_vir = np.exp(-1.86)*len(np.where((rmags < -19.55) & (rdata < 1.0) & (np.abs(vdata) < 3500))[0])**0.51

		# Find color of Red Sequence, measured as SDSS_g-SDSS_r vs. SDSS_r absolute magnitude
		color_data = gmags-rmags
		color_cut = np.where(color_data > np.mean(gmags-rmags))[0]
		size = len(color_cut)
		hist1 = mp.hist(color_data[color_cut],bins=size/25.,normed=True,histtype='step')
		hist2 = mp.hist(color_data[color_cut],bins=size/30.,normed=True,histtype='step')
		hist3 = mp.hist(color_data[color_cut],bins=size/35.,normed=True,histtype='step')
		mp.close()
		RS_color = np.mean([hist1[1][np.where(hist1[0]==hist1[0].max())][0],hist2[1][np.where(hist2[0]==hist2[0].max())][0],hist3[1][np.where(hist3[0]==hist3[0].max())][0]])
		RS_sigma = astats.biweight_midvariance(color_data)
	
		# Measure Sweet-Spot Richness
		# Sweet Spot Calculation when:
		# v_disp*2 < vdata < v_disp*3, r < r_virial, all colors, rmag_absolute < -19.0
		SS_richness = len(np.where((np.abs(vdata) < vel_disp*2)&(np.abs(vdata) > vel_disp*1)&(rdata <= r_vir)&(rmags < -19.0))[0])
		background = len(np.where((np.abs(vdata) < vel_disp*3)&(np.abs(vdata) > vel_disp*2)&(rdata <= r_vir*6)&(rdata >= r_vir*3)&(color_data<(RS_color+RS_sigma))&(color_data>(RS_color-RS_sigma))&(rmags<-19))[0])

		return SS_richness - background
	
	def vel_disp_mass_mixing(self,*args,**kwargs):
		''' Velocity Dispersion estimator combined with mass scatter calculation for millennium galaxy clusters '''
		pass

	def luminosity_mass_mixing(self,*args,**kwargs):
		''' Galaxy Luminosity estimator combined with mass scatter calculation for millennium galaxy clusters '''
		pass



























