## caustic_class_stack2D.py
'''
#########
- This script contains the class 'selfstack', used by caustic_mass_stack2D.py.
- Parts of the program are either modelled or directly taken from D. Gifford's CausticMass.py code. 
- This is the most up-to-date stacking code.
#########

self_stack:
	functions:
	attributes:
'''


## Import Modules ##
import numpy as np
import numpy.random as npr
import astStats

## Program ##

class Stack:

	def __init__(self,varib,U,C,CausticSurface,MassCalc):
		''' Initial function for class selfstack '''
		# Adding dictionary varib to class namespace
		self.__dict__.update(varib)
		# Initializing previously instanced classes 
		self.U = U		
		self.C = C
		self.CausticSurface = CausticSurface
		self.MassCalc = MassCalc


	def kernel_caustic_masscalc(self,r,v,halodata,bindata,derived_hvd,k,l=None):
		'''
		- This function runs the basic procedure of the Caustic Technique from CausticMass.py by using:
		- Caustic.gaussian_kernel()
		- CausticSurface.main()
		- MassCalc.main()
		'''

		## Unpack HaloData Array into Namespace
		# Define Cluster Variables Initially as Individual Halo Variables (Used for self stack and for LOS mass estimates)
		m_crit200,r_crit200,z,srad,esrad,hvd = halodata
		# If working on ensembles (not LOS) for bin stack, define r200 and hvd as bin r200 and hvd
		if self.self_stack == False and l == None:
			m_crit200,r_crit200,hvd = bindata

		## Print Details
		if l == None:
			text = '## Working On Cluster #'+str(k)+''
		else:
			text = '## Working On Cluster #'+str(k)+'\n## Line of Sight #'+str(l)+''
		self.U.print_separation(text,type=2)

		richness = len( np.where(r<r_crit200)[0])
		print 'Richness Within R200 =',richness
		print 'Calculated HVD =',np.floor(derived_hvd)

		## Mirror Phase Space Velocity Data if mirror == True
		mirror = True
		if mirror == True:
			rvalues,vvalues = np.append(r,r),np.append(v,-v)		
		else:
			rvalues,vvalues = r,v

		## Perform Kernel Density Estimation
		# Set Kernel Density Phase Space limits in Mpc and km/s
		xmax = 6.0			# Mpc
		ymax = 5000.0			# km/s
		self.C.gaussian_kernel(rvalues,vvalues,r_crit200,normalization=self.H0,scale=self.q,xmax=xmax,ymax=ymax,xres=200,yres=220)
		self.C.img_tot = self.C.img/np.max(np.abs(self.C.img))
 		self.C.img_grad_tot = self.C.img_grad/np.max(np.abs(self.C.img_grad))
 		self.C.img_inf_tot = self.C.img_inf/np.max(np.abs(self.C.img_inf))
		
		## Define Beta Profile
		beta_profile = False		# beta profile as a funct of radius?
		if beta_profile == True:
			xbeta,abeta = np.loadtxt(''+str(root)+'/nkern/Caustic/average_betaprofile.tab',dtype='float',usecols=(0,1),unpack=True)
			fit = np.polyfit((xbeta*r_crit200)[xbeta<4],abeta[xbeta<4],6)
			self.beta = fit[0]*self.C.x_range**6+fit[1]*self.C.x_range**5+fit[2]*self.C.x_range**4+fit[3]*self.C.x_range**3+fit[4]*self.C.x_range**2+fit[5]*self.C.x_range+fit[6]	
		else:				# constant beta w/ radius
			self.beta = np.ones(self.C.x_range.shape)*self.beta

		## Caustic Surface Calculation
		# This function takes RA, DEC and Z as first input, for now it is empty, all relavent files go to class dictionary
		if self.small_set == True:
			self.CS = self.CausticSurface(np.zeros(0),self.C.x_range,self.C.y_range,self.C.img_tot,r200=r_crit200,halo_scale_radius=srad,halo_scale_radius_e=esrad,halo_vdisp=derived_hvd,beta=self.beta)
		else:
			self.CS = self.CausticSurface(np.zeros(0),self.C.x_range,self.C.y_range,self.C.img_tot,r200=r_crit200,halo_vdisp=derived_hvd,beta=self.beta)

		## Mass Calculation, leave clus_z blank for now
		self.MC = self.MassCalc(self.C.x_range,self.CS.Ar_finalD,derived_hvd,0,r200=r_crit200,beta=self.beta,fbr=self.fbeta,H0=self.H0)

		return self.MC.M200,self.MC.M200_est,self.CS.Ar_finalD,self.CS.vesc_fit	


	def build_ensemble(self,r,v,gmags,rmags,imags,halodata,l):
		''' 
		- This function selects specific galaxies per line of sight using a sepcified method of stacking
		- The current method of interloper treatment is using CausticMass.py's ShiftGapper technique
		- Certain Measures were taken when using the ShiftGapper to ensure gal_num galaxies within r200
		- Interloper treatment is always done for the LOS, and can be done for ensemble_los if desired
		'''

		##
		# For ease
		gal_num = self.gal_num	
	
		## Program
		# Unpack halodata array into local namespace depending on SS or BS
		m_crit200,r_crit200,z,srad,esrad,hvd = halodata
	
		# Sort galaxies by r Magnitude
		bright = np.argsort(rmags)
		r,v,gmags,rmags,imags = r[bright],v[bright],gmags[bright],rmags[bright],imags[bright]

		if self.method_num == 0:
			en_r,en_v,en_gmags,en_rmags,en_imags,ln_r,ln_v,ln_gmags,ln_rmags,ln_imags = self.build_method_0(r,v,gmags,rmags,imags,r_crit200)

		elif self.method_num == 1:
			en_r,en_v,en_gmags,en_rmags,en_imags,ln_r,ln_v,ln_gmags,ln_rmags,ln_imags,samp_size = self.build_method_1(r,v,gmags,rmags,imags,r_crit200)
		
#		elif self.method_num == 2:
#			en_r,en_v,en_m,ln_r,ln_v,ln_m = self.build_method_2(r,v,mags,r_crit200)

		return en_r,en_v,en_gmags,en_rmags,en_imags,ln_r,ln_v,ln_gmags,ln_rmags,ln_imags



	def build_method_0(self,r,v,gmags,rmags,imags,r_crit200):
		'''Picking top brightest galaxies, such that there are gal_num galaxies within r200'''

		gal_num = self.gal_num

		## ------------------------------- ##
		## Build Ensemble (en => ensemble) ##
		## ------------------------------- ##
	
		# define indicies of galaxies within r200
		within = np.where(r<r_crit200)[0]
		# pick out gal_num 'th index in list, (include extra to counter shiftgapper's possible dimishement of richness)
		if gal_num < 10:
			excess = gal_num * 3.0 / 5.0				# galaxy excess to counter shiftgapper
		else:
			excess = gal_num / 5.0
		end = within[:gal_num + excess + 1][-1]		# instead of indexing I am slicing b/c of chance of not enough gals existing...	
		# Choose Ngals within r200
		if self.clean_ens == True:
			excess *= 2.0				# make excess a bit larger than previously defined
			end = within[:gal_num + excess + 1][-1]
			r2,v2,gmags2,rmags2,imags2 = self.C.shiftgapper(np.vstack([r[:end],v[:end],gmags[:end],rmags[:end],imags[:end]]).T).T # Shiftgapper inputs and outputs data as transpose...
			within = np.where(r2<r_crit200)[0]	# re-calculate within array with new sample
			excess = gal_num / 5.0
			end = within[:gal_num + excess + 1][-1]
			# Append to ensemble array
			en_r,en_v,en_gmags,en_rmags,en_imags = r2[:end],v2[:end],gmags2[:end],rmags2[:end],imags2[:end]
		else:
			en_r,en_v,en_gmags,en_rmags,en_imags = r[0:end],v[0:end],gmags[0:end],rmags[0:end],imags[0:end]


		## ----------------------------------------- ##
		## Build Line of Sight (ln => line of sight) ##
		## ----------------------------------------- ##

		# define indicies of galaxies within r200
		within = np.where(r<r_crit200)[0]
		# pick out gal_num 'th index in list, (include extra to counter shiftgapper's possible dimishement of richness)
		if gal_num < 10:
			excess = gal_num * 3.0 / 5.0				# galaxy excess to counter shiftgapper
		else:
			excess = gal_num / 5.0
		end = within[:gal_num + excess + 1][-1]		# instead of indexing I am slicing b/c
		# shiftgapper on line of sight
		r2,v2,gmags2,rmags2,imags2 = self.C.shiftgapper(np.vstack([r[:end],v[:end],gmags[:end],rmags[:end],imags[:end]]).T).T
		within = np.where(r2<r_crit200)[0]		# re-calculate within array with new sample
		# Sort by rmags
		sort = np.argsort(rmags2)
		r2,v2,gmags2,rmags2,imags2 = r2[sort],v2[sort],gmags2[sort],rmags2[sort],imags2[sort]
		# Now feed ln arrays correct gal_num richness within r200
		end = within[:gal_num + 1][-1]
		ln_r,ln_v,ln_gmags,ln_rmags,ln_imags = r2[:end],v2[:end],gmags2[:end],rmags2[:end],imags2[:end]	

		# Done! Now we have en_r and ln_r arrays, which will either be stacked (former) or put straight into Caustic technique (latter)
		return en_r,en_v,en_gmags,en_rmags,en_imags,ln_r,ln_v,ln_gmags,ln_rmags,ln_imags


	def build_method_1(self,r,v,gmags,rmags,imags,r_crit200):
		'''Randomly choosing bright galaxies until gal_num galaxies are within r200'''
		gal_num = self.gal_num
		
		# reduce size of sample to something reasonable within magnitude limits
		sample = gal_num * 25				# arbitrary coefficient, see sites page post Apr 24th, 2013 for more info
		r,v,gmags,rmags,imags = r[:sample],v[:sample],gmags[:sample],rmags[:sample],imags[:sample]
		samp_size = len(r)				# actual size of sample (might be less than gal_num*25)
		self.samp_size = samp_size

		# create random numbered array for galaxy selection
		if gal_num < 10:				# when gal_num < 10, has trouble hitting gal_num richness inside r200
			excess = gal_num * 4.0 / 5.0
		else:
			excess = gal_num * 2.0 / 5.0

		samp_num = gal_num + excess			# sample size of randomly generated numbers, start too low on purpose, then raise in loop
		loop = True					# break condition
		while loop == True:				# see method 0 comments on variables such as 'excess' and 'within' and 'end'
			for j in range(3):			# before upping sample size, try to get a good sample a few times
				rando = npr.randint(0,samp_size,samp_num)
				within = np.where(r[rando]<=r_crit200)[0]
				if len(within) >= gal_num + excess:
					loop = False
			if len(within) < gal_num + excess:
				samp_num += 2
		### Build Ensemble
		if self.clean_ens == True:

			r2,v2,gmags2,rmags2,imags2 = self.C.shiftgapper(np.vstack([r[rando],v[rando],gmags[rando],rmags[rand],imags[rando]]).T).T
			within = np.where(r2<r_crit200)[0]
			excess = gal_num / 5.0
			end = within[:gal_num + excess + 1][-1]
			# Append to ensemble array
			en_r,en_v,en_gmags,en_rmags,en_imags = r2[:end],v2[:end],gmags2[:end],rmags2[:end],imags2[:end]
		else:
			excess = gal_num / 5.0
			end = within[:gal_num + excess + 1][-1]
			en_r,en_v,en_gmags,en_rmags,en_imags = r[rando][:end],v[rando][:end],gmags[rando][:end],rmags[rando][:end],imags[rando][:end]

		### Build LOS
		if gal_num < 10:
			excess = gal_num * 4 / 5.0
		else:
			excess = gal_num / 5.0
		try:
			end = within[:gal_num + excess + 1][-1]
			# shiftgapper
			r2,v2,gmags2,rmags2,imags2 = self.C.shiftgapper(np.vstack([r[rando][:end],v[rando][:end],gmags[rando][:end],rmags[rando][:end],imags[rando][:end]]).T).T
			# sort by rmags
			sort = np.argsort(rmags2)
			r2,v2,gmags2,rmags2,imags2 = r2[sort],v2[sort],gmags2[sort],rmags2[sort],imags2[sort]
			# feed back gal_num gals within r200
			within = np.where(r2<r_crit200)[0]
			end = within[:gal_num + 1][-1]
			richness = len(within)
		except IndexError:
			print '****RAISED INDEX ERROR on LOS Building****'
			richness = 0		

		# Make sure gal_num richness inside r200
		j = 0
		while richness < gal_num:
			j += 1
			loop = True
			while loop == True:				
				for j in range(3):			
					rando = npr.randint(0,samp_size,samp_num)
					within = np.where(r[rando]<=r_crit200)[0]
					if len(within) >= gal_num + excess:
						loop = False
				if len(within) < gal_num + excess:
					samp_num += 2
			try:
				end = within[:gal_num + excess + 1][-1]
				r2,v2,gmags2,rmags2,imags2 = self.C.shiftgapper(np.vstack([r[rando][:end],v[rando][:end],gmags[rando][:end],rmags[rando][:end],imags[rando][:end]]).T).T
				within = np.where(r2<r_crit200)[0]
				end = within[:gal_num + 1][-1]
				richness = len(within)
			except IndexError:
				richness = 0

			if j >= 100:
				break

		ln_r,ln_v,ln_gmags,ln_rmags,ln_imags = r2[:end],v2[:end],gmags2[:end],rmags2[:end],imags2[:end]
		# Done! Now we have en_r and ln_r arrays (ensemble and line of sight arrays)
		
		return en_r,en_v,en_gmags,en_rmags,en_imags,ln_r,ln_v,ln_gmags,ln_rmags,ln_imags,samp_size


	def build_method_2(self,r,v,mags,r_crit200):
		'''Ordered Set of galaxies with respect to magnitude'''
		### NOT FINSHED YET !! ###
		### Build Ensemble 
		within = np.where(r<r_crit200)[0]
		excess = gal_num / 5.0
		end = gal_num + excess
		if clean_end_los == True:
			pass		

		end = within[gal_num*line_num + line_num]
		en_r,en_v,en_mags,en_gpx3d,en_gpy3d,en_gpz3d,en_gvx3d,en_gvy3d,en_gvz3d = r[:end][l::line_num],v[:end][l::line_num],mags[:end][l::line_num],gpx3d[:end][l::line_num],gpy3d[:end][l::line_num],gpz3d[:end][l::line_num],gvx3d[:end][l::line_num],gvy3d[:end][l::line_num],gvz3d[:end][l::line_num]	
		### Build LOS
		end = within[:gal_num*line_num + 5*line_num][-1] + 1
		rt,vt,magst,gpx3dt,gpy3dt,gpz3dt,gvx3dt,gvy3dt,gvz3dt = U.shiftgapper(vstack((r[:end][l::line_num],v[:end][l::line_num],mags[:end][l::line_num],gpx3d[:end][l::line_num],gpy3d[:end][l::line_num],gpz3d[:end][l::line_num],gvx3d[:end][l::line_num],gvy3d[:end][l::line_num],gvz3d[:end][l::line_num])).T).T
		within = where(rt<r_crit200)[0]
		end = within[:gal_num][-1] + 1
		r,v,mags,gpx3d,gpy3d,gpz3d,gvx3d,gvy3d,gvz3d = rt[:end],vt[:end],magst[:end],gpx3dt[:end],gpy3dt[:end],gpz3dt[:end],gvx3dt[:end],gvy3dt[:end],gvz3dt[:end] 
		###########################
		return en_r,en_v,en_m,ln_r,ln_v,ln_m





class SelfStack:

	def __init__(self,varib,U,C,S,CausticSurface,MassCalc):
		''' Initial function for class selfstack '''
		# Adding dictionary varib to class namespace
		self.__dict__.update(varib)
		# Initializing previously instanced classes 
		self.U = U		
		self.C = C
		self.S = S
		self.CausticSurface = CausticSurface
		self.MassCalc = MassCalc


	def self_stack_clusters(self,HaloID,HaloData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j):
		''' Building Ensemble Cluster and Calculating Property Statistics '''
		## Unpack HaloData array 
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData

		## Define Arrays for Building Ensemble and LOS
		# Ensemble Arrays:	[Successive Ensemble Number][Data]
		# Line of Sight Arrays:	[Line Of Sight][Data]
		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags,ens_hvd = [],[],[],[],[],[]
		ens_caumass, ens_caumass_est, ens_causurf, ens_nfwsurf = [], [], [], []
		los_r,los_v,los_gmags,los_rmags,los_imags,los_hvd = [],[],[],[],[],[]
		los_caumass, los_caumass_est, los_causurf, los_nfwsurf = [], [], [], []
		sample_size,pro_pos = [],[]

		## Loop over lines of sight
		for l in range(self.line_num):
			if self.light_cone == True:
				# Configure RA, DEC and Z into cluster-centric radius and velocity
				pass
			else:
				# Line of Sight Calculation for naturally 3D data
				r, v, projected_pos = self.U.line_of_sight(Gal_P[j],Gal_V[j],Halo_P[k],Halo_V[k],k)

			# Limit Data in Phase Space
			r,v,gmags,rmags,imags,samp_size = self.U.limit_gals(r,v,G_Mags[j],R_Mags[j],I_Mags[j],R_crit200[k],HVD[k])

			# Build LOS and Ensemble, with given method of stacking
			en_r,en_v,en_gmags,en_rmags,en_imags,ln_r,ln_v,ln_gmags,ln_rmags,ln_imags = self.S.build_ensemble(r,v,gmags,rmags,imags,HaloData.T[k],l)	

			# Build Ensemble Arrays
			ens_r.extend(en_r)
			ens_v.extend(en_v)
			ens_gmags.extend(en_gmags)
			ens_rmags.extend(en_rmags)
			ens_imags.extend(en_imags)
			
			# Calculate LOS HVD (this is after shiftgapper) if run_los == True
			if self.run_los == True:
				ln_within = np.where(ln_r<R_crit200[k])[0]
				gal_count = len(ln_within)
				if gal_count <= 3:
					'''biweightScale can't take less than 4 elements'''
					# Calculate hvd with numpy std of galaxies within r200 (b/c this is quoted richness)
					ln_hvd = np.std(np.copy(ln_v)[ln_within])
				else:
					# Calculate hvd with astStats biweightScale (see Beers 1990)
					try:
						ln_hvd = astStats.biweightScale(np.copy(ln_v)[ln_within],9.0)
					# Sometimes divide by zero error in biweight function for low gal_num
					except ZeroDivisionError:
						print 'ZeroDivisionError in biweightfunction, line 140 in caustic_class_stack2D'
						print 'ln_v[ln_within]=',ln_v[ln_within]
						ln_hvd = np.std(np.copy(ln_v)[ln_within])
			else:
				ln_hvd = []

			# Run Caustic Technique for LOS mass estimation if run_los == True
			if self.run_los == True:
				ln_caumass,ln_caumass_est,ln_causurf,ln_nfwsurf = self.S.kernel_caustic_masscalc(ln_r,ln_v,HaloData.T[k],np.zeros(2),ln_hvd,k,l)
				print 'j = '+str(j)+', k = '+str(k)+', l = '+str(l)
			else:
				ln_caumass,ln_caumass_est,ln_causurf,ln_nfwsurf = [],[],[],[]
	
			# Append LOS Data Arrays
			los_r.append(ln_r)
			los_v.append(ln_v)
			los_gmags.append(ln_gmags)
			los_rmags.append(ln_rmags)
			los_imags.append(ln_imags)
			los_hvd.append(ln_hvd)
			los_caumass.append(ln_caumass)
			los_caumass_est.append(ln_caumass_est)
			los_causurf.append(ln_causurf)
			los_nfwsurf.append(ln_nfwsurf)
			sample_size.append(samp_size)
			pro_pos.append(projected_pos)

		# Shiftgapper for Ensemble Interloper treatment
		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags = self.C.shiftgapper(np.vstack([ens_r,ens_v,ens_gmags,ens_rmags,ens_imags]).T).T

		# Sort by R_Mag
		sort = np.argsort(ens_rmags)
		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags = ens_r[sort],ens_v[sort],ens_gmags[sort],ens_rmags[sort],ens_imags[sort]

		# Reduce system to gal_num richness within r200
		within = np.where(ens_r <= R_crit200[k])[0]
		end = within[:self.gal_num*self.line_num + 1][-1]
		ens_r = ens_r[:end]
		ens_v = ens_v[:end]
		ens_gmags = ens_gmags[:end]
		ens_rmags = ens_rmags[:end]
		ens_imags = ens_imags[:end]

		# Calculate HVD
		en_hvd = astStats.biweightScale(np.copy(ens_v)[np.where(ens_r<=R_crit200[k])],9.0)

		# Run Caustic Technique! 
		en_caumass,en_caumass_est,en_causurf,en_nfwsurf = self.S.kernel_caustic_masscalc(ens_r,ens_v,HaloData.T[k],np.zeros(2),en_hvd,k)

		# Append Ensemble Data Arrays
		ens_hvd.append(en_hvd)
		ens_caumass.append(en_caumass)
		ens_caumass_est.append(en_caumass_est)

		# Turn into numpy arrays
		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags = np.array(ens_r),np.array(ens_v),np.array(ens_gmags),np.array(ens_rmags),np.array(ens_imags)
		ens_hvd,ens_caumass,ens_caumass_est = np.array(ens_hvd),np.array(ens_caumass),np.array(ens_caumass_est)
		ens_causurf,ens_nfwsurf = np.array(en_causurf),np.array(en_nfwsurf)
		los_r,los_v,los_gmags,los_rmags,los_imags = np.array(los_r),np.array(los_v),np.array(los_gmags),np.array(los_rmags),np.array(los_imags)
		los_hvd,los_caumass,los_caumass_est = np.array(los_hvd),np.array(los_caumass),np.array(los_caumass_est)
		los_causurf,los_nfwsurf = np.array(los_causurf),np.array(los_nfwsurf)
		sample_size,pro_pos = np.array(sample_size),np.array(pro_pos)

		return ens_r,ens_v,ens_gmags,ens_rmags,ens_imags,ens_hvd,ens_caumass,ens_caumass_est,ens_causurf,ens_nfwsurf,los_r,los_v,los_gmags,los_rmags,los_imags,los_hvd,los_caumass,los_caumass_est,los_causurf,los_nfwsurf,self.C.x_range,sample_size,pro_pos	





class BinStack:

	def __init__(self,varib,U,C,S,CausticSurface,MassCalc):
		''' Initial function for class BinStack '''
		# Adding dictionary varib to class namespace
		self.__dict__.update(varib)
		# Initializing previously instanced classes 
		self.U = U		
		self.C = C
		self.S = S
		self.CausticSurface = CausticSurface
		self.MassCalc = MassCalc



	def bin_stack_clusters(self,HaloID,HaloData,BinData,Halo_P,Halo_V,Gal_P,Gal_V,G_Mags,R_Mags,I_Mags,k,j):
		''' Building Ensemble Cluster and Calculating Property Statistics '''
		## Unpack HaloData array 
		M_crit200,R_crit200,Z,SRAD,ESRAD,HVD = HaloData
		BIN_M200,BIN_R200,BIN_HVD = BinData

		## Define Arrays for Building Ensemble and LOS
		# Ensemble Arrays:	[Successive Ensemble Number][Data]
		# Line of Sight Arrays:	[Line Of Sight][Data]
		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags,ens_hvd = [],[],[],[],[],[]
		ens_caumass, ens_caumass_est, ens_causurf, ens_nfwsurf = [], [], [], []
		los_r,los_v,los_gmags,los_rmags,los_imags,los_hvd = [],[],[],[],[],[]
		los_caumass, los_caumass_est, los_causurf, los_nfwsurf = [], [], [], []
		sample_size,pro_pos = [],[]

		## Loop over lines of sight (different clusters)
		for [l,s] in zip(np.arange(self.line_num*j,self.line_num*(j+1)),self.stack_range[j*self.line_num:(j+1)*self.line_num]):

			# Define index for to-be-stacked halo as cluster index (k) + line of sight index (l)

			if self.light_cone == True:
				# Configure RA, DEC and Z into cluster-centric radius and velocity
				pass
			else:
				# Line of Sight Calculation for naturally 3D data
				r, v, projected_pos = self.U.line_of_sight(Gal_P[l],Gal_V[l],Halo_P[s],Halo_V[s],s)

			# Limit Data in Phase Space
			r,v,gmags,rmags,imags,samp_size = self.U.limit_gals(r,v,G_Mags[l],R_Mags[l],I_Mags[l],R_crit200[s],HVD[s])

			# Build LOS and Ensemble, with given method of stacking
			en_r,en_v,en_gmags,en_rmags,en_imags,ln_r,ln_v,ln_gmags,ln_rmags,ln_imags = self.S.build_ensemble(r,v,gmags,rmags,imags,HaloData.T[s],l)	

			# If Scale data before stack is desired
			if self.scale_data == True:
				en_r = self.U.scale_gals(en_r,R_crit200[s])				
	
			# Build Ensemble Arrays
			ens_r.extend(en_r)
			ens_v.extend(en_v)
			ens_gmags.extend(en_gmags)
			ens_rmags.extend(en_rmags)
			ens_imags.extend(en_imags)
			
			# Calculate LOS HVD (this is after shiftgapper) if run_los == True
			if self.run_los == True:
				# Pick out gals within r200
				ln_within = np.where(ln_r<R_crit200[s])[0]
				gal_count = len(ln_within)
				if gal_count <= 3:
					'''biweightScale can't take less than 4 elements'''
					# Calculate hvd with numpy std of galaxies within r200 (b/c this is quoted richness)
					ln_hvd = np.std(np.copy(ln_v)[ln_within])
				else:
					# Calculate hvd with astStats biweightScale (see Beers 1990)
					try:
						ln_hvd = astStats.biweightScale(np.copy(ln_v)[ln_within],9.0)
					# Sometimes divide by zero error in biweight function for low gal_num
					except ZeroDivisionError:
						print 'ZeroDivisionError in biweightfunction, line 140 in caustic_class_stack2D'
						print 'ln_v[ln_within]=',ln_v[ln_within]
						ln_hvd = np.std(np.copy(ln_v)[ln_within])
			else:
				ln_hvd = []

			# Run Caustic Technique for LOS mass estimation if run_los == True
			if self.run_los == True:
				ln_caumass,ln_caumass_est,ln_causurf,ln_nfwsurf = self.S.kernel_caustic_masscalc(ln_r,ln_v,HaloData.T[s],BinData.T[k],ln_hvd,k,l)
				print 'j = '+str(j)+', k = '+str(k)+', l = '+str(l)+', s = '+str(s)
			else:
				ln_caumass,ln_caumass_est,ln_causurf,ln_nfwsurf = [],[],[],[]
		
	
			# Append LOS Data Arrays
			los_r.append(ln_r)
			los_v.append(ln_v)
			los_gmags.append(ln_gmags)
			los_rmags.append(ln_rmags)
			los_imags.append(ln_imags)
			los_hvd.append(ln_hvd)
			los_caumass.append(ln_caumass)
			los_caumass_est.append(ln_caumass_est)
			los_causurf.append(ln_causurf)
			los_nfwsurf.append(ln_nfwsurf)
			sample_size.append(samp_size)
			pro_pos.append(projected_pos)

		# If scale data == True, re-scale by ensemble r200
		if self.scale_data == True:
			ens_r = np.array(ens_r)*BIN_R200[k]

		# Shiftgapper for Ensemble Interloper treatment
		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags = self.C.shiftgapper(np.vstack([ens_r,ens_v,ens_gmags,ens_rmags,ens_imags]).T).T

		# Sort by R_Mag
		sort = np.argsort(ens_rmags)
		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags = ens_r[sort],ens_v[sort],ens_gmags[sort],ens_rmags[sort],ens_imags[sort]

		# Reduce system to gal_num richness within r200
		within = np.where(ens_r <= BIN_R200[k])[0]
		end = within[:self.gal_num*self.line_num + 1][-1]
		ens_r = ens_r[:end]
		ens_v = ens_v[:end]
		ens_gmags = ens_gmags[:end]
		ens_rmags = ens_rmags[:end]
		ens_imags = ens_imags[:end]

		# Calculate HVD
		en_hvd = astStats.biweightScale(np.copy(ens_v)[np.where(ens_r<=BIN_R200[k])],9.0)

		# Run Caustic Technique!
		en_caumass,en_caumass_est,en_causurf,en_nfwsurf = self.S.kernel_caustic_masscalc(ens_r,ens_v,HaloData.T[k],BinData.T[k],en_hvd,k)

		# Append Ensemble Data Arrays
		ens_hvd.append(en_hvd)
		ens_caumass.append(en_caumass)
		ens_caumass_est.append(en_caumass_est)

		# Turn into numpy arrays
		ens_r,ens_v,ens_gmags,ens_rmags,ens_imags = np.array(ens_r),np.array(ens_v),np.array(ens_gmags),np.array(ens_rmags),np.array(ens_imags)
		ens_hvd,ens_caumass,ens_caumass_est = np.array(ens_hvd),np.array(ens_caumass),np.array(ens_caumass_est)
		ens_causurf,ens_nfwsurf = np.array(en_causurf),np.array(en_nfwsurf)
		los_r,los_v,los_gmags,los_rmags,los_imags = np.array(los_r),np.array(los_v),np.array(los_gmags),np.array(los_rmags),np.array(los_imags)
		los_hvd,los_caumass,los_caumass_est = np.array(los_hvd),np.array(los_caumass),np.array(los_caumass_est)
		los_causurf,los_nfwsurf = np.array(los_causurf),np.array(los_nfwsurf)
		sample_size,pro_pos = np.array(sample_size),np.array(pro_pos)

		return ens_r,ens_v,ens_gmags,ens_rmags,ens_imags,ens_hvd,ens_caumass,ens_caumass_est,ens_causurf,ens_nfwsurf,los_r,los_v,los_gmags,los_rmags,los_imags,los_hvd,los_caumass,los_caumass_est,los_causurf,los_nfwsurf,self.C.x_range,sample_size,pro_pos	



			







