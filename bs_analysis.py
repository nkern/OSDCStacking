
## To check for completeness of runs
import numpy as np
for (j,k) in zip(np.arange(1,50),[2,5,10,15,25,50,100]*7):
	for i in range(2100/k):
		try:
			f = open('bs_m0_run'+str(j)+'/Ensemble_'+str(i)+'_Data.pkl','rb')
		except:
			print 'j =',j,',i =',i


## Make Ensemble Mass 1-1 plot
mp.loglog([BINM200[0],BINM200[-1]],[BINM200[0],BINM200[-1]],'b')
mp.loglog(BINM200,ENS_CAUMASS,'ko',alpha=.7)
mp.xlim(6e13, 1.3e15)
mp.ylim(3.5e13,1.3e15)
mp.xlabel('Median of Bin, Binned on Table Mass',fontsize=15)
mp.ylabel('Ensemble Mass Estimate',fontsize=15)
mp.title('1-1 Correlation, Ngal='+str(gal_num)+', LOS='+str(line_num)+', ENS Richness='+str(line_num*gal_num))
mp.figtext(.65,.15,'Cell #'+str(cell_num)+'\nMethod #'+str(varib['method_num'])+'\nBias = '+str(np.around(ens_mbias,2))+'\nScatter = '+str(np.around(ens_mscat,2)))
mp.show()

## Make Ensemble Mass_Est 1-1 plot
mp.loglog([BIN_M200[0],BIN_M200[-1]],[BIN_M200[0],BIN_M200[-1]],'b')
mp.loglog(BIN_M200,ENS_CAUMASS_EST,'ro',alpha=.7)
mp.xlim(6e13, 1.3e15)
mp.ylim(3.5e13,1.3e15)
mp.xlabel('Median of Bin, Binned on Table Mass',fontsize=15)
mp.ylabel('Ensemble Mass Estimate',fontsize=15)
mp.title('1-1 Correlation Using R200 Estimation, Ngal='+str(gal_num)+', LOS='+str(line_num)+', ENS Richness='+str(line_num*gal_num))
#mp.figtext(.65,.15,'Cell #'+str(cell_num)+'\nMethod #'+str(varib['method_num'])+'\nBias = '+str(np.around(ens_mbias,2))+'\nScatter = '+str(np.around(ens_mscat,2)))
mp.show()


# Make Ensemble HVD 1-1 plot
mp.loglog([BIN_HVD[0],BIN_HVD[-1]],[BIN_HVD[0],BIN_HVD[-1]],'b')
mp.loglog(BIN_HVD,ENS_HVD,'ko')
mp.xlim(BIN_HVD[-1]-50,BIN_HVD[0]+100)
mp.ylim(BIN_HVD[-1]-50,BIN_HVD[0]+100)
mp.xlabel('Median of Bin, Binned on Table HVD',fontsize=15)
mp.ylabel('Ensemble HVD Estimate',fontsize=15)
mp.title('1-1 Correlation, Ngal='+str(gal_num)+', LOS='+str(line_num)+', ENS Richness='+str(line_num*gal_num))
mp.figtext(.7,.2,'Cell #'+str(cell_num)+'\nMethod #'+str(varib['method_num'])+'\nBias = '+str(np.around(ens_vbias,2))+'\nScatter = '+str(np.around(ens_vscat,2)))
mp.show()

# Make LOS Mass 1-1 Plot

mp.loglog([M_crit200[0],M_crit200[-1]],[M_crit200[0],M_crit200[-1]],'b')
mp.loglog(M_crit200[0:halo_num],LOS_CAUMASS.ravel(),'ko',alpha=.7)
mp.xlim(6e13, 1.3e15)
mp.ylim(3.5e13,1.3e15)
mp.xlabel('Median of Bin, Binned on Table Mass',fontsize=15)
mp.ylabel('Ensemble Mass Estimate',fontsize=15)
mp.title('1-1 Correlation, Ngal='+str(gal_num)+', LOS='+str(line_num)+', ENS Richness='+str(line_num*gal_num))
mp.figtext(.65,.15,'Cell #'+str(cell_num)+'\nMethod #'+str(varib['method_num'])+'\nBias = '+str(np.around(los_mbias,2))+'\nScatter = '+str(np.around(los_mscat,2)))
mp.show()


# Make Phase Space w/ Members plot
i = 0
mems = np.array(ENS_MEM[i],bool)
p1, = mp.plot(ENS_R[i][mems],ENS_V[i][mems],'k.')
p2, = mp.plot(ENS_R[i][~mems],ENS_V[i][~mems],'c.',alpha=.4)
p3, = mp.plot(x_range,ENS_CAUSURF[i],'b')
p4 = mp.axvline(BINR200[i],ymax=.2,c='r')
mp.xlabel('cluster-centric radius (Mpc)',fontsize=16)
mp.ylabel('line of sight velocity (km/s)',fontsize=16)
mp.xlim(0,BINR200[i]*2)
mp.legend([p1,p2,p3,p4],["members","non-members","caustic surf","r200"],prop={'size':11},loc=4)
mp.title('Phase Space, Ngal='+str(gal_num)+', Nclus='+str(line_num))

### Bootstrap 1-1 plots ###
#cell_num can be 9,13,23,27
cell_num = 23
ENS_CAUMASS = np.array(map(np.median,d['CAUMASS_cell'+str(cell_num)].T))
BIN_M200 = np.array(map(np.median,d['BINM200_cell'+str(cell_num)].T))
ENS_HVD = np.array(map(np.median,d['HVD_cell'+str(cell_num)].T))
BIN_HVD = np.array(map(np.median,d['BINHVD_cell'+str(cell_num)].T))
ENS_CAUMASS_ERR = np.array(map(np.std,d['CAUMASS_cell'+str(cell_num)].T))
ENS_HVD_ERR = np.array(map(np.std,d['HVD_cell'+str(cell_num)].T))
BIN_M200_ERR = np.array(map(np.median,d['BINM200_STD_cell'+str(cell_num)].T))
BIN_HVD_ERR = np.array(map(np.median,d['BINHVD_STD_cell'+str(cell_num)].T))

# Mass 1-1
p1 = mp.errorbar(BIN_M200,ENS_CAUMASS,yerr=ENS_CAUMASS_ERR,xerr=BIN_M200_ERR,fmt='.',color='k',alpha=.4)
p2, = mp.plot([BIN_M200[0],BIN_M200[-1]],[BIN_M200[0],BIN_M200[-1]],'b',linewidth=2)
mp.xscale('log')
mp.yscale('log')
mp.xlim(6e13, 1.3e15)
mp.ylim(3.5e13,1.3e15)

# HVD 1-1
p1 = mp.errorbar(BIN_HVD,ENS_HVD,yerr=ENS_HVD_ERR,xerr=BIN_HVD_ERR,fmt='.',color='k',alpha=.4)
p2, = mp.plot([BIN_HVD[0],BIN_HVD[-1]],[BIN_HVD[0],BIN_HVD[-1]],'b',linewidth=2)
mp.xlim(BIN_HVD[-1]-50,BIN_HVD[0]+100)
mp.ylim(BIN_HVD[-1]-50,BIN_HVD[0]+100)
mp.xscale('log')
mp.yscale('log')














