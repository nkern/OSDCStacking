import cPickle as pkl
import matplotlib.pyplot as mp
import numpy.ma as ma
import numpy as np

file = open('table_analysis_binstack_run_table2.pkl','rb')
input = pkl.Unpickler(file)
data = input.load()

(ENS_MBIAS,ENS_MSCAT,ENS_VBIAS,ENS_VSCAT,LOS_MBIAS,LOS_MSCAT,LOS_VBIAS,LOS_VSCAT,RUN_NUM,GAL_NUM,LINE_NUM,RICH_NUM) = data

## add LOS = 1 table data
los1_mscat = [1.17,0.78,0.63,0.48,0.38,0.31,0.30]

## Gifford et al. 2013 Guo 2010 Individual Cluster Stats
IND_RICH 		= np.array([10,15,25,35,50,75,100,150],int)
IND_MSCAT		= np.array([.94,.73,.53,.46,.41,.36,.33,.31])
IND_MBIAS 		= np.array([-.64,-.39,-.16,-.10,-.03,.00,.01,.01])
IND_MBIAS_ERR 		= np.array([.09,.08,.05,.05,.04,.03,.03,.03])
##

raise NameError

mp.errorbar(IND_RICH,IND_MBIAS,yerr=IND_MBIAS_ERR,lw=2,c='DarkBlue',alpha=.7)

## ENSEMBLE MASS SCATTER
fig,axes = mp.subplots(1,1)
p1, = axes.plot(RICH_NUM[0],ENS_MSCAT[0],c='DarkBlue',marker='o',alpha=.8)
axes.plot(RICH_NUM[0],ENS_MSCAT[0],'DarkBlue',alpha=.8)
p2, = axes.plot(RICH_NUM[1],ENS_MSCAT[1],c='DarkGreen',marker='o',alpha=.8)
axes.plot(RICH_NUM[1],ENS_MSCAT[1],'DarkGreen',alpha=.8)
p3, = axes.plot(RICH_NUM[2],ENS_MSCAT[2],c='FireBrick',marker='o',alpha=.8)
axes.plot(RICH_NUM[2],ENS_MSCAT[2],'FireBrick',alpha=.8)
p4, = axes.plot(RICH_NUM[3],ENS_MSCAT[3],c='Indigo',marker='o',alpha=.8)
axes.plot(RICH_NUM[3],ENS_MSCAT[3],'Indigo',alpha=.8)
p5, = axes.plot(RICH_NUM[4],ENS_MSCAT[4],c='DarkSlateGrey',marker='o',alpha=.8)
axes.plot(RICH_NUM[4],ENS_MSCAT[4],'DarkSlateGrey',alpha=.8)
p6, = axes.plot(IND_RICH,IND_MSCAT,'ko',ms=8,alpha=.8)
axes.plot(IND_RICH,IND_MSCAT,'k',lw=3,alpha=.4)
axes.set_xlabel('Ensemble Richness, N, within R200',fontsize=14)
axes.set_ylabel('Mass Estimate Scatter',fontsize=14)
axes.set_title('Bin Stacked Ensemble Mass Scatter (2100 Halos) 2/20/14')
axes.legend([p1,p2,p3,p4,p5,p6],["N=5, 2<LOS<100","N=10, 2<LOS<100","N=15, 2<LOS<100","N=25, 2<LOS<100","N=50, 2<LOS<100","Gifford 2013"],loc=1)
axes.set_xlim(-10,1000)
axes.set_ylim(0,1.0)
axes.locator_params(axis='both',nbins=12)
##



## ENSEMBLE MASS BIAS
fig,axes = mp.subplots(1,1)
p1, = axes.plot(RICH_NUM[0],ENS_MBIAS[0],c='DarkBlue',marker='o',alpha=.7)
axes.plot(RICH_NUM[0],ENS_MBIAS[0],'DarkBlue',alpha=.7)
p2, = mp.plot(RICH_NUM[1],ENS_MBIAS[1],c='DarkGreen',marker='o',alpha=.7)
axes.plot(RICH_NUM[1],ENS_MBIAS[1],'DarkGreen',alpha=.7)
p3, = axes.plot(RICH_NUM[2],ENS_MBIAS[2],c='FireBrick',marker='o',alpha=.7)
axes.plot(RICH_NUM[2],ENS_MBIAS[2],'FireBrick',alpha=.7)
p4, = mp.plot(RICH_NUM[3],ENS_MBIAS[3],c='Indigo',marker='o',alpha=.7)
axes.plot(RICH_NUM[3],ENS_MBIAS[3],'Indigo',alpha=.7)
p5, = axes.plot(RICH_NUM[4],ENS_MBIAS[4],c='DarkSlateGrey',marker='o',alpha=.7)
axes.plot(RICH_NUM[4],ENS_MBIAS[4],'DarkSlateGrey',alpha=.7)
p6 = mp.errorbar(IND_RICH,IND_MBIAS,yerr=IND_MBIAS_ERR,lw=2,c='k',alpha=.6)
axes.set_xlabel('Ensemble Richness within R200',fontsize=14)
axes.set_ylabel('Mass Estimate Bias',fontsize=14)
axes.set_title('Bin Stacked Ensemble Mass Bias (2100 Halos) 2/20/14')
axes.legend([p1,p2,p3,p4,p5,p6],["N=5, 2<LOS<100","N=10, 2<LOS<100","N=15, 2<LOS<100","N=25, 2<LOS<100","N=50, 2<LOS<100","Gifford 2013"],loc=7)
axes.set_ylim(-.9,.1)
axes.set_xlim(0,300)
axes.locator_params(axis='both',nbins=12)
##


## LOS MASS SCATTER
fig,axes = mp.subplots(1,1)
p1, = axes.plot(IND_RICH,IND_MSCAT,'ko',ms=8,alpha=.8)
axes.plot(IND_RICH,IND_MSCAT,'k',lw=3,alpha=.4)
p2, = axes.plot([5,10,15,25,50,100,150],los1_mscat,c='r',marker='o')
axes.plot([5,10,15,25,50,100,150],los1_mscat,c='r')
axes.set_xlabel('Line of Sight Sampling Richness',fontsize=14)
axes.set_ylabel('Mass Estimate Scatter',fontsize=14)
axes.set_title('Averaged Systems Mass Scatter (2124 Halos) 12/23/13')
axes.legend([p1,p2],["Gifford 2013","LOS=1"],loc=1)
axes.set_xlim(0,160)
axes.set_ylim(-.05,1.3)
axes.locator_params(axis='both',nbins=12)
##


## LOS MASS BIAS
fig,axes = mp.subplots(1,1)
p8 = axes.errorbar(IND_RICH,IND_MBIAS,yerr=IND_MBIAS_ERR,lw=2,c='k',alpha=.5)
p1, = axes.plot(GAL_NUM[0],LOS_MBIAS[0],c='DarkBlue',marker='o',alpha=.7)
axes.plot(GAL_NUM[0],LOS_MBIAS[0],'DarkBlue',alpha=.3)
p2, = axes.plot(GAL_NUM[1],LOS_MBIAS[1],c='DarkGreen',marker='o',alpha=.7)
axes.plot(GAL_NUM[1],LOS_MBIAS[1],'DarkGreen',alpha=.3)
p3, = axes.plot(GAL_NUM[2],LOS_MBIAS[2],c='FireBrick',marker='o',alpha=.7)
axes.plot(GAL_NUM[2],LOS_MBIAS[2],'FireBrick',alpha=.3)
p4, = axes.plot(GAL_NUM[3],LOS_MBIAS[3],c='Indigo',marker='o',alpha=.7)
axes.plot(GAL_NUM[3],LOS_MBIAS[3],'Indigo',alpha=.3)
p5, = axes.plot(GAL_NUM[4],LOS_MBIAS[4],c='DarkSlateGrey',marker='o',alpha=.7)
axes.plot(GAL_NUM[4],LOS_MBIAS[4],'DarkSlateGrey',alpha=.3)
p6, = axes.plot(GAL_NUM[5],LOS_MBIAS[5],c='Magenta',marker='o',alpha=.7)
axes.plot(GAL_NUM[5],LOS_MBIAS[5],c='Magenta',alpha=.3)
p7, = axes.plot(GAL_NUM[6],LOS_MBIAS[6],c='SpringGreen',marker='o',alpha=.7)
axes.plot(GAL_NUM[6],LOS_MBIAS[6],c='SpringGreen',alpha=.3)
axes.set_xlabel('Line of Sight Sampling Richness',fontsize=14)
axes.set_ylabel('Mass Estimate Bias',fontsize=14)
axes.set_title('Averaged Systems Mass Bias (2124 Halos) 12/23/13')
axes.legend([p1,p2,p3,p4,p5,p6,p7,p8],["N=5, 2<LOS<100","N=10, 2<LOS<100","N=15, 2<LOS<100","N=25, 2<LOS<100","N=50, 2<LOS<100","N=100, 2<LOS<100","N=150, 2<LOS<100","Gifford 2013"],loc=7)
axes.set_xlim(0,160)
axes.set_ylim(-1.,.1)
axes.locator_params(axis='both',nbins=12)
##

### ENS_MBIAS Density Map
mp.figure()
mp.imshow(ENS_MBIAS,cmap='jet')
cbar = mp.colorbar()
cbar.set_label('Mass Bias',fontsize=12)
mp.grid()
mp.xticks([0,1,2,3,4,5,6],[2,5,10,15,25,50,100])
mp.yticks([0,1,2,3,4,5,6],[5,10,15,25,50,100,150])
mp.xlabel('Clusters per Bin (LOS)',fontsize=15)
mp.ylabel('Galaxies per Cluster (Ngal)',fontsize=15)
mp.title('Run Table Ensemble Mass Bias Density Plot')

### ENS_MSCAT Density Map
mp.figure()
mp.imshow(ENS_MSCAT,cmap='jet')
cbar = mp.colorbar()
cbar.set_label('Mass Scatter',fontsize=12)
mp.grid()
mp.xticks([0,1,2,3,4,5,6],[2,5,10,15,25,50,100])
mp.yticks([0,1,2,3,4,5,6],[5,10,15,25,50,100,150])
mp.xlabel('Clusters per Bin (LOS)',fontsize=15)
mp.ylabel('Galaxies per Cluster (Ngal)',fontsize=15)
mp.title('Run Table Ensemble Mass Scatter Density Plot')




mp.show()
















