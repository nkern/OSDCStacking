# Fit Gaussian to Histogram

from scipy.stats import norm
import matplotlib.mlab as mlab


fig = mp.figure(figsize=(12,9))
fig.suptitle("Self Stack Ensemble Velocity Dispersion Histograms, Method 1, Ngal="+str(varib['gal_num'])+", LOS="+str(varib['line_num']))
j = 0
for i in np.arange(0,140,35):
	j += 1
	ax = fig.add_subplot('22'+str(j))	
	x = ENS_V[i]
	n,bins,patches = ax.hist(x,range=(-x.max(),x.max()),color='SlateGrey',alpha=.7,align='left',normed=True)
	p1 = mp.Rectangle((0,0),0,0,fc='SlateGrey',alpha=.7)
	ax.add_patch(p1)
	(mu,sigma) = norm.fit(x)
	y = mlab.normpdf(bins,mu,sigma)
	p2, = ax.plot(bins,y,'r',linewidth=3)
	ax.set_xlim(-x.max(),x.max())
	ax.set_ylim(0,np.max(y)+np.max(y)/2)
	ax.legend([p1,p2],["Pec. VD Hist.","Gaussian Fit"],prop={'size':10})
	ax.set_xlabel('differential velocity (km/s)',fontsize=10)
	ax.set_ylabel('p(v)',fontsize=10)
	ax.set_title( "Ensemble "+str(i),fontsize=11)

mp.show()

