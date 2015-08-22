from __future__ import division
import scipy as sp
import healpy as hp
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/home/io/Dropbox/Projekter/Hubble/VelocityField/code')
import hubble_classes as hc
import cPickle

path = '/home/io/Dropbox/Projekter/Hubble/VelocityField/cases/sim16/'
skyfraction = ''
parameters = sp.load(path+'parameters'+skyfraction+'.save')
#parameters = cPickle.load(path+'parameters.save')
observers = sp.load(path+'observers'+skyfraction+'.npy')



bin_number = 0
color = 'b'
plt.figure()
for mind_bin, maxd_bin in zip(parameters.bindistances[0:-1],parameters.bindistances[1:]):
    d = sp.mean((mind_bin,maxd_bin))
    ls = observers[0].ls[bin_number]
    cls = observers[0].cls[bin_number]
    bin_number = bin_number+1

    cls_normed = sp.sqrt(ls*(ls+1) * cls)
#    alpha = 1/(1+bin_number)
    alpha=1
    lw=10/bin_number
    plt.plot(ls, cls_normed,'b',label="d = %.2f" % d,alpha=alpha,lw=lw,color=color)

    
plt.xlabel('$l$'); 
plt.ylabel('$\sqrt{l(l+1)C_l}$ [km/s]'); 
plt.axis([0,20,0,800])
plt.grid()
if skyfraction == '1.0':
    plt.title("Full sky")
else:
    plt.title("Cut sky")
#plt.xscale('log')
#plt.yscale('log')

plt.savefig('/home/io/Dropbox/SharedStuff/LSSTVelocityField/observer_all_bins_'+skyfraction+'.pdf')

