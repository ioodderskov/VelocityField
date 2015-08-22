from __future__ import division
import scipy as sp
import sys
sys.path.insert(0, '/home/io/Dropbox/Projekter/Hubble/VelocityField/code')
import hubble_classes as hc
import hubble_functions as hf
#from gplot import Plot
#plt = Plot('PRL')
import matplotlib.pyplot as plt

path = '/home/io/Dropbox/Projekter/Hubble/VelocityField/cases/Planck512/'
parameters = sp.load(path+'parameters.save')
observers = sp.load(path+'observers.npy')


counts = sp.zeros_like(parameters.bindistances[0:-1])
mean_dists = sp.zeros_like(parameters.bindistances[0:-1])
widths = sp.zeros_like(parameters.bindistances[0:-1])
binnumber = 0
for mind_bin, maxd_bin in zip(parameters.bindistances[0:-1],parameters.bindistances[1:]):
    observed_halos_in_bin = [chosen_halo for chosen_halo in observers[0].chosen_halos\
                            if (mind_bin < chosen_halo.r) & (chosen_halo.r < maxd_bin) ]
    counts[binnumber] = len(observed_halos_in_bin)
    mean_dists[binnumber] = sp.mean((mind_bin,maxd_bin))
    widths[binnumber] = maxd_bin-mind_bin-0.8
    binnumber = binnumber+1                                

plt.figure(figsize=(15,6))
plt.bar(mean_dists,counts/1000,width=widths,align='center',color='grey',edgecolor='black',linewidth=2)
#plt.plot(mean_dists,counts)
ax = plt.gca()
ax.set_xticks(mean_dists)

def tick_function(ds):
    return ["%.0f" % d for d in ds]

ax.set_xticklabels(tick_function(mean_dists))

plt.axis([parameters.bindistances[0]-1,parameters.bindistances[-1]+1,0,3e1])
plt.xlabel('Bin distance [Mpc/h]')
plt.ylabel('Number of SNe [thousands]')
#plt.finalize()
plt.savefig('/home/io/Dropbox/SharedStuff/LSSTVelocityField/bindistances.pdf')
plt.show()