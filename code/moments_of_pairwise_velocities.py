from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import sys

calculate_moments = 1

if len(sys.argv) != 2:
    print "Wrong number of arguments"

case = sys.argv[1]
path = '../cases/'+case+'/'

if calculate_moments:

    
    pairwise_velocities_file = path+'pairwise_velocities.npy'
    radial_distances_file = path+'radial_distances.npy'
    pair_masses_file = path+'pair_masses.npy'
    
    pairwise_velocities = sp.load(pairwise_velocities_file)
    radial_distances = sp.load(radial_distances_file)
    pair_masses = sp.load(pair_masses_file)
    print "pair_masses = ", pair_masses
    print "min(pair_masses) = ", sp.amin(pair_masses)
    print "max(pair_masses) = ", sp.amax(pair_masses)
#    print "radial_distances = ", radial_distances
#    print "pairwise_velocities = ", pairwise_velocities

#    print "sqrt(mean(pairwise_velocities)**2) = ", sp.sqrt(sp.mean(pairwise_velocities**2))
#    print "mean(pairwise_velocities) = ", sp.mean(pairwise_velocities)
    
    
    
    dr = 1 # Mpc/h
    rmin, rmax = sp.amin(radial_distances), sp.amax(radial_distances) 
    rrange = rmax-rmin 
    N = int(sp.ceil(rrange/dr))
    rs = sp.linspace(rmin,rmax,N)
    v12_of_r = [[] for index in range(N)]
    
    for r,v12,pair_mass in zip(radial_distances,pairwise_velocities,pair_masses):

        #if any(pair_mass <= 1e11) | any(pair_mass >= 1e15):
#        if any(pair_mass <= 1e11):
#            continue

#        if pair_mass[0] <= 1e11:
#            print "Something is wrong"
#        if pair_mass[1] <= 1e11:
#            print "Something is wrong"
#        if pair_mass[0] >= 1e15:
#            print "Something is wrong"
#        if pair_mass[1] >= 1e15:
#            print "Something is wrong"
        index = int(sp.floor((r-rmin)/dr))
        v12_of_r[index].append(v12)
            
    
    sigma_12s = sp.zeros(N)
    v12_means = sp.zeros(N)
    for index in range(len(sigma_12s)):
        v12_of_r_index = sp.array(v12_of_r[index])
        print "number of counts in the", index,"th bin:", len(v12_of_r_index)
        sigma_12 = sp.sqrt(sp.mean(v12_of_r_index**2))
        v12_mean = -sp.mean(v12_of_r_index)
        sigma_12s[index] = sigma_12
        v12_means[index] = v12_mean

    sp.save(path+'rs.npy',rs)
    sp.save(path+'sigma_12s.npy',sigma_12s)
    sp.save(path+'v12_means.npy',v12_means)

else:

	rs = sp.load(path+'rs.npy')
	sigma_12s = sp.load(path+'sigma_12s.npy')
	v12_means = sp.load(path+'v12_means.npy')

plt.figure(figsize=(6,5))
ax1=plt.subplot(111)

ax1.plot(rs,sigma_12s,'b',label='$\sigma_{12}$')
ax1.plot(rs,v12_means,'g',label='$|v_{12}|$')
plt.xlabel('r [Mpc/h]')
plt.ylabel('[km/s]')
plt.xscale('log')
plt.xlim([0.1,10])
#plt.axis([0.1,100,0,900])
plt.legend()
plt.savefig(path+'v12_of_r.pdf')



