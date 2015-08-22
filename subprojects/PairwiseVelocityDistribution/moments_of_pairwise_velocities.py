from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import sys


def plot_pairwise_velocities_mass(cases,color):

    #central_halo_masses = ['3.5e11']
    #central_halo_masses = ['3.50e+11','9.98e+11','5.00e+12','2.50e+13','5.20e+14']
    central_halo_masses = ['2.00e+11','1.08e+12','6.50e+12','8.00e+13','5.75e+14']
    double_central_halo_masses = [sp.double(central_halo_mass) for central_halo_mass in central_halo_masses]
    #path = '../cases/'+case+'/ROCKSTAR_'
    #path = '../cases/'+case+'/'
    path = '../cases/'+case
    
    Rs = [1,5]
    dRs = [1,0.2]
    
    
    round = 0
    subplots = [221,222]
    for R,dR,subplot in zip(Rs,dRs,subplots):
        v12_of_masses = []
        sigma_pp_of_masses = []
        Rmin, Rmax = R-dR/2, R+dR/2
        for central_halo_mass in central_halo_masses:    
            pairwise_velocities_file = path+'pairwise_velocities_'+central_halo_mass+'.npy'
            radial_distances_file = path+'radial_distances_'+central_halo_mass+'.npy'
            
            pairwise_velocities = sp.load(pairwise_velocities_file)
            radial_distances = sp.load(radial_distances_file)
    
            if round == 0:
                if not 'all_pairwise_velocities' in locals():
                    all_pairwise_velocities = pairwise_velocities
                else:
                    all_pairwise_velocities = sp.hstack((all_pairwise_velocities,pairwise_velocities))
                if not 'all_radial_distances' in locals():
                    all_radial_distances = radial_distances
                else:
                    all_radial_distances = sp.hstack((all_radial_distances,radial_distances))
    
            pairwise_velocities_R = sp.array([pairwise_velocity\
                    for pairwise_velocity,radial_distance in zip(pairwise_velocities,radial_distances)\
                    if (Rmin < radial_distance) & (radial_distance < Rmax)])
    
            print "len(pairwise_velocities_R) = ", len(pairwise_velocities_R)
            v12 = -sp.mean(pairwise_velocities_R)
            sigma_pp = sp.sqrt(sp.mean(pairwise_velocities_R**2))
    
            v12_of_masses.append(v12)
            sigma_pp_of_masses.append(sigma_pp)
    
        plt.subplot(subplot) 
        plt.plot(double_central_halo_masses,sigma_pp_of_masses,'.-',color=color,label=case)
        if subplot == 221:
            plt.ylabel('$\sigma_{||}$ [km/s]')
            plt.title('R=1Mpc/h')
            plt.legend(loc=2,prop={'size':8})
        if subplot == 222:
            plt.title('R=5Mpc/h')
        plt.axis([1e11,1e15,0,600])
        plt.xscale('log')
    
        plt.subplot(subplot+2) 
        plt.plot(double_central_halo_masses,v12_of_masses,'.-',color=color,label=case)
        if (subplot == 221) | (subplot == 222):
            plt.xlabel('$M_{200}$ [$M_{sun}$/h]')
        if subplot == 221:
            plt.ylabel('$-v_{12}$ [km/s]')
        plt.axis([1e11,1e15,0,600])
        plt.xscale('log')
    
        round = round+1

        return all_radial_distances, all_pairwise_velocities

    


def plot_pairwise_velocities_r(case,color,all_radial_distances,all_radial_velocities):
    dr = 0.3 # Mpc/h
    rmin, rmax = sp.amin(all_radial_distances), sp.amax(all_radial_distances) 
    rrange = rmax-rmin 
    N = int(sp.ceil(rrange/dr))
    rs = sp.linspace(rmin,rmax,N)
    v12_of_r = [[] for index in range(N)]
    
    for r,v12 in zip(all_radial_distances,all_pairwise_velocities):
    
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
    
    
    plt.plot(rs,sigma_12s,color=color,label='$\sigma_{12}$')
    plt.plot(rs,v12_means,color=color,label='$|v_{12}|$')
    plt.xlabel('r [Mpc/h]')
    plt.ylabel('[km/s]')
    plt.xscale('log')
    plt.axis([0.5,100,0,600])
    

#if len(sys.argv) != 2:
#    print "Wrong number of arguments"
#
#case = sys.argv[1]

print "Plotting -v12 and sigma_pp as a function of pair masses"
plt.figure()
#cases = ['CoDECS_LCDM/','CoDECS_EXP001/','CoDECS_EXP003/']
cases = ['LyA_WMAP/ROCKSTAR_','LyA_Planck/ROCKSTAR_']
#cases = ['CoDECS_LCDM/','CoDECS_LCDM/ROCKSTAR_']
#colors = ['black','b','g']
colors = ['b','g']
for case,color,number in zip(cases,colors,range(len(cases))):
    if number == 0:
        all_radial_distances_0, all_pairwise_velocities_0 = plot_pairwise_velocities_mass(case,color)
    if number == 1:
        all_radial_distances_1, all_pairwise_velocities_1 = plot_pairwise_velocities_mass(case,color)
    if number == 2:
        all_radial_distances_2, all_pairwise_velocities_2 = plot_pairwise_velocities_mass(case,color)

plt.savefig('pairwise_velocities_mass.pdf')



print "Plotting -v12 and sigma_pp as a function of r"
plt.figure(figsize=(6,5))
for case,color,number in zip(cases,colors,range(len(cases))):
    if number == 0:
        all_radial_distances, all_pairwise_velocities = all_radial_distances_0, all_pairwise_velocities_0
    if number == 1:
        all_radial_distances, all_pairwise_velocities = all_radial_distances_1, all_pairwise_velocities_1
    if number == 2:
        all_radial_distances, all_pairwise_velocities = all_radial_distances_2, all_pairwise_velocities_2
    plot_pairwise_velocities_r(case,color,all_radial_distances,all_pairwise_velocities)

plt.savefig('pairwise_velocities_r.pdf')
