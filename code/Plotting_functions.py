from __future__ import division
import numpy as sp
import sys
import matplotlib.pyplot as plt
sys.path.insert(0,'/home/io/Dropbox/Projekter/Hubble/VelocityField/SNdata')

plt.rc('font',family = 'serif')








def plot_redshiftdistribution(histogram, alpha, color):
    
    Wz = histogram[0]
    zbins = histogram[1]
    width = zbins[1]-zbins[0]
    center = (zbins[:-1]+zbins[1:])/2
    plt.bar(center, Wz, align = 'center', width=width, alpha=alpha, color=color)

    plt.xlabel('Redshift',fontsize=16)
    plt.ylabel('Fraction of SN Ia',fontsize=16)
    plt.axis([0,0.1,0,0.2])



def plot_hubblediagram(radial_distances,radial_velocities):

    
    radial_distances = sp.array(radial_distances)
    radial_velocities = sp.array(radial_velocities)

    c = 3e5
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(radial_distances, radial_velocities/c,'bo') 


    plt.axis([0,300,0,0.1])
    plt.xlabel('$r$ [Mpc/h]',fontsize=16)
    plt.ylabel('Redshift',fontsize=16)


plot_redshiftdistribution_on=0
plot_hubblediagram_on=1    
    
if plot_redshiftdistribution_on:
    
    plot_redshiftdistribution(sp.load('histogram_table.npy'),0.7, 'green')
    plot_redshiftdistribution(sp.load('histogram_mock.npy'),0.3, 'blue')

#    plt.savefig('/home/io/Dropbox/SharedStuff/hubble2013/redshiftdistribution.pdf')

    

if plot_hubblediagram_on:
    
    radial_distances = sp.load('radial_distances.npy')
    radial_velocities = sp.load('radial_velocities.npy')
    
    plot_hubblediagram(radial_distances,radial_velocities)
    
#    plt.savefig('/home/io/Dropbox/SharedStuff/hubble2013/hubblediagram.pdf')



