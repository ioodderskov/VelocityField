from __future__ import division
import numpy as sp
import sys
sys.path.insert(0,'/home/io/Dropbox/Projekter/Hubble/VelocityField/SNdata')








def make_hubbleplot(radial_distances,radial_velocities,mind,maxd,boxsize):

    
    radial_distances = sp.array(radial_distances)
    radial_velocities = sp.array(radial_velocities)

    c = 3e5
    #import matplotlib.pyplot as plt
    #plt.rc('font',family = 'serif')
    #
    #fig = plt.figure()
    #ax = fig.add_subplot(1,1,1)
    #ax.plot(radial_distances, radial_velocities/c,'bo') 


    #plt.axis([0,300,0,0.1])
    #plt.xlabel('$r[Mpc/h]$',fontsize=16)
    #plt.ylabel('Redshift',fontsize=16)
    
