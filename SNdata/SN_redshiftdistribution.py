from __future__ import division
import numpy as sp
# Chosen plot options
#from gplot import Plot 
#plt = Plot('latex_full_hubble')
import matplotlib.pyplot as mplt
mplt.rc('font',family = 'serif')


zmin = 0.023
zmax = 0.1
MV = -19.504
Nbins = 20


def plot_SNdata(z,mB):
    log_cz = sp.log10(3e5*z)
    mplt.plot(0.2*mB,log_cz,'ro')
    mplt.axis([2.6,4.0,3.3,4.7])



def get_z_distribution():

    table_folder = '../SNdata/'
    N_tot = 0
    z_tot = sp.array([])

   
    z_all = sp.loadtxt(table_folder+'Hicken2009a_galaxies',usecols=(5,))
    z = sp.array([z_all[i] for i in range(len(z_all)) if z_all[i] < zmax and z_all[i] > zmin])
    N_tot = N_tot+len(z)
    z_tot = sp.concatenate((z_tot,z)) 

    z_all = sp.loadtxt(table_folder+'Jha2007_galaxies',usecols=(2,))
    z = sp.array([z_all[i] for i in range(len(z_all)) if z_all[i] < zmax and z_all[i] > zmin])
    N_tot = N_tot+len(z)
    z_tot = sp.concatenate((z_tot,z)) 
    

    
    print "The number of supernovae is", N_tot

    weights=sp.ones_like(z_tot)/len(z_tot)
    z_distribution = mplt.hist(z_tot,sp.linspace(zmin,zmax,Nbins), alpha = 0.7, color='green', weights=weights)
    Wz = z_distribution[0]
    zbins = z_distribution[1]


    
    return Wz, zbins, N_tot


def make_histograms(radial_velocities):

    Wz_table, zbins_table, N_table = get_z_distribution()    

    c = 3e5 # km/s
    z_mock = sp.array(radial_velocities)/c
    weights = sp.ones_like(z_mock)/len(z_mock)
    z_dist = mplt.hist(z_mock,zbins_table, alpha = 0.3, color='blue', weights=weights)
    Wz_mock = z_dist[0]
    z_mock = z_dist[1]
    mplt.xlabel('Redshift',fontsize=16)
    mplt.ylabel('Fraction of SN Ia',fontsize=16)
    mplt.axis([0,0.1,0,0.2])

    return Wz_mock, z_mock, len(z_mock)

get_z_distribution()
mplt.xlabel('Redshift',fontsize=16)
mplt.ylabel('Fraction of SN Ia',fontsize=16)
mplt.axis([0,0.1,0,0.2])
#mplt.figure(figsize=(1,2))

#import matplotlib.pyplot as mplt
#fig = mplt.figure()
#ax = fig.add_subplot(1,1,1)
#Wz, zbins = get_z_distribution()
#mplt.xlabel('Redshift')
#mplt.ylabel('Fraction of SN Ia')
#mplt.axis([0.01,0.09,0,0.2])

