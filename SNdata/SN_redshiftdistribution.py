from __future__ import division
import numpy as sp



zmin = 0.01
zmax = 0.1
MV = -19.504
Nbins = 20
bins = sp.linspace(zmin,zmax,Nbins)


def plot_SNdata(z,mB):
    log_cz = sp.log10(3e5*z)
    #mplt.plot(0.2*mB,log_cz,'ro')
    #mplt.axis([2.6,4.0,3.3,4.7])



def get_table_distribution():

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

    histogram_table = sp.histogram(z_tot,bins, weights=weights)
    sp.save('histogram_table.npy',histogram_table)

    Wz = histogram_table[0]
    zbins = histogram_table[1]
    
    return Wz, zbins, N_tot



def get_mock_distribution(radial_velocities):

    c = 3e5 # km/s
    z_mock = radial_velocities/c
    N_tot = len(z_mock)
    weights = sp.ones_like(z_mock)/N_tot

    histogram_mock = sp.histogram(z_mock,bins, weights=weights)
    sp.save('histogram_mock.npy',histogram_mock)

    N_tot = len(z_mock)
    Wz = histogram_mock[0]
    zbins = histogram_mock[1]
    
    return Wz, zbins, N_tot


