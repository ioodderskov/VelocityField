from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt



zmin = 0.023
zmax = 0.1
MV = -19.504


def plot_SNdata(z,mB):
    log_cz = sp.log10(3e5*z)
    plt.plot(0.2*mB,log_cz,'ro')
    plt.axis([2.6,4.0,3.3,4.7])



def get_z_distribution():

    table_folder = '/home/io/Dropbox/Projekter/Hubble/VelocityField/SNdata/'
    N_tot = 0
    z_tot = sp.array([])

   
    z_all = sp.loadtxt(table_folder+'Galaxies_Hicken2009a_z',usecols=(5,))
    z = sp.array([z_all[i] for i in range(len(z_all)) if z_all[i] < zmax and z_all[i] > zmin])
    N_tot = N_tot+len(z)
    z_tot = sp.concatenate((z_tot,z)) 

    

    
    print "The number of supernovae is", N_tot
    
#    plt.figure()
    
    z_distribution = plt.hist(z_tot,sp.linspace(zmin,zmax,20))
    Wz = z_distribution[0]
    zbins = z_distribution[1]
    plt.show()
    
    return Wz, zbins


