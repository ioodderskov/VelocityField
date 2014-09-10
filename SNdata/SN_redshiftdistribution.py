from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt



zmin = 0.023
zmax = 0.1
MV = -19.504

def load_from_table(table,col1,col2):

    data = sp.loadtxt(table,usecols=(col1,col2))
    z_all = data[:,0]
    m_all = data[:,1]
    z = sp.array([z_all[i] for i in range(len(z_all)) if z_all[i] < zmax and z_all[i] > zmin])
    m = sp.array([m_all[i] for i in range(len(z_all)) if z_all[i] < zmax and z_all[i] > zmin])

    print "The number of supernovae in",table,"is",len(z)
    return z,m,len(z)

def plot_SNdata(z,mB):
    log_cz = sp.log10(3e5*z)
    plt.plot(0.2*mB,log_cz,'ro')



def get_z_distribution():

    table_folder = '/home/io/Dropbox/Projekter/Hubble/VelocityField/SNdata/'
    N_tot = 0
    z_tot = sp.array([])

#    z_tab0, mV_tab0, N_tab0 = load_from_table(table_folder+'Hicken2009a_table8.txt',1,4)
#    plot_SNdata(z_tab0, mV_tab0)
#    N_tot = N_tot+N_tab0
#    z_tot = sp.concatenate((z_tot,z_tab0)) 
   
#    z_tab1, mB_tab1, N_tab1 = load_from_table(table_folder+'Hickens_SNe_table1',1,3)
#    plot_SNdata(z_tab1, mB_tab1)
#    N_tot = N_tot+N_tab1
#    z_tot = sp.concatenate((z_tot,z_tab1))
#    
#    z_tab2, mB_tab2, N_tab2 = load_from_table(table_folder+'Hickens_SNe_table2',1,3)
#    plot_SNdata(z_tab2, mB_tab2)
#    N_tot = N_tot+N_tab2
#    z_tot = sp.concatenate((z_tot,z_tab2))
#    #
#    z_tab3, muA_tab3, N_tab3 = load_from_table(table_folder+'Hickens_SNe_table3',1,5)
#    mA_tab3 = muA_tab3+MV
#    plot_SNdata(z_tab3, mA_tab3)
#    N_tot = N_tot+N_tab3
#    z_tot = sp.concatenate((z_tot,z_tab3))
##    #
#    z_tab4, muA_tab4, N_tab4 = load_from_table(table_folder+'Hickens_SNe_table4',1,5)
#    mA_tab4 = muA_tab4+MV
#    plot_SNdata(z_tab4, mA_tab4)
#    N_tot = N_tot+N_tab4
#    z_tot = sp.concatenate((z_tot,z_tab4))
   
    z_all = sp.loadtxt(table_folder+'Galaxies_Hicken2009a_z',usecols=(5,))
    z = sp.array([z_all[i] for i in range(len(z_all)) if z_all[i] < zmax and z_all[i] > zmin])
    N_tot = N_tot+len(z)
    z_tot = sp.concatenate((z_tot,z)) 

    
    plt.axis([2.6,4.0,3.3,4.7])
    
    print "The number of supernovae is", N_tot
    
    plt.figure()
    
    z_distribution = plt.hist(z_tot,sp.linspace(zmin,zmax,20))
    Wz = z_distribution[0]
    zbins = z_distribution[1]
    plt.show()
    
    return Wz, zbins


#
#f = open('SNdistribution.txt','w')
#f.write("%s\t" % bl)
