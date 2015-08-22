from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import scipy.special
import mpmath
import pdb

def spherical_bessel(l,x):
    factor = sp.sqrt(sp.pi/(2*x))
    besselj = mpmath.besselj(l+0.5,x)
    return factor*besselj


def velocitypowerspectrum(MPS_file,R,h,omega_m,d_shell):
    #Stuff to read in from the parameterfile
    #camb_transfer_function = '../../cases/Planck512/TF_Planck_BAO_matterpower_z0.dat'
#    camb_matterpowerspectrum = 'test/test_z0_matterpower.dat'
#
#    with open('params.ini', 'r') as file:
#        # read a list of lines into data
#        data = file.readlines()

#    pdb.set_trace()
#    h = sp.double(data[39].split()[2])/100
#    omega_b = sp.double(data[35].split()[2])/h**2
#    omega_cdm = sp.double(data[36].split()[2])/h**2
#    omega_m = omega_b+omega_cdm
    lmax = 100
#    x = 150/h #Mpc

    
    
    end = 596
    data = sp.loadtxt(MPS_file)
    ks, Pmk = data[:end,0], data[:end,1] # De sidste 15 vaerdier af ks giver nan i besselfunc.
    
    #norm = TF[0]
    #TF = TF/norm
    
#    plt.plot(ks,Pmk)
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.xlabel('k [Mpc^-1]')
#    plt.title('Matterpowerspectrum')
    #The powerspectrum looks fine (compared to Binney and Tremaine p. 730)
    #It looks as if the units of CAMB are /Mpc, not h/Mpc (by comparison)
    
    #Pmk = ks*TF**2
    
    
    
    
    Wk2 = 1/(1+ks*R) 
#    Wk2 = 1 # for testing
    
    H = 100 # km/s/Mpc
    f = omega_m**0.6
    
    Pvk = H**2*f**2*ks**-2*Pmk*Wk2
#    pdb.set_trace()
    
    B = sp.zeros(lmax+1)
    C = sp.zeros(lmax)
    spherical_bessel_kx_l = sp.zeros((len(ks),lmax+1))
#    special_spherical_bessel_kx_l = sp.zeros((len(ks),lmax+1))
    
    for row, k in enumerate(ks):
        spherical_bessel_kx_l[row] = sp.special.sph_jn(lmax,k*d_shell)[0]
    #    factor = sp.sqrt(sp.pi/(2*k*x))
    #    for column, l in enumerate(range(lmax)):
    #        spherical_bessel_kx_l[row][column] = spherical_bessel(l,k*x)
    
        
    for l in range(lmax+1):
        integrand = ks**2/(2*sp.pi)**3*Pvk*spherical_bessel_kx_l[:,l]**2
#        pdb.set_trace()
        B[l] = 4*sp.pi*sp.trapz(integrand,ks)
        
    for l in range(1,lmax):
        C[l] = 4*sp.pi/(2*l+1)*(l*B[l-1]+(l+1)*B[l+1])
    
     
    #plt.plot(ks,Pmk)
#    plt.xscale('log')
#    plt.yscale('log')
    #plt.xlabel('$k [Mpc^{-1}]$')
    #plt.ylabel('P(k) $[Mpc^{3}]$')
    #plt.axis([0.001, 10, 1e-4, 1e5])
    
#    plt.figure()
    ls = sp.array(range(1,lmax+1))
    scaled_Cl = sp.sqrt(ls*(ls+1)*C)
    plt.plot(ls[1:],scaled_Cl[1:],label='R = %.2f Mpc/h' % R)
#    plt.xscale('log')
#    plt.yscale('log')
    plt.title('Velocitypowerspectrum')
    plt.xlabel('$l$')
    plt.ylabel('$\sqrt{l(l+1)C_l}$ [km/s]')
    plt.xlim([0,20])
#    plt.ylim([0,800])
#    pdb.set_trace()
    plt.legend(frameon=0)

#h = 0.7
#
#camb_matterpowerspectrum = 'test/test_z0_matterpower.dat'
#MPS_file = camb_matterpowerspectrum
#
#with open('params.ini', 'r') as file:
#    # read a list of lines into data
#    data = file.readlines()
#
#h = sp.double(data[39].split()[2])/100
#omega_b = sp.double(data[35].split()[2])/h**2
#omega_cdm = sp.double(data[36].split()[2])/h**2
#omega_m = omega_b+omega_cdm
#d_shell = 150/h #Mpc
#
#
#velocitypowerspectrum(MPS_file,0,h,omega_m,d_shell)
##velocitypowerspectrum(0.03) #Mpc. Roughly the diameter of the Milky Way
#velocitypowerspectrum(MPS_file,4/h,h,omega_m,d_shell) #Mpc. 8 neighbours (from healpy.get_all_neighbours) and eq. 8 in paper 2
#velocitypowerspectrum(MPS_file,150/h*0.13,h,omega_m,d_shell) # Roughly the diameter of the largest hole
