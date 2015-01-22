from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import scipy.special
import math

def spherical_bessel(l,x):
    factor = sp.sqrt(sp.pi/(2*x))
    besselj = mpmath.besselj(l+0.5,x)
    return factor*besselj


def velocitypowerspectrum(R):
    #Stuff to read in from the parameterfile
    #camb_transfer_function = '../../cases/Planck512/TF_Planck_BAO_matterpower_z0.dat'
    camb_matterpowerspectrum = '../../cases/Planck512/TF_Planck_BAO_matterpower_z0.dat'
    omega_m = 0.258
    lmax = 100
    h = 0.678
    x = 150/h #Mpc
    
    
    end = 596
    data = sp.loadtxt(camb_matterpowerspectrum)
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
    #Wk2 = 1 # for testing
    
    H = 100*h
    f = omega_m**0.6
    
    Pvk = H**2*f**2*ks**-2*Pmk*Wk2
    
    B = sp.zeros(lmax+1)
    C = sp.zeros(lmax)
    spherical_bessel_kx_l = sp.zeros((len(ks),lmax+1))
    special_spherical_bessel_kx_l = sp.zeros((len(ks),lmax+1))
    
    for row, k in enumerate(ks):
        spherical_bessel_kx_l[row] = sp.special.sph_jn(lmax,k*x)[0]
    #    factor = sp.sqrt(sp.pi/(2*k*x))
    #    for column, l in enumerate(range(lmax)):
    #        spherical_bessel_kx_l[row][column] = spherical_bessel(l,k*x)
    
        
    for l in range(lmax+1):
        integrand = ks**2/(2*sp.pi)**3*Pvk*spherical_bessel_kx_l[:,l]**2
        B[l] = 4*sp.pi*sp.trapz(integrand,ks)
        
    for l in range(1,lmax):
        C[l] = 4*sp.pi/(2*l+1)*(l*B[l-1]+(l+1)*B[l+1])
    
     
    #plt.plot(ks,Pmk)
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlabel('$k [Mpc^{-1}]$')
    #plt.ylabel('P(k) $[Mpc^{3}]$')
    #plt.axis([0.001, 10, 1e-4, 1e5])
    
#    plt.figure()
    ls = sp.array(range(1,lmax+1))
    scaled_Cl = sp.sqrt(ls*(ls+1)*C)
    plt.plot(ls,scaled_Cl)
#    plt.xscale('log')
#    plt.yscale('log')
    plt.title('Velocitypowerspectrum')
    plt.xlabel('$l$')
    plt.ylabel('$\sqrt{l(l+1)C_l}$ [km/s]')
    plt.xlim([1,20])

h = 0.678
velocitypowerspectrum(0)
velocitypowerspectrum(0.03) #Mpc. Roughly the diameter of the Milky Way
velocitypowerspectrum(4/h) #Mpc. 8 neighbours (from healpy.get_all_neighbours) and eq. 8 in paper 2
velocitypowerspectrum(150/h*0.13) # Roughly the diameter of the largest hole
