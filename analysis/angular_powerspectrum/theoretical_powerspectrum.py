from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import scipy.special


#Stuff to read in from the parameterfile
camb_transfer_function = '../../cases/Planck512/TF_Planck_BAO_matterpower_z0.dat'
omega_m = 0.3
lmax = 100
x = 150


data = sp.loadtxt(camb_transfer_function)
ks, TF = data[:,0]*0.678, data[:,1]
norm = TF[0]
TF = TF/norm

Pmk = ks*TF**2

Wk2 = 1 # for testing

H = 100
f = omega_m**0.6

Pvk = H**2*f**2*ks**-2*Pmk*Wk2

B = sp.zeros(lmax+1)
C = sp.zeros(lmax)
spherical_bessel_kx_l = sp.zeros((len(ks),lmax+1))

for row, k in enumerate(ks):
    spherical_bessel_kx_l[row] = sp.special.sph_jn(lmax,k*x)[0]
    
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

ls = sp.array(range(1,lmax+1))
scaled_Cl = sp.sqrt(ls*(ls+1)*C)
plt.plot(ls,scaled_Cl)
plt.xscale('log')
plt.yscale('log')
#plt.xlim([1,100])