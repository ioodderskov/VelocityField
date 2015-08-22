from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import scipy.special


#Stuff to read in from the parameterfile
#camb_transfer_function = '../../cases/Planck512/TF_Planck_BAO_matterpower_z0.dat'
camb_matterpowerspectrum = '../../cases/Planck512/TF_Planck_BAO_matterpower_z0.dat'
omega_m = 0.25
lmax = 100
x = 150
R = 2 # omtrendt radius af Maelkevejen

data = sp.loadtxt(camb_matterpowerspectrum)
ks, Pmk = data[:,0], data[:,1]
#norm = TF[0]
#TF = TF/norm


#Pmk = ks*TF**2

H = 100
f = omega_m**0.6

C = sp.zeros(lmax)



for l in range(1,lmax):
    k_star = 2/sp.exp(1)*(1-sp.log(2)/(2*l))*l/x
    indices_term1 = (k_star < ks) & (ks < 1./R)
    indices_term2 = ks >= 1./R
#    if k_star < 1./R:
#        indices_term1 = (k_star < ks) & (ks < 1./R)
#        factor_term1 = 1
#    else:
#        indices_term1 = (1/R < ks) & (ks < k_star)
#        factor_term1 = -1
        
    integrand_term1 = 1./(ks[indices_term1]*x)**2*Pmk[indices_term1]
    integrand_term2 = 1./((ks[indices_term2]*x)**2*(ks[indices_term2]*R))*Pmk[indices_term2] 

    integral_term1 = sp.trapz(integrand_term1,ks[indices_term1])
    integral_term2 = sp.trapz(integrand_term2,ks[indices_term2])    
    C[l] = H**2*f**2/sp.pi*(integral_term1 + integral_term2)
#    
#for l in range(1,lmax):
#    C[l] = 4*sp.pi/(2*l+1)*(l*B[l-1]+(l+1)*B[l+1])
#
# 
##plt.plot(ks,Pmk)
##plt.xscale('log')
##plt.yscale('log')
##plt.xlabel('$k [Mpc^{-1}]$')
##plt.ylabel('P(k) $[Mpc^{3}]$')
##plt.axis([0.001, 10, 1e-4, 1e5])
#
ls = sp.array(range(1,lmax+1))
scaled_Cl = sp.sqrt(ls*(ls+1)*C)
plt.plot(ls,scaled_Cl)
plt.xscale('log')
plt.yscale('log')
#plt.axis([1,100,50,1000])
#plt.xlim([1,100])
#
