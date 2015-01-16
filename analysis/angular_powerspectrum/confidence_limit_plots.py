from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import pdb


def confidenslimit(h,conf):

#    delta = 1e-4
    delta = 1
    N = len(h)
    middel = h.mean(axis=0)
    
    sigmaU = 0;
    sigmaL = 0;
    
    # How big is the deviation for each observer?
    dev = h-middel;
    
    # If the allowed deviation is zero, everybody are deviators.
    frac_of_deviators = 1;
    
    # Now the allowed deviation is increased until a percentage equal to "conf"
    # of the observers are measuring a deviation less then or equal to the one allowed.
    # Half of the deviators will be above, and half below the allowed values.
    
    # Upper limit is calculated
    while frac_of_deviators > (1-conf)/2:
        deviators = dev > sigmaU;
        frac_of_deviators = 1.*sp.sum(deviators)/N
                
        sigmaU = sigmaU+delta
        
    # Lower limit is calculated
    
    frac_of_deviators = 1;
    
    while frac_of_deviators > (1-conf)/2:
        deviators = dev < -sigmaL
        frac_of_deviators = 1.*sp.sum(deviators)/N
        
        sigmaL = sigmaL+delta
        
        
    return(sigmaL,sigmaU)



def calc_confidence_intervals(H):

    number_of_bins = len(H[1,:]);

    sigma68 = sp.zeros((number_of_bins,2));
    sigma95 = sp.zeros((number_of_bins,2));
    sigma99 = sp.zeros((number_of_bins,2));
    
    for b in range(0,number_of_bins):
        
        sl68, su68 = confidenslimit(sp.take(H,[b],axis=1),0.683);
        sl95, su95 = confidenslimit(sp.take(H,[b],axis=1),0.954);
        sl99, su99 = confidenslimit(sp.take(H,[b],axis=1),0.997);

        # Foerste soejle indeholder alle de nedre graenser, anden soejle alle de oevre
        sigma68[b,:] = sl68, su68
        sigma95[b,:] = sl95, su95
        sigma99[b,:] = sl99, su99
     
    
    return(sigma68,sigma95,sigma99)

def scale(ls,Cls):
    return sp.array(sp.sqrt(Cls*ls*(ls+1)))    
    
def plot_patch(ls,H,sigma68,sigma95,sigma99):
    print("Plotting confidence intervals")
    x = sp.append(ls, sp.flipud(ls))
    
    # Initiate figure
#    plt.clf()
#    fig = plt.figure(1)
    
    # Columnwise mean of H:
    Hmean = H.mean(axis=0)

    y68_l = scale(ls,-sigma68[:,0]+Hmean)
    y68_u = scale(ls,sigma68[:,1]+Hmean)
    y68 = sp.append(y68_l,sp.flipud(y68_u))
    
    y95_l = scale(ls,-sigma95[:,0]+Hmean)
    y95_u = scale(ls,sigma95[:,1]+Hmean)
    y95 = sp.append(y95_l,sp.flipud(y95_u))
    
    y99_l = scale(ls,-sigma99[:,0]+Hmean)
    y99_u = scale(ls,sigma99[:,1]+Hmean)
    y99 = sp.append(y99_l,sp.flipud(y99_u))
    
    p68=Polygon(zip(x,y68),alpha=0.3,lw=0)
    p95=Polygon(zip(x,y95),alpha=0.2,lw=0)
    p99=Polygon(zip(x,y99),alpha=0.1,lw=0)
    
#    plt.figure(1)
#    plt.hold(True)

#    plt.xlabel('$r_{max} [Mpc/h]$')
#    plt.axis([64, 256, 0.90, 1.10])    


    Hmean = scale(ls,Hmean)
    
    plt.gca().add_artist(p68)
    plt.gca().add_artist(p95)
    plt.gca().add_artist(p99)
    plt.plot(ls,Hmean)

    
    return(0)



case = 'Planck512'

data = sp.loadtxt('../../cases/' + case + '/powerspectra_lightcone.txt')
ls = data[0,2:]
Cls = data[1:,2:]

sigma68, sigma95, sigma99 = calc_confidence_intervals(Cls)

mean_Cls = sp.mean(Cls,axis=0)
scaled_Cls = scale(ls,mean_Cls)
plt.figure()
#plt.plot(ls,scaled_Cls)
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel('l')
plt.ylabel('$\sqrt{C_{l}\cdot l(l+1)}$ [km/s]')
plt.xlim([0,20])
plt.ylim([0,800])   

#scaled_sigma68 = sp.array(sp.sqrt(sigma68[:,0]*ls*(ls+1)))],[sp.array(sp.sqrt(sigma68[:,1]*ls*(ls+1)))]
plot_patch(ls,Cls,sigma68,sigma95,sigma99)
plt.savefig('../../cases/' + case + '/powerspectrum_lightcone.pdf')