# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 13:50:55 2014

@author: io
"""
from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
#from gplot import Plot
from matplotlib.patches import Polygon

#plt = Plot('latex_full_wojtak')

def min_test(a):
    print("Nu er jeg i funktionen")
    return a



#x = sp.linspace(0,3,100)
#y = sp.sin(x)
#
#plt.plot(x,y)
#
#plt.savefig("min_sinus_figur.pdf ")

def load_Hubbleconstants(fil,a,N_observers):
    print "Loading Hubbleconstants"
    print fil
    data = sp.loadtxt(fil);
    rmax = data[0,1:];

    omegam = 0.306843; # Kun relevant for lightcone

        
    hubble = 100* sp.sqrt(omegam/(a*a*a) + 1. - omegam);


    Htot = data[1:,1:]

    # Frasorterer nan-results
    Htot_nonan = Htot[~sp.isnan(Htot).any(axis=1)]
    print "The total number of nan-results was", len(Htot)-len(Htot_nonan)

    H = Htot_nonan[0:N_observers,:]/hubble; # Indeholer H/H0 for alle observatører 

    return(rmax,H)
    
def mean_sys_error(h,rmax):
    dev = h-sp.ones(sp.shape(h));
    error = sp.sqrt(sp.sum(sp.sum(dev**2,axis=0))/(len(rmax)*len(h))) ;
    print error   
    return error;
    

def confidenslimit(h,conf):

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
                
        sigmaU = sigmaU+1e-4
        
    # Lower limit is calculated
    
    frac_of_deviators = 1;
    
    while frac_of_deviators > (1-conf)/2:
        deviators = dev < -sigmaL
        frac_of_deviators = 1.*sp.sum(deviators)/N
        
        sigmaL = sigmaL+1e-4
        
        
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

        # Første søjle indeholder alle de nedre grænser, anden søjle alle de øvre
        sigma68[b,:] = sl68, su68
        sigma95[b,:] = sl95, su95
        sigma99[b,:] = sl99, su99
     
    
    return(sigma68,sigma95,sigma99)


def mu_and_sigma(rmax,R,H):
    diff = sp.absolute(rmax-R)
    index = sp.where(diff == diff.min())
    dev1 = H[:,index]-1
    mu = dev1.mean()
    dev2 = (dev1-mu)**2
    sigma68 = sp.sqrt( dev2.mean() )
    # Convert to percent and round
    mu_pc = round(mu*100,1)
    sigma68_pc = round(sigma68*100,1)
    return(mu_pc, sigma68_pc)
    

def plot_patch(rmax,H,sigma68,sigma95,sigma99):
    print("Plotting confidence intervals")
    x = sp.append(rmax, sp.flipud(rmax))
    
    # Initiate figure
#    plt.clf()
#    fig = plt.figure(1)
    
    # Columnwise mean of H:
    Hmean = H.mean(axis=0)

    y68_l = -sigma68[:,0]+Hmean
    y68_u = sigma68[:,1]+Hmean
    y68 = sp.append(y68_l,sp.flipud(y68_u))
    
    y95_l = -sigma95[:,0]+Hmean
    y95_u = sigma95[:,1]+Hmean
    y95 = sp.append(y95_l,sp.flipud(y95_u))
    
    y99_l = -sigma99[:,0]+Hmean
    y99_u = sigma99[:,1]+Hmean
    y99 = sp.append(y99_l,sp.flipud(y99_u))
    
    p68=Polygon(zip(x,y68),alpha=0.3,lw=0)
    p95=Polygon(zip(x,y95),alpha=0.2,lw=0)
    p99=Polygon(zip(x,y99),alpha=0.1,lw=0)
    
#    plt.figure(1)
#    plt.hold(True)

#    plt.xlabel('$r_{max} [Mpc/h]$')
#    plt.axis([64, 256, 0.90, 1.10])    
    
    plt.gca().add_artist(p68)
    plt.gca().add_artist(p95)
    plt.gca().add_artist(p99)
    plt.plot(rmax,Hmean)
    

    
    return(0)
    
    
def plot_lines(rmax,H,sigma68,sigma95,sigma99):
    
    print("Plotting lines for comparison")

    # Columnwise mean of H:
    Hmean = H.mean(axis=0) 
    
    
    y68_l = -sigma68[:,0]+Hmean
    y68_u = sigma68[:,1]+Hmean
    
    y95_l = -sigma95[:,0]+Hmean
    y95_u = sigma95[:,1]+Hmean
    
    y99_l = -sigma99[:,0]+Hmean
    y99_u = sigma99[:,1]+Hmean

    
    for y in [Hmean, y68_l, y68_u, y95_l, y95_u, y99_l, y99_u]:
        plt.plot(rmax,y,'k:')
        
#    plt.plot(rmax,Hmean,':')
    