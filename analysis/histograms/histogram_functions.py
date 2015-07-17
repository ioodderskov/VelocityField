# -*- coding: utf-8 -*- 

from __future__ import division
import scipy as sp
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import pdb
import matplotlib.pyplot as plt

def make_normalised_histogram(nbins,data):
#
#    weights = sp.ones_like(data)/len(data)
#    histogram = sp.histogram(data,bins=nbins,weights=weights)
    histogram = sp.histogram(data,bins=nbins,normed=True)
    fractions_in_bins = histogram[0]
    bin_locations = histogram[1]
    bin_centers = (bin_locations[:-1]+bin_locations[1:])/2
    
    bins = bin_centers
    bars = fractions_in_bins
#    print bins, bars
    
    return bins, bars
    
def fit_distribution_to_normalised_histogram(bins,bars,f):
    
#    dist = getattr(stats,dist_name)
#    pdb.set_trace()
    p0 = [-0.1,0.03]
#    plt.plot(bins,f(bins,*p0),'r')
#    print bins, bars
    popt, pcov = curve_fit(f,bins,bars,p0=p0)
    print "p0=", p0
    print "popt=", popt
    mu = sp.exp(popt[0])
    sigma = popt[1]*mu
    print "for x: mu = ", mu , "sigma = ", sigma
#    pdb.set_trace()
    
    bins_fine = sp.linspace(bins[0],bins[-1],100)
    
    # I make a finer spacing of the bins in order to make the fits look smooth
    fitted_bars = f(bins_fine,*popt)
    
    deviation = 8.8
    
    print "The probability of measuring a Hubble constant", deviation,"% or more too low (high) is",
    prob_low = integrate.quad(lambda x: f(x,*popt), 0, 1-deviation/100.)[0]
    prob_high = integrate.quad(lambda x: f(x,*popt), 1+deviation/100,20)[0]
    print prob_low*100,"%",
    print "(",prob_high*100,"%)"
    
    return bins_fine, fitted_bars
    
    
#På den måde kan man i curvefit skrive
#curve_fit(f,,xdata, ydata, p)
#og når man skal plotte det kan man så skrive:
#p_optimized = 
#
#curve_fit(f,,xdata, ydata, p)
#f(t, *p_optimized)
#Det er i hvert fald noget i den retning jeg bruger.
#Så er det eneste man skal holde styr på sit første gæt.