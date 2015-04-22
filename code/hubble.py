from __future__ import division
#from guppy import hpy
import sys
import hubble_functions as hf
import hubble_classes as hc
import powerspectrum_functions as pf
import parallel_processing as pp
import multiprocessing
from functools import partial
import scipy as sp
#import cPickle
import assignment_to_grid as ag
import resource
from scipy import stats # Import the scipy.stats module
from scipy.optimize import curve_fit # Import the curve fitting module
import matplotlib.pyplot as plt
import os
import scipy.integrate as integrate



print "Before beginning. The total memory used is:" 
print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000, 'Mb'


# There is one argument, namely the parameterfile
if len(sys.argv) != 2:
    print "Wrong number of arguments"
    
parameterfile = sys.argv[1]
parameters = hc.Parameters(parameterfile)


if parameters.snapshot:
    hf.read_snapshot(parameters) # technically, its particles, not halos, in this case. But never mind.
    observers = hf.initiate_observers(parameters,[])

elif parameters.use_lightcone:
    observers = hf.initiate_observers(parameters,[])
    
else:
    halocatalogue = hf.load_halocatalogue(parameters,parameters.halocatalogue_file)
    hf.initiate_halos(parameters,halocatalogue)
    observers = hf.initiate_observers(parameters)

particles = []
if parameters.use_snapshot_for_background:
    particles = hf.read_snapshot(parameters)

partial_observe_and_analyse = partial(pp.observe_and_analyse,parameters=parameters,particles=particles)

if parameters.parallel_processing:
    pool = multiprocessing.Pool(maxtasksperchild=32)
    observers = pool.map(partial_observe_and_analyse,observers)
    pool.close()
    pool.join()
else:
    observers = map(partial_observe_and_analyse,observers)


if parameters.assign_to_grid:    
    ag.create_density_and_velocity_grid(parameters)


if parameters.calculate_hubble_constants:
    hf.print_hubbleconstants_to_file(parameters,observers)
    
if parameters.calculate_powerspectra:
    pf.print_powerspectra_to_file(parameters,observers)

sp.save(parameters.path+'observers.npy',observers)

if parameters.correct_for_peculiar_velocities:

    print "------ observers --------"
    local_velocities = sp.array([observer.local_velocity for observer in observers if len(observer.local_velocity) != 0] )
    local_velocity_corrections = sp.array([observer.local_velocity_correction for observer in observers if len(observer.local_velocity) != 0])
    
    
    print "correlation coefficients:"
    print "x:", sp.corrcoef(local_velocities[:,0],local_velocity_corrections[:,0])
    print "y:", sp.corrcoef(local_velocities[:,1],local_velocity_corrections[:,1])
    print "z:", sp.corrcoef(local_velocities[:,2],local_velocity_corrections[:,2])
    print "abs(local_velocities)/abs(local_velocity_corrections)",\
    sp.mean(sp.absolute(local_velocities),axis=0)/sp.mean(sp.absolute(local_velocity_corrections),axis=0)
    
    bulk_flows = sp.array([observer.local_velocity-observer.local_velocity_correction for observer in observers])
    
    print "------ observed halos --------"
    observed_velocities = sp.array([observer.chosen_halos[0].observed_velocity for observer in observers if len(observer.chosen_halos) == 1])
    velocity_corrections = sp.array([observer.chosen_halos[0].velocity_correction for observer in observers if len(observer.chosen_halos) == 1])
    print "correlation coefficients:"
    print "x:", sp.corrcoef(observed_velocities[:,0],velocity_corrections[:,0])
    print "y:", sp.corrcoef(observed_velocities[:,1],velocity_corrections[:,1])
    print "z:", sp.corrcoef(observed_velocities[:,2],velocity_corrections[:,2])
    print "abs(observed_velocities)/abs(velocity_corrections)",\
    sp.mean(sp.absolute(observed_velocities),axis=0)/sp.mean(sp.absolute(velocity_corrections),axis=0)
    
    
    print "The total memory used is:" 
    print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000, 'Mb'

if parameters.calculate_pairwise_velocities:

    data = []
    for observer in observers:
        for velocity in observer.observed_radial_peculiar_velocities:
            data.append(velocity)

#        data = sp.array([velocity for velocity in observer.observed_radial_peculiar_velocities for observer in observers])
 
    data = sp.array(data)
    cut_tails = sp.absolute(data) < 200 
    data_cut = data[cut_tails]

    number_of_bins = 20
    plt.figure(figsize=(6,5))
    ax1=plt.subplot(111)
    
    dist = getattr(stats,'norm')
    
    xx = data
    hist, bin_edges = sp.histogram(xx, bins=number_of_bins, density=True) # Calculate histogram
    x_hist = bin_edges[1:] # Take the upper edges of the bins
    y_hist = hist.cumsum()/hist.cumsum().max()  # Normalise the cumulative sum   

    histogram = sp.histogram(data,bins=number_of_bins,normed=True)
    numbers_in_bins = histogram[0]
    fractions_in_bins = numbers_in_bins/sp.sum(numbers_in_bins)
    bin_locations = histogram[1]
    bin_width = bin_locations[1]-bin_locations[0]
    bin_centers = (bin_locations[:-1]+bin_locations[1:])/2
    factor = sp.sum(numbers_in_bins)

    plt.bar(bin_centers, numbers_in_bins/factor, align = 'center', width=bin_width, alpha=0.5, color='g')
    
    def cdf(xdata,loc,scale):
        return stats.norm.cdf(xdata,loc=loc,scale=scale)
    def pdf(xdata,loc,scale):
        return stats.norm.pdf(xdata,loc=loc,scale=scale)

        
    (mu_cdf, sigma_cdf), pcov = curve_fit(cdf, x_hist, y_hist,p0=[0,100])
    (mu_pdf, sigma_pdf), pcov = curve_fit(pdf, bin_centers, numbers_in_bins,p0=[0,100])    
    mu, sigma = dist.fit(data)
    mu_cut, sigma_cut = dist.fit(data_cut)
    
    x_resolution = 300
   
    x = sp.linspace(-500, 500, num=x_resolution) # values for x-axis
    fitted_distribution = dist.pdf(x, loc=mu,scale=sigma)
    fitted_distribution_pdf = dist.pdf(x, loc=mu_pdf,scale=sigma_pdf)
    fitted_distribution_cdf = dist.pdf(x, loc=mu_cdf,scale=sigma_cdf)
    fitted_distribution_cut = dist.pdf(x, loc=mu_cut,scale=sigma_cut)
    
    ax1.plot(x,fitted_distribution/factor, 'b', lw=2,
             label='All data: mu=%0.1f, sigma=%0.1f' %(mu,sigma))
    ax1.plot(x,fitted_distribution_pdf/factor, 'r', lw=2,
             label='Fit to pdf: mu=%0.1f, sigma=%0.1f' %(mu_pdf,sigma_pdf))
    ax1.plot(x,fitted_distribution_cdf/factor, 'black', lw=2,
             label='Fit to cdf: mu=%0.1f, sigma=%0.1f' %(mu_cdf,sigma_cdf))
    ax1.plot(x,fitted_distribution_cut/factor, 'y', lw=2,
             label='Cut data: mu=%0.1f, sigma=%0.1f' %(mu_cut,sigma_cut))
    plt.legend()

    plt.xlabel('Pairwise velocities',fontsize=16)
    plt.ylabel('P(pairwise velocity)',fontsize=16)

    plt.plot([1,1],[0,0.4],'k--',linewidth=1.5)
#    print "shape = ", shape_out
    
    def fitted_norm(x,mu,sigma):
        return 1./(sp.sqrt(2*sp.pi)*sigma)*sp.exp(-(x-mu)**2/(2*sigma**2))

#    sigma = shape_out
    
#    print "mu, sigma = ", mu,sigma
##    print "exp(mu) = ", sp.exp(mu)
#
    print "Checking that the histogram is normalised"
    print "sum(fractions_in_bins) = ", sp.sum(fractions_in_bins)
    
    print "Checking that the fitted distributions are normalised"
    mus = [mu, mu_pdf, mu_cdf, mu_cut]
    sigmas = [sigma, sigma_pdf, sigma_cdf, sigma_cut]
    for m,s in zip(mus,sigmas):
        print integrate.quad(lambda x: fitted_norm(x,m,s),-500,500)

#    sp.save(parameters.path+'pairwise_velocities.npy',data)
#    plt.title("mu = %f and sigma = %f" % (mu,sigma))
    plt.show()
#    plt.savefig(parameters.path+'pairwise_velocity_distribution.pdf')

