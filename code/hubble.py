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
    pool = multiprocessing.Pool()
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
    
    ## FIT THE DISTRIBUTION
#    (loc_out, scale_out), pcov = curve_fit(
#                lambda bin_centers, loc, scale: dist.pdf(bin_centers, loc, scale=scale),
#                bin_centers, numbers_in_bins/factor)         
#
#    mu = loc_out
#    sigma = scale_out
    
    mu, sigma = dist.fit(data)
    
    x_resolution = 300
   
    x = sp.linspace(-500, 500, num=x_resolution) # values for x-axis
    fitted_distribution = dist.pdf(x, loc=mu,scale=sigma)
    ax1.plot(x,fitted_distribution/factor, 'b', lw=2, label='Fitted distribution')

    plt.xlabel('Pairwise velocities',fontsize=16)
    plt.ylabel('P(pairwise velocity)',fontsize=16)

    plt.plot([1,1],[0,0.3],'k--',linewidth=1.5)
#    print "shape = ", shape_out
    
    def fitted_norm(x,mu,sigma):
        return 1./(sp.sqrt(2*sp.pi)*sigma)*sp.exp(-(x-mu)**2/(2*sigma**2))

#    sigma = shape_out
    
    print "mu, sigma = ", mu,sigma
#    print "exp(mu) = ", sp.exp(mu)
    
    print "Checking that the fitted distribution is normalised"
    print integrate.quad(lambda x: fitted_norm(x,mu,sigma),-500,500)
    print "Checking that the histogram is normed"
    print "sum(fractions_in_bins) = ", sp.sum(fractions_in_bins)
    sp.save(parameters.path+'pairwise_velocities.npy',data)
    plt.title("mu = %f and sigma = %f" % (mu,sigma))

    plt.savefig(parameters.path+'pairwise_velocity_distribution.pdf')


