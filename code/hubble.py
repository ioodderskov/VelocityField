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
    print "observers = ", observers

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

if parameters.observer_choice == 'all':
    print "Not saving the observers - there are too many!"
else:
    sp.save(parameters.path+'observers.npy',observers)

if parameters.calculate_pairwise_velocities:

    pairwise_velocities = []
    radial_distances = []
    pair_masses = []
    for observer in observers:
        for velocity, radial_distance, pair_mass in zip(observer.observed_radial_peculiar_velocities,observer.radial_distances,observer.pair_masses):
            pairwise_velocities.append(velocity)
            radial_distances.append(radial_distance)
            pair_masses.append(pair_mass)
    
    sp.save(parameters.path+'pairwise_velocities.npy',sp.array(pairwise_velocities))
    sp.save(parameters.path+'radial_distances.npy',sp.array(radial_distances))
    sp.save(parameters.path+'pair_masses.npy',sp.array(pair_masses))

if parameters.correct_for_peculiar_velocities:

    if parameters.use_local_velocity:

        print "------ observers --------"
        local_velocities = sp.array([observer.local_velocity for observer in observers if len(observer.local_velocity) != 0] )
        local_velocity_corrections = sp.array([observer.local_velocity_correction for observer in observers if len(observer.local_velocity) != 0])
        
        
        print "correlation coefficients:"
        print "x:", sp.corrcoef(local_velocities[:,0],local_velocity_corrections[:,0])
        print "y:", sp.corrcoef(local_velocities[:,1],local_velocity_corrections[:,1])
        print "z:", sp.corrcoef(local_velocities[:,2],local_velocity_corrections[:,2])
        print "abs(local_velocities)/abs(local_velocity_corrections)",\
        sp.mean(sp.absolute(local_velocities),axis=0)/sp.mean(sp.absolute(local_velocity_corrections),axis=0)
        
    
    print "------ observed halos --------"
    observed_velocities = sp.array([observer.chosen_halos[0].observed_velocity for observer in observers if len(observer.chosen_halos) == 1])
    velocity_corrections = sp.array([observer.chosen_halos[0].velocity_correction for observer in observers if len(observer.chosen_halos) == 1])
    print "correlation coefficients:"
    print "x:", sp.corrcoef(observed_velocities[:,0],velocity_corrections[:,0])
    print "y:", sp.corrcoef(observed_velocities[:,1],velocity_corrections[:,1])
    print "z:", sp.corrcoef(observed_velocities[:,2],velocity_corrections[:,2])
    print "abs(observed_velocities)/abs(velocity_corrections)",\
    sp.mean(sp.absolute(observed_velocities),axis=0)/sp.mean(sp.absolute(velocity_corrections),axis=0)
