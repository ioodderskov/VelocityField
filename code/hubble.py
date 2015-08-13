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
import cPickle
import assignment_to_grid as ag
import resource
from scipy import stats # Import the scipy.stats module
from scipy.optimize import curve_fit # Import the curve fitting module
import matplotlib.pyplot as plt
import os
import scipy.integrate as integrate
import pdb
#import pickle

# There is one argument, namely the parameterfile
if len(sys.argv) != 2:
    print "Wrong number of arguments"
    
parameterfile = sys.argv[1]
parameters = hc.Parameters(parameterfile)


if parameters.data_type == "snapshot":
    hf.read_snapshot(parameters) # technically, its particles, not halos, in this case. But never mind.
    
else:
    halocatalogue = hf.load_halocatalogue(parameters,parameters.halocatalogue_file)
    hf.initiate_halos(parameters,halocatalogue)
    if parameters.use_HOD:
        hf.identify_galaxies(parameters)


if parameters.assign_to_grid:     
    ag.create_density_and_velocity_grid(parameters)
    
if parameters.data_to_observe == 'grid':
    hf.initiate_grid(parameters)

observers = hf.initiate_observers(parameters)

particles = []
if parameters.use_snapshot_for_background:
    particles = hf.read_snapshot(parameters)
#else:
#    particles = parameters.halos

partial_observe_and_analyse = partial(pp.observe_and_analyse,parameters=parameters,particles=particles)

if parameters.parallel_processing:
    pool = multiprocessing.Pool()
    observers = pool.map(partial_observe_and_analyse,observers)
    pool.close()
    pool.join()
else:
    observers = map(partial_observe_and_analyse,observers)

    
if parameters.calculate_hubble_constants:
    hf.print_hubbleconstants_to_file(parameters,observers)
    
#if parameters.calculate_powerspectra:
#    pf.print_powerspectra_to_file(parameters,observers)

if parameters.calculate_pairwise_velocities:
    hf.calculate_pairwise_velocities(parameters,observers)

if len(observers) <= 2000:
    sp.save(parameters.path+'observers.npy',observers)


if parameters.correct_for_peculiar_velocities:
    hf.calculate_velocity_correlation_coefficients(parameters,observers)


#parameters.halos = []
f = open(parameters.path+'parameters.save','w')    
cPickle.dump(parameters,f)
f.close()



    
