from __future__ import division
import sys
import hubble_functions as hf
import hubble_classes as hc
import powerspectrum_functions as pf
import parallel_processing as pp
import multiprocessing
from functools import partial

# There is one argument, namely the parameterfile
if len(sys.argv) != 2:
    print "Wrong number of arguments"
    
parameterfile = sys.argv[1]
parameters = hc.Parameters(parameterfile)

halocatalogue = hf.load_halocatalogue(parameters)

halos = hf.initiate_halos(parameters,halocatalogue)
observers = hf.initiate_observers(parameters,halos)

partial_observe_and_analyse = partial(pp.observe_and_analyse,parameters=parameters,halos=halos)

if parameters.parallel_processing:
    pool = multiprocessing.Pool()
    observers = pool.map(partial_observe_and_analyse,observers)
    pool.close()
    pool.join()
else:
    observers = map(partial_observe_and_analyse,observers)

 

if parameters.calculate_hubble_constants:
    hf.print_hubbleconstants_to_file(parameters,observers)
    
if parameters.calculate_powerspectra:
    pf.print_powerspectra_to_file(parameters,observers)