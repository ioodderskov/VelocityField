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




# There is one argument, namely the parameterfile
if len(sys.argv) != 2:
    print "Wrong number of arguments"
    
parameterfile = sys.argv[1]
parameters = hc.Parameters(parameterfile)


if parameters.snapshot:
    halos = hf.read_snapshot(parameters) # technically, its particles, not halos, in this case. But never mind.
    observers = hf.initiate_observers(parameters,[])

elif parameters.use_lightcone:
    observers = hf.initiate_observers(parameters,[])
    
else:
    halocatalogue = hf.load_halocatalogue(parameters,parameters.halocatalogue_file)
    halos = hf.initiate_halos(parameters,halocatalogue)
    observers = hf.initiate_observers(parameters,halos)

particles = []
if parameters.use_snapshot_for_background:
    particles = hf.read_snapshot(parameters)

partial_observe_and_analyse = partial(pp.observe_and_analyse,parameters=parameters,halos=halos,particles=particles)

if parameters.parallel_processing:
    pool = multiprocessing.Pool()
    observers = pool.map(partial_observe_and_analyse,observers)
    pool.close()
    pool.join()
else:
    observers = map(partial_observe_and_analyse,observers)


if parameters.assign_to_grid:    
    ag.create_density_and_velocity_grid(parameters,halos)


if parameters.calculate_hubble_constants:
    hf.print_hubbleconstants_to_file(parameters,observers)
    
if parameters.calculate_powerspectra:
    pf.print_powerspectra_to_file(parameters,observers)





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






#f = file(parameters.path+'parameters.save', 'wb')
#cPickle.dump(parameters, f, protocol=cPickle.HIGHEST_PROTOCOL)
#f.close()
#
#sp.save(parameters.path+'halos',halos)

