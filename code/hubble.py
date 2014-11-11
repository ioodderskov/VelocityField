from __future__ import division
import sys
import hubble_functions as hf
import hubble_classes as hc


# There is one argument, namely the parameterfile
if len(sys.argv) != 2:
    print "Wrong number of arguments"
    
parameterfile = sys.argv[1]
parameters = hc.Parameters(parameterfile)

halocatalogue = hf.load_halocatalogue(parameters)

halos = hf.initiate_halos(parameters,halocatalogue)
observers = hf.initiate_observers(parameters,halos)

for ob in observers:
    
    ob.observe(parameters,halos)
    if parameters.calculate_hubble_constants:
        ob.do_hubble_analysis(parameters)

hf.print_hubbleconstants_to_file(parameters,observers)