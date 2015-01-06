from __future__ import division
import sys
import hubble_functions as hf
import hubble_classes as hc
import powerspectrum_functions as pf


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

    if parameters.calculate_powerspectra:
        ob.calculate_powerspectra(parameters)        

if parameters.calculate_hubble_constants:
    hf.print_hubbleconstants_to_file(parameters,observers)
    
if parameters.calculate_powerspectra:
    pf.print_powerspectra_to_file(parameters,observers)