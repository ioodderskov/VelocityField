from __future__ import division

def observe_and_analyse(observer,parameters,halos):


    if parameters.use_lightcone:
        observer.observe(parameters,[])
    else:
        observer.observe(parameters,halos)
    
    if parameters.calculate_hubble_constants:
        observer.do_hubble_analysis(parameters)

    if parameters.calculate_powerspectra:
        observer.calculate_powerspectra(parameters)  
        
    return observer
        