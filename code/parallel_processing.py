from __future__ import division


def observe_and_analyse(observer,parameters,particles):


    if parameters.use_lightcone:
        observer.observe(parameters,particles)
    elif parameters.calculate_pairwise_velocities:
        observer.calculate_pairwise_velocities(parameters)
    else:
        observer.observe(parameters,particles)

    
    if parameters.calculate_hubble_constants:
        observer.do_hubble_analysis(parameters)

    if parameters.calculate_powerspectra:
        observer.calculate_powerspectra(parameters)  
        
    return observer
        
