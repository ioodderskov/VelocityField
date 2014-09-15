from __future__ import division
import numpy as sp
import sys
sys.path.insert(0,'/home/io/Dropbox/Projekter/Hubble/VelocityField/SNdata')
#import SN_redshiftdistribution as SN
import hubble_functions as hf
import multiprocessing
from functools import partial


def distance(z):
    c = 3e5 # km/s
    H0 = 100 # km/s/Mpc/h
    d = z*c/H0 # Mpc/h
    return d

def mean_sys_error(h,Wz):
    dev = h-sp.ones(sp.shape(h));
    norm_Wz = Wz/sum(Wz)
    error = sp.sqrt(sp.sum(norm_Wz*dev**2)) 
    return error;
    


def error_for_observer(obs, observer_list, Wz, zbins, halo_list, boxsize, number_of_cones, skyfraction):   

    radial_distances_tot = list([])
    radial_velocities_tot = list([])
    bindistances = sp.zeros(2)
    Hubbleconstants = sp.NaN*sp.ones(len(zbins))
    
    
    
    for i in range(len(zbins)-1):
        zmin = zbins[i]
        zmax = zbins[i+1]
        number_of_SNe = int(Wz[i])
    
    
        bindistances[0] = distance(zmin)
        bindistances[1] = distance(zmax)
        sys.stdout.flush()
        if bindistances[1] > boxsize/2:
            break
    
         # Calculate the Hubbleconstants for all observers and all distances (or number of SNe)
        radial_distances, radial_velocities = hf.calculate_hubble_constants_for_all_observers(obs,observer_list, halo_list, number_of_SNe, bindistances, boxsize, number_of_cones, skyfraction)
        radial_distances_tot = radial_distances_tot+radial_distances
        radial_velocities_tot = radial_velocities_tot+radial_velocities
        
#        print "i=",i, ": H0 = ", observer_list[obs].Hubbleconstants[1], "zmin,zmax =", zmin,zmax,"dmin,dmax = ",bindistances[0],bindistances[1], "#SN = ",number_of_SNe
        Hubbleconstants[i+1] = observer_list[obs].Hubbleconstants[1]
        
    bindistances = distance(zbins)
    observer_list[obs].Hubbleconstants = Hubbleconstants
    # Print the results to a file
    #hf.print_hubbleconstants(hubblefile,bindistances,observer_list)
    
    
    Hubbleconstants[ sp.isnan(Hubbleconstants) ] = 100
    
    error = mean_sys_error(Hubbleconstants/100,bindistances)
    print "obs = ",obs,": The std of the theoretical probability distribution is", error*100, "%"
    
    return error


    
   
def calculate_mean_deviation_over_surveyrange(observer_list,Wz, zbins,  halo_list, boxsize, number_of_cones, skyfraction):    

    

    N_observers = len(observer_list)    
#    errors = sp.zeros(N_observers)
    print "The number of observers is", N_observers
    print "Wz = ", Wz
    print "zbins = ", zbins
    
    partial_error_for_observer = partial(error_for_observer,observer_list=observer_list, Wz=Wz, zbins=zbins, halo_list=halo_list, boxsize=boxsize, number_of_cones=number_of_cones, skyfraction=skyfraction)
    
    pool = multiprocessing.Pool()
    errors = pool.map(partial_error_for_observer,range(0,N_observers))
    pool.close()
    pool.join()
#    for obs in range(0,N_observers):
#        errors[obs] = 100*error_for_observer(obs, observer_list, Wz, zbins, halo_list, boxsize, number_of_cones, skyfraction) # %
    
    return sum(errors)*100/N_observers
