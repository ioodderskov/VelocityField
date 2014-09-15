from __future__ import division
import scipy as sp
import sys
sys.path.insert(0,'/home/io/Dropbox/Projekter/Hubble/VelocityField/SNdata')
import SN_redshiftdistribution as SN
import hubble_functions as hf


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
    
   
def calculate_mean_deviation_over_surveyrange(observer_list, halo_list, boxsize, number_of_cones, skyfraction):    

    Wz, zbins = SN.get_z_distribution()
    #Wz, zbins = sp.array([13,30,22,26,26,23,10,5,4,3,3,2,3,2,0,1,0,0]), sp.arange(0.01,0.1,0.005)
    #Wz = Wz*100
    

    N_observers = len(observer_list)    
    errors = sp.zeros(N_observers)
    print "The number of observers is", N_observers
    print "Wz = ", Wz
    print "zbins = ", zbins
    
    sys.stdout.flush()
    
    for obs in range(0,N_observers):
        
        print "obs = ", obs
        
        sys.stdout.flush()
    
        # Calculate the bin-distances
        #bindistances = hf.calculate_bindistances(mind, maxd, width)
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
            print bindistances[1]
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
        errors[obs] = error*100 # %
        print "obs = ",obs,": The std of the theoretical probability distribution is", error*100, "%"
        sys.stdout.flush()
    
    return sum(errors)/N_observers
