from __future__ import division
import scipy as sp
import random
import scipy.linalg as la
#import multiprocessing
#from functools import partial
#import pdb
import scipy.ndimage as ndi
import hubble_classes as hc
from scipy import interpolate
import pdb






        


# This function calculate the bin distances, so that every shell has the same volume
def calculate_bindistances(mind, maxd, width):
    mind_bin = mind    
    maxd_bin = mind+width   
    bindistances = [mind]
    while maxd_bin <= maxd:
        bindistances.append(maxd_bin)
        maxd_bin_tmp = (2*maxd_bin**3-mind_bin**3)**(1/3)
        mind_bin = maxd_bin
        maxd_bin = maxd_bin_tmp
    
    return bindistances
    
    
# This function calculates the local Hubble parameters
def calculate_H(rvsum,r2sum,halo_count):

    rvmean = rvsum/halo_count
    r2mean = r2sum/halo_count
    Hloc = rvmean/r2mean
    
    return Hloc
    
    
    
# Load the halo catalogue
def load_halocatalogue(halocatalogue_file):
    
    halocatalogue_unsorted = sp.loadtxt(halocatalogue_file)
    halocatalogue = sp.array(sorted(halocatalogue_unsorted, 
                                    key=lambda halocatalogue_unsorted: halocatalogue_unsorted[2]))
                                    
    return halocatalogue
    
    
def initiate_halos(parameters, halocatalogue):
    n_halos = len(halocatalogue)
    halos = [None]*n_halos
    
    for h in range(n_halos):
        position = halocatalogue[h,[8,9,10]]
        velocity = halocatalogue[h,[11,12,13]]
        mass = halocatalogue[h,2]
        ID = int(halocatalogue[h,0])
        ID_host = int(halocatalogue[h,33])
        
        halo = hc.Halo(position,velocity,mass,ID,ID_host,h)
        halos[h] = halo
        
    return halos
    

def initiate_observers(parameters,halos):
    
    if parameters.use_lightcone:
        observers = initiate_observers_from_file(parameters)
        
    else:
    
        if parameters.observer_choice == 'from_file':
            observers = initiate_observers_from_file(parameters)
    
        if parameters.observer_choice == 'subhalos':
            observers = initiate_observers_subhalos(parameters,halos)
    
        if parameters.observer_choice == 'random_halos':
            observers = initiate_observers_random_halos(parameters,halos)
    
        if parameters.observer_choice == 'random_positions':
            observers = initiate_observers_random_positions(parameters)
    
    
    return observers
    
    
def initiate_observers_from_file(parameters):
    observer_positions = sp.array(sp.loadtxt(parameters.observerfile))
    observers = [None]*len(observer_positions)
    
    for ob in range(len(observer_positions)):
        
        position = observer_positions[ob,[0,1,2]]/1000
        observers[ob] = hc.Observer(ob,position)
        
    return observers[0:parameters.number_of_observers]
    
    
def initiate_observers_random_halos(parameters,halos):
    random_indices = sp.random.random_integers(0,len(halos)-1,parameters.number_of_observers)
    observers = [None]*len(len(random_indices))
    
    for ob_index, ob_number in zip(random_indices,range(len(random_indices))):
        halo = halos[ob_index]
        position = halo.position
        
        observers[ob_number] = hc.Observer(position)
        
    return observers[0:parameters.number_of_observers]
    
def initiate_observers_subhalos(parameters,halos):
    masses = [halo.mass for halo in halos]
    
    localgroup_indices = (parameters.sub_min_m < masses) & (masses < parameters.sub_max_m)
    virgo_indices = (parameters.host_min_m < masses) & (masses < parameters.host_max_m)
    
    ID_hosts = sp.array([halo.ID_host for halo in halos])
    IDs = sp.array([halo.ID for halo in halos])
    virgo_IDs = IDs[virgo_indices]
    
    observer_indices = [localgroup_index for localgroup_index in localgroup_indices \
                        if ID_hosts[localgroup_index] in virgo_IDs]
                            
    observers = [None]*len(observer_indices)
    
    for ob_index, ob_number in zip(observer_indices, range(len(observer_indices))):
        halo = halos[ob_index]
        position = halo.position
        observers[ob_number] = hc.Observer(position)
    
    return observers[0:parameters.number_of_observers]

def initiate_observers_random_positions(parameters):

    sp.random.seed(0)
    observer_positions = parameters.boxsize*sp.rand(parameters.number_of_observers,3)
    
    observers = [None]*len(observer_positions)
        
    for ob in range(len(observer_positions)):
        position = [observer_positions[0,1,2]]
        observers[ob] = hc.Observer(position)
    
    
    return observers[0:parameters.number_of_observers]
    
    
def periodic_boundaries(parameters,xobs,yobs,zobs,xop,yop,zop):
    
    x,y,z = xop-xobs,yop-yobs,zop-zobs
    
    def periodic_coordinate(coordinate,parameters):
        
        if coordinate > parameters.boxsize/2:
            coordinate = coordinate-parameters.boxsize
        if coordinate < -parameters.boxsize/2:
            coordinate = coordinate+parameters.boxsize
            
        return coordinate
            
      
    x = periodic_coordinate(x,parameters)+xobs
    y = periodic_coordinate(y,parameters)+yobs    
    z = periodic_coordinate(z,parameters)+zobs
    
    return [x,y,z]
    
    
    
def spherical_coordinates(parameters,xobs,yobs,zobs,xop,yop,zop):

    rvec = sp.array([xop-xobs,yop-yobs,zop-zobs])
    r = la.norm(rvec)

    # Just to prevent problems if the observer is on top of the halo
    if r == 0:
        r = 1e-15         
        
    theta = sp.arccos((zop-zobs)/r)
    phi = sp.arctan2(yop-yobs,xop-xobs)+sp.pi

    return r,theta, phi
    
#    if parameters.distances_from_perturbed_metric:
#
#        # No reason to do a lengthy calculation of the perturbed distance, if 
#        # the halo won't be used anyway
#        max_distance_correction = sp.sqrt(1-2*parameters.potential_min)
#        if (ro*max_distance_correction < parameters.mind) or (ro > parameters.maxd):
#            return ro,theta
#        
#        rper = calculate_distance_in_perturbed_metric(parameters,xobs,yobs,zobs,xop,yop,zop,ro)
#        return rper, theta
#        
#    else:
     
    
    
def distance_correction_from_perturbed_metric(parameters,xobs,yobs,zobs,xop,yop,zop):
    
    res = 10
    f = sp.linspace(0,1,res)
    
    d = sp.array([f*(xop-xobs)+xobs, f*(yop-yobs)+yobs, f*(zop-zobs)+zobs])
#    pdb.set_trace()

#    d[d<0] = d[d<0]+parameters.boxsize
    Ng = len(parameters.potential)    
    d_grid = d*Ng/parameters.boxsize-1/2.
    # For the sake of the periodic boundaries
    d_grid[d_grid < 0] = d_grid[d_grid < 0]+Ng
    d_grid[d_grid > Ng] = d_grid[d_grid > Ng]-Ng
#    pdb.set_trace()
 
    psi = ndi.map_coordinates(parameters.potential,d_grid,mode='nearest')
#    pdb.set_trace()
    tck = interpolate.splrep(f,sp.sqrt(1-2*psi),s=0)
#    pdb.set_trace()
    psi_int = interpolate.splint(0,1,tck)
#    rper = psi_int*ro
#    pdb.set_trace()
    # This is only to check the interpolation of the potential onto the vector
#    sp.save("../cases/sim16/line_of_sight_vector",d*parameters.boxsize)
#    sp.save("../cases/sim16/psi_along_the_line_of_sight",psi)
#    sp.save("../cases/sim16/potential_array",parameters.potential)
#
#    import matplotlib.pyplot as plt
#    x = d[0,:]
##    x = d[0,:]*parameters.boxsize    
#    potential_along_x_axis = parameters.potential[:,0,0]
#    x_pot = sp.linspace(1,15,8)
#
#    plt.plot(x_pot,potential_along_x_axis)
#    plt.plot(x_pot-parameters.boxsize,potential_along_x_axis)
#    plt.plot(x_pot+parameters.boxsize,potential_along_x_axis)
#    plt.plot(x,psi)
#    plt.show()

    
    return psi_int

def select_candidates(parameters,halos,candidates):

    if parameters.observed_halos == 'all':
        selected_candidates = range(len(candidates))    

    if parameters.observed_halos == 'mass_weighted':
        selected_candidates = mass_weighted_selection_of_halos(parameters,halos,candidates)

    if parameters.observed_halos == 'random':
        selected_candidates = random_selection_of_halos(parameters,halos,candidates)
        
    return selected_candidates

   

    
def random_selection_of_halos(parameters,halos,candidates):
    
    selected_candidates = []
    random.seed(0)
    for n in range(parameters.number_of_SNe):
        rnd_candidate_number = random.randint(0,len(candidates)-1) # For some reason, the endpoint is included in the interval.
#        rnd_candidate_number = candidates[rnd_int]
        selected_candidates.append(rnd_candidate_number)
        
 
    return selected_candidates

def mass_weighted_selection_of_halos(parameters,halos,candidates):
    
    total_mass = sp.sum([c.mass for c in candidates])
    
    cum_mass = 0
    cum_masses = []
    
    for candidate in candidates:
        mass = candidate.mass
        cum_mass = cum_mass+mass
        cum_masses.append(cum_mass/total_mass)
        
    
    selected_candidates = []
    random.seed(0)
    for n in range(parameters.number_of_SNe):
        rnd = random.random()
        for candidate_number in range(len(candidates)):
            if cum_masses[candidate_number] >= rnd:
                selected_candidates.append(candidate_number)
                break
            
    return selected_candidates
        
        
def print_hubbleconstants_to_file(parameters,observers):

    if parameters.vary_skyfraction:
        
        for row, skyfraction in enumerate(parameters.skyfractions):
            
            f = open(parameters.hubblefile+str(skyfraction),'w')
            f.write("085\t")

            for bl in parameters.bindistances:
                f.write("%s\t" % bl)
                
            for ob_number, ob in enumerate(observers):
               f.write("\n%s\t" % ob_number)
               
               for b in range(1,len(parameters.bindistances)):
                   f.write("%s\t" % ob.Hubbleconstants[row][b])
                   
               f.write("\n")
               
            f.close()
            
    elif parameters.vary_number_of_SNe:

        f = open(parameters.hubblefile,'w')
        f.write("085\t")
        
        for number_of_SNe in parameters.numbers_of_SNe:
            f.write("%s\t" % number_of_SNe)
            
        for ob_number, ob in enumerate(observers):
            f.write("\n%s\t" % ob_number)
            
            for n in range(1,len(parameters.numbers_of_SNe)):
                f.write("%s\t" % ob.Hubbleconstants[n][-1])
                
            f.write("\n")
            
        f.close()


        
    else:
    
        f = open(parameters.hubblefile,'w')
        f.write("085\t")
        
        for bl in parameters.bindistances:
            f.write("%s\t" % bl)
            
        for ob_number, ob in enumerate(observers):
           f.write("\n%s\t" % ob_number)
           
           for b in range(1,len(parameters.bindistances)):
               f.write("%s\t" % ob.Hubbleconstants[b])
               
           f.write("\n")
           
        f.close()
    

    
    

        
        
        
def calculate_Hs_for_these_observed_halos(parameters,list_of_halos):
    
    rvsum = sp.zeros(len(parameters.bindistances))
    r2sum = sp.zeros(len(parameters.bindistances))
    halo_counts = sp.zeros(len(parameters.bindistances))
    
    for observed_halo in list_of_halos: 
        r = observed_halo.r
        vr = observed_halo.vr
        b = len(parameters.bindistances)-1
        rb = parameters.bindistances[b]
        
        while r < rb:
            
            rvsum[b] = rvsum[b]+r*vr
            r2sum[b] = r2sum[b]+r**2
            halo_counts[b] = halo_counts[b]+1
            b = b-1
            rb = parameters.bindistances[b]
            
    
    Hubbleconstants = [None]*len(parameters.bindistances)
    for b in range(len(parameters.bindistances)):
        
        if halo_counts[b] == 0:
            continue
        
        else:
            Hubbleconstants[b] = calculate_H(rvsum[b],r2sum[b],halo_counts[b])
        
    return Hubbleconstants
    

def calculate_Hs_for_varying_skyfractions(parameters,observed_halos):
    
      
    if parameters.number_of_cones == 1:
        theta_max_values = sp.arccos(1-2*parameters.skyfractions)
        
    if parameters.number_of_cones == 2:
        theta_max_values = sp.arccos(1-parameters.skyfractions)
        
    skyfraction_Hubbleconstants_array = sp.zeros((len(parameters.skyfractions),len(parameters.bindistances)))
    
    for row, theta_max in enumerate(theta_max_values):
        
        if parameters.number_of_cones == 1:
            observed_halos_in_cone = [oh for oh in observed_halos if oh.theta < theta_max]
            
        if parameters.number_of_cones == 2:
            observed_halos_in_cone = [oh for oh in observed_halos if oh.theta < theta_max or sp.pi-oh.theta < theta_max]
        
        skyfraction_Hubbleconstants_array[row] = calculate_Hs_for_these_observed_halos(parameters,observed_halos_in_cone)

        
        
    return skyfraction_Hubbleconstants_array
    
def calculate_Hs_for_varying_number_of_SNe(parameters,observed_halos):
    
    number_of_SNe_Hubbleconstants_array = sp.zeros((len(parameters.numbers_of_SNe),len(parameters.bindistances)))

    
    for row, number_of_SNe in enumerate(parameters.numbers_of_SNe):
        observed_sample_of_halos = random.sample(observed_halos,number_of_SNe)
        
        number_of_SNe_Hubbleconstants_array[row] = calculate_Hs_for_these_observed_halos(parameters,observed_sample_of_halos)
        
    return number_of_SNe_Hubbleconstants_array

    
    
   
                    
      
