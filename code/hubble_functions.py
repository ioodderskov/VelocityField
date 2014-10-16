from __future__ import division
import scipy as sp
import random
import scipy.linalg as la
import multiprocessing
from functools import partial
import pdb
import scipy.ndimage as ndi





class Halos:
    def __init__(self,x,y,z,vx,vy,vz,mass,ID,ID_host):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.mass = mass
        self.ID = ID
        self.ID_host = ID_host

    
class Observers:
    def __init__(self,x,y,z,selected_halos,Hubbleconstants,thetas,phis,rs,vrs):
        self.x = x
        self.y = y
        self.z = z
        self.selected_halos = selected_halos
        self.Hubbleconstants = Hubbleconstants
        self.thetas = thetas
        self.phis = phis
        self.rs = rs
        self.vrs = vrs
        


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
    
    
    
    
    
################## Read the halo data ##############################################
def read_halo_file(mass_sorted_data):

    n_halos = len(mass_sorted_data)
    halo_list = [None]*n_halos
    
    for i in range(n_halos):
        [x,y,z] = mass_sorted_data[i,[8,9,10]]
        [vx,vy,vz] = mass_sorted_data[i,[11,12,13]]
        mass = mass_sorted_data[i,2]
        ID = int(mass_sorted_data[i,0])
        ID_host = int(mass_sorted_data[i,33])
        halo = Halos(x,y,z,vx,vy,vz,mass,ID,ID_host)
        halo_list[i] = halo
    
    return halo_list

######################################################################################




######## Calculating the Hubbleconstants for all observers ########################



# This function calculates the Hubble constants for a given observer
def find_hubble_constants_for_observer(observer_number,observer_list,halo_list,observed_halos,bindistances,boxsize,number_of_cones,skyfraction,distances_from_perturbed_metric, potential_file):

    observer = observer_list[observer_number]
    [x,y,z] = [observer.x, observer.y, observer.z]
    mind = bindistances[0]
    maxd = bindistances[-1]
    selected_halos = select_halos(x,y,z,halo_list,mind,maxd,observed_halos,boxsize,number_of_cones,skyfraction, distances_from_perturbed_metric, potential_file)
    Hubbleconstants, radial_distances, radial_velocities = Hubble(x,y,z,halo_list, selected_halos, bindistances, boxsize, distances_from_perturbed_metric, potential_file)

    return Hubbleconstants, radial_distances, radial_velocities



def calculate_hubble_constants_for_all_observers(obs,observer_list, halo_list, observed_halos, bindistances, boxsize, number_of_cones, skyfraction, distances_from_perturbed_metric, potential_file):
    if isinstance(obs,int):    
        obs = list([obs])


    partial_find_hubble_constants_for_observer = partial(find_hubble_constants_for_observer,observer_list=observer_list,halo_list=halo_list,observed_halos=observed_halos,bindistances=bindistances,boxsize=boxsize,number_of_cones=number_of_cones,skyfraction=skyfraction, distances_from_perturbed_metric=distances_from_perturbed_metric, potential_file=potential_file)
#    pdb.set_trace()    
    pool = multiprocessing.Pool()
    out = pool.map(partial_find_hubble_constants_for_observer,obs)
#    out = map(partial_find_hubble_constants_for_observer,obs)
    pool.close()
    pool.join()

    # Writing the Hubbleconstants to the observer_list. For some reason, it does not work to do so
    # using parallel processing :-/
    for observer_number in obs:
#        observer.selected_halos = selected_halos
        observer = observer_list[observer_number]
        Hubbleconstants = out[observer_number][0]
        observer.Hubbleconstants = Hubbleconstants
    
    radial_distances = sp.array(out[-1][1])
    radial_velocities = sp.array(out[-1][2])
    return radial_distances, radial_velocities
        
##############################################################################



    
    
    
    
    


######################### Identification of observers #########################


#This function finds the halos that are have the characteristics specified in the parameterfile
#Or, if find_observers = 0, the positions are simple read from a file.
def find_observers(observer_choice,number_of_observers,number_of_SNe,boxsize,observerfile,halo_list,sub_min_m,sub_max_m,host_min_m,host_max_m):
    
    dummy = sp.zeros(number_of_SNe)
   
   ###################### Reading observers from file ################################
    if observer_choice == 'from_file':
        observer_positions = sp.loadtxt(observerfile)  

        observer_list = [None]*len(observer_positions)
    
        # Reading coordinates of the observers, and saving them in an observer list    
        for observer_number in range(len(observer_list)):
            [x,y,z] = [observer_positions[observer_number,0],observer_positions[observer_number,1],observer_positions[observer_number,2]]
            [x,y,z] = sp.array([x,y,z])/1000
            observer = Observers(x,y,z,dummy,dummy,dummy,dummy,dummy,dummy)
            observer_list[observer_number] = observer
   ########################################################################################                    
        



    
   ###################### Identifying subhalos as observers ################################
    if observer_choice == 'subhalos':

        masses = [halo.mass for halo in halo_list]
        
        # Identifying halos with masses corresponding to the Virgo Supercluster or the Local Group
        localgroup_indices = (sub_min_m < masses) & (masses < sub_max_m)
        virgo_indices = (host_min_m < masses) & (masses < host_max_m)

        #
        ID_hosts = sp.array([halo.ID_host for halo in halo_list])
        IDs = sp.array([halo.ID for halo in halo_list])
        virgo_IDs = IDs[virgo_indices]
    
        observers = [i for i, elem in enumerate(ID_hosts[localgroup_indices]) if elem in virgo_IDs]
        all_indices = sp.array(range(len(halo_list)))
        observer_indices = sp.array(all_indices[localgroup_indices])[observers]
        
        observer_list = [None]*len(observer_indices)
   
        # Reading coordinates of the observers, and saving them in an observer list    
        for halo_index,observer_number in zip(observer_indices,range(len(observer_indices))):
            observer = halo_list[halo_index]
            [x,y,z] = [observer.x, observer.y, observer.z]
            observer = Observers(x,y,z,dummy,dummy,dummy,dummy,dummy,dummy)
            observer_list[observer_number] = observer
   ########################################################################################            
            
            


            
   ###################### Choosing random positions as observers ################################         
    if observer_choice == 'random_positions':
        sp.random.seed(0)
        observer_positions = boxsize*sp.rand(number_of_observers,3)
        observer_list = [None]*len(observer_positions)
        for observer_number in range(number_of_observers):
            [x,y,z] = [observer_positions[observer_number,0],observer_positions[observer_number,1],observer_positions[observer_number,2]]
            observer = Observers(x,y,z,dummy,dummy,dummy,dummy,dummy,dummy)
            observer_list[observer_number] = observer
   ########################################################################################            





            
            
   ###################### Choosing random halos as observers ################################         
    if observer_choice == 'random_halos':
        random_halo_indices = sp.random.random_integers(0,len(halo_list)-1,number_of_observers)
        observer_list = [None]*len(random_halo_indices)
        
        for halo_index,observer_number in zip(random_halo_indices,range(len(random_halo_indices))):
            observer = halo_list[halo_index]
            [x,y,z] = [observer.x,observer.y,observer.z]
            observer = Observers(x,y,z,dummy,dummy,dummy,dummy,dummy,dummy)
            observer_list[observer_number] = observer
   ########################################################################################                    
        




    return observer_list
    



# This function select the halos observed by a given observer, and saves the observations
def observations(observer_number,observer_list,halo_list,mind,maxd,number_of_SNe,boxsize,number_of_cones,skyfraction,distances_from_perturbed_metric, potential_file):
    observer = observer_list[observer_number]
    [x,y,z] = [observer.x, observer.y, observer.z]

    total_mass = 0
    halo_candidates = []
    hubble = 100
    
    for i in range(len(halo_list)):
        
        halo = halo_list[i]
        #Putting the observer in origo
        [xo,yo,zo] = [halo.x - x, halo.y - y,  halo.z - z]
        [xo,yo,zo] = periodic_boundaries(xo,yo,zo,boxsize)
        
        rvec = sp.array([xo,yo,zo])
        r = la.norm(rvec)
        
        # We are using the CMB frame of reference for the velocities
        [vx,vy,vz] = [halo.vx, halo.vy, halo.vz]
#        vr = (xo*vx+yo*vy+zo*vz)/r+r*hubble
        
        # The halo is only chosen if it is within the chosen distance range
        if r < mind or r > maxd:
            continue
        
        theta = sp.arccos(zo/r)

                
        # The halo is only chosen if it is within the chosen patch
        if number_of_cones == 1:
            theta_max = sp.arccos(1-2*skyfraction)
            if theta > theta_max:
                continue
            
        if number_of_cones == 2:
            theta_max = sp.arccos(1-skyfraction)
            if theta > theta_max and sp.pi-theta > theta_max:
                continue

        halo_candidates.append(i) # Saving the index for the halo
        total_mass = total_mass+halo.mass

        
    selected_halos = mass_weighted_selection_of_halos(halo_list, halo_candidates, number_of_SNe,total_mass)

    thetas = sp.zeros(number_of_SNe) 
    phis = sp.zeros(number_of_SNe)
    rs = sp.zeros(number_of_SNe)
    vprs = sp.zeros(number_of_SNe)
    
    for count, i in enumerate(selected_halos):
        halo = halo_list[i]
        
        #Putting the observer in origo
        [xo,yo,zo] = [halo.x - x, halo.y - y,  halo.z - z]
        [xo,yo,zo] = periodic_boundaries(xo,yo,zo,boxsize)
        
        rvec = sp.array([xo,yo,zo])
        ro = la.norm(rvec)
        if distances_from_perturbed_metric == 1:
            r = calculate_distance_in_perturbed_metric(x,y,z,potential_file,xo,yo,zo,boxsize)
        else:
            r = ro        
        # We are using the CMB frame of reference for the velocities
        [vx,vy,vz] = [halo.vx, halo.vy, halo.vz]
        vpr = (xo*vx+yo*vy+zo*vz)/ro
        
        theta = sp.arccos(zo/ro)
        if xo > 0:
            phi = sp.arctan(yo/xo)
        elif xo < 0:
            phi = sp.arctan(yo/xo)+sp.pi
        elif xo == 0:
            phi = sp.sign(yo)*sp.pi/2
        else:
            print "Something must be wrong!"

        thetas[count] = theta
        phis[count] = phi
        rs[count] = r
        vprs[count] = vpr

        
    return thetas, phis,rs,vprs

def calculate_distance_in_perturbed_metric(x,y,z,potential_file,xo,yo,zo,boxsize):
    # Load potential.
    potential = sp.loadtxt(potential_file)

    # Place the observer in origo
    # This can be done by substracting the observer positions (x,y,z) from the 
    # x,y and z coordinates of the potential
    
    # No, that won't work! I would need the potential in all the new positions
    # But I only need to apply the periodicity to the d-vector
    
    # Create a set of equidistant point along the line of sight from the observer (0,0,0)
    # to the observed halo    
    res = 100
    f = sp.linspace(0,1,res)
    rvec = [xo,yo,zo]
    ro = la.norm(rvec)
    d = sp.array([f*xo+x,f*yo+y,f*zo+z])
    d[d<0] = d[d<0]+boxsize
#    pdb.set_trace()
    
#    print d    
    # The loaded potential needs to be arranged in a 3D array.
    # The indices could be given by the coordinates.
    grid = int(sp.ceil(len(potential)**(1./3)))
    pot_grid = sp.zeros((grid,grid,grid))
    row = 0
    for i in range(grid):
        for j in range(grid):
            for k in range(grid):
                pot_grid[i][j][k] = potential[row][3]
                row += 1
#    print pot_grid  
#    print "grid = ",grid              
    # map_coordinates assumes that the indices along each dimension are 0:length(dim)
    # dp can be scaled to fit this by
    d = d/boxsize # Now each coordinate is between 0 and 1
    d = d*(grid-1) # Now the coordinates are scaled to the grid dimensions
    # The -1 in (grid-1) is because the indices ranges from 0 to 7
    
    # Calculate the quantity (1-2*psi) (maybe times Delta_r?) along the line of sight
    
    import scipy.ndimage as ndi
    
    psi = ndi.map_coordinates(pot_grid,d)
#    print "d = ",d
#    print "psi = ",psi

#    x = np.arange(0,2*np.pi+np.pi/4,2*np.pi/8)
#    y = np.sin(x)
    from scipy import interpolate
#    pdb.set_trace()

    tck = interpolate.splrep(f,1-2*psi,s=0)

#    xnew = np.arange(0,2*np.pi,np.pi/50)
#    phi_new = interpolate.splev(d,tck,der=0)
    

    psi_int = interpolate.splint(0,1,tck)
#    print psi_int
    
    r = psi_int*ro
    return r   
    
    
#################################### selection of halos #############################    
    
# This function select the halos observed by a given observer
def select_halos(x,y,z,halo_list,mind,maxd,observed_halos,boxsize,number_of_cones,skyfraction, distances_from_perturbed_metric, potential_file):
    total_mass = 0
    halo_candidates = []


    
    for i in range(len(halo_list)):
        
        halo = halo_list[i]
        #Putting the observer in origo
        [xo,yo,zo] = [halo.x - x, halo.y - y,  halo.z - z]
        [xo,yo,zo] = periodic_boundaries(xo,yo,zo,boxsize)

        rvec = sp.array([xo,yo,zo])
        ro = la.norm(rvec) # This is the distance without taking metric perturbations into account
        if distances_from_perturbed_metric == 1:
            r = calculate_distance_in_perturbed_metric(x,y,z,potential_file,xo,yo,zo,boxsize)
        else:
            r = ro
        
        
        
        # The halo is only chosen if it is within the chosen distance range
        if r < mind or r > maxd:
            continue

        
        theta = sp.arccos(zo/ro)

                
        # The halo is only chosen if it is within the chosen patch
        if number_of_cones == 1:
            theta_max = sp.arccos(1-2*skyfraction)
            if theta > theta_max:
                continue
            
        if number_of_cones == 2:
            theta_max = sp.arccos(1-skyfraction)
            if theta > theta_max and sp.pi-theta > theta_max:
                continue

        halo_candidates.append(i) # Saving the index for the halo
        total_mass = total_mass+halo.mass

        
    selected_halos = mass_weighted_selection_of_halos(halo_list, halo_candidates, observed_halos,total_mass)
    
        
    return selected_halos







# This function choses the wished number of observed halos from the candidates
def mass_weighted_selection_of_halos(halo_list, halo_candidates, observed_halos, total_mass):
    cum_mass = 0
    halo_masses_cum = []

    for i in halo_candidates:
        halo = halo_list[i]
        mass = halo.mass
        cum_mass = cum_mass+mass
        halo_masses_cum.append(cum_mass/total_mass)

    # The halos are selected
    selected_halos = []        
    random.seed(0)
    for i in range(observed_halos):
        rnd = random.random()
        for j, index in zip(range(len(halo_candidates)),halo_candidates):
            if halo_masses_cum[j] >= rnd:
                selected_halos.append(index)
                break
                
    return selected_halos


######################################################################################





################################# Calculation of Hubble constants ######################

     
# This function bins and calculates Hubble constants from an observer position and a
# list of observed halos        
def Hubble(x,y,z,halo_list, selected_halos, bindistances,boxsize,distances_from_perturbed_metric,potential_file):
    hubble = 100 # Distances are in Mpc/h
    b = 1

    [rvsum, r2sum] = [0,0]
    radial_distances = []
    radial_velocities = []

    Hubbleconstants = sp.nan*sp.ones(len(bindistances))
    rvsum = sp.zeros(len(bindistances))
    r2sum = sp.zeros(len(bindistances))
    halo_counts = sp.zeros(len(bindistances))

    for i in selected_halos:
        
        halo = halo_list[i]
        #Calculate the distance to the halo
        [xo,yo,zo] = [halo.x - x, halo.y - y,  halo.z - z]
        [xo,yo,zo] = periodic_boundaries(xo,yo,zo,boxsize)
        rvec = sp.array([xo,yo,zo])
        r = la.norm(rvec)
        
        # We are using the CMB frame of reference for the velocities
        [vx,vy,vz] = [halo.vx, halo.vy, halo.vz]
        vr = (xo*vx+yo*vy+zo*vz)/r+r*hubble
                
        # Returning some values for making plots
        radial_distances.append(r)
        radial_velocities.append(vr)

        b = len(bindistances)-1
        rb = bindistances[b]
#        print "before while: b = ",b,"rb = ",rb, "and r = ",r

        while r < rb:
#            print "in while: b = ",b,"rb = ",rb, "and r = ",r
            rvsum[b] = rvsum[b]+r*vr
            r2sum[b] = r2sum[b]+r**2
            halo_counts[b] = halo_counts[b]+1
            b = b-1
            rb = bindistances[b]


    # Hubbleconstant for the innermost distance is -1, since there is no 
    #observed halos within this distance
    for b in range(len(bindistances)):
#        print "b = ",b,"and halo_counts[b] = ",halo_counts[b]
        
        if halo_counts[b] == 0:
            continue
        else:
            Hubbleconstants[b] = calculate_H(rvsum[b],r2sum[b],halo_counts[b])

    return Hubbleconstants, radial_distances, radial_velocities

def calculate_H(rvsum,r2sum,halo_count):
    
        rvmean = rvsum/halo_count
        r2mean = r2sum/halo_count
        Hloc = rvmean/r2mean
        
        return Hloc


















#################################################################################    
    
# This function applies periodic boundaries
def periodic_boundaries(x,y,z,boxsize):
        if x > boxsize/2: 
            x = x-boxsize
        if x < -boxsize/2: 
            x = x+boxsize            
   
        if y > boxsize/2: 
            y = y-boxsize
        if y < -boxsize/2: 
            y = y+boxsize       
    
        if z > boxsize/2: 
            z = z-boxsize
        if z < -boxsize/2: 
            z = z+boxsize
        
        return [x,y,z]
    
#################################################################################
    
    
    
################## Write Hubbleconstants and bindistances/number of SNe to file ##############
    
#This function write the results to file    
def print_hubbleconstants(hubblefile,bindistances,observer_list):
    f = open(hubblefile,'w')

    for bl in bindistances:
        f.write("%s\t" % bl)
        
    #The number of bins is one less than the number of bin-distances, since the 
    # minimal distance does not count
    number_of_bins = len(bindistances)
    
    for i in range(len(observer_list)):
        observer = observer_list[i]
        f.write("\n%s\t" % i)

        for b in range(1,number_of_bins):
            f.write("%s\t" % observer.Hubbleconstants[b])
            
        f.write("\n")
#################################################################################
        

                    
      
