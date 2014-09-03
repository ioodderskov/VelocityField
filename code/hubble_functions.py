from __future__ import division
import scipy as sp
import random





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
    def __init__(self,x,y,z,selected_halos,H):
        self.x = x
        self.y = y
        self.z = z
        self.selected_halos = selected_halos
        self.H = H


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
    


######################### Identification of observers #########################


#This function finds the halos that are have the characteristics specified in the parameterfile
#Or, if find_observers = 0, the positions are simple read from a file.
def find_observers(find_observers,observerfile,halo_list,masses,sub_min_m,sub_max_m,host_min_m,host_max_m):
   
    # If find_observers = 0, the observer posistion are read from a file
    if find_observers == 0:
        observer_positions = sp.loadtxt(observerfile)  

        observer_list = [None]*len(observer_positions)
    
        # Reading coordinates of the observers, and saving them in an observer list    
        for observer_number in range(len(observer_list)):
            [x,y,z] = [observer_positions[observer_number,0],observer_positions[observer_number,1],observer_positions[observer_number,2]]
            [x,y,z] = sp.array([x,y,z])/1000
            observer = Observers(x,y,z,[],[])
            observer_list[observer_number] = observer
        
        



    
    # The observer positions are read or found:
    if find_observers == 1:

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
            observer = Observers(x,y,z,[],[])
            observer_list[observer_number] = observer
        





    return observer_list
    
    

   
    
    
    
    

# This function calculates the Hubble constants for a given observer
def find_hubble_constants_for_observer(x,y,z,halo_list,mind,maxd,observed_halos,bindistances,boxsize,number_of_cones,sky_cover):
    selected_halos = select_halos(x,y,z,halo_list,mind,maxd,observed_halos,boxsize,number_of_cones,sky_cover)
    Hubbleconstants, radial_distances, radial_velocities = Hubble(x,y,z,halo_list, selected_halos, bindistances, boxsize)
    return Hubbleconstants, radial_distances, radial_velocities, selected_halos 
    
    
    
    
#################################### selection of halos #############################    
    
# This function select the halos observed by a given observer
def select_halos(x,y,z,halo_list,mind,maxd,observed_halos,boxsize,number_of_cones,skyfraction):
    total_mass = 0
    halo_candidates = []


    
    for i in range(len(halo_list)):
        
        halo = halo_list[i]
        #Putting the observer in origo
        [xo,yo,zo] = [halo.x - x, halo.y - y,  halo.z - z]
        [xo,yo,zo] = periodic_boundaries(xo,yo,zo,boxsize)
        
            
        r = sp.sqrt(xo*xo+yo*yo+zo*zo)
        
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
def Hubble(x,y,z,halo_list, selected_halos, bindistances,boxsize):
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
        r = sp.sqrt(xo*xo+yo*yo+zo*zo)
        
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
    
    
    
################## Write Hubbleconstants and bindistances to file ##############
    
#This function write the results to file    
def print_hubbleconstants(hubblefile,bindistances,observer_list):
    f = open(hubblefile,'w')

    for bl in bindistances:
        f.write("%s\t" % bl)
        
    #The number of bins is one less than the number of bin-distances, since the 
    # minimal distance does not count
    number_of_bins = len(bindistances)
    
    for i, observer in enumerate(observer_list):
        f.write("\n%s\t" % i)

        for b in range(1,number_of_bins):
            f.write("%s\t" % observer.Hubbleconstants[b])
            
        f.write("\n")
            
      