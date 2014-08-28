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
    
    
    
    
    
    
    
    
    
    

# This function calculates the Hubble constants for a given observer
def find_hubble_constants_for_observer(halo_index,x,y,z,halo_list,mind,maxd,observed_halos,bindistances,boxsize,number_of_cones,sky_cover):
    selected_halos = select_halos(halo_index,x,y,z,halo_list,mind,maxd,observed_halos,boxsize,number_of_cones,sky_cover)
    selected_and_sorted_halos = sort_halos(x,y,z,halo_list, selected_halos,boxsize)
    Hubbleconstants, radial_distances, radial_velocities, xo, yo, zo = Hubble(x,y,z,halo_list, selected_and_sorted_halos, bindistances, boxsize)
    return Hubbleconstants, xo, yo, zo, radial_distances, radial_velocities, selected_halos 
    
    
    
    
#################################### selection of halos #############################    
    
# This function select the halos observed by a given observer
def select_halos(halo_index,x,y,z,halo_list,mind,maxd,observed_halos,boxsize,number_of_cones,skyfraction):
    total_mass = 0
    halo_candidates = []


    
    for i in range(len(halo_list)):
        # Checking if the halo is the observer. In that case, skip.
        if i == halo_index:
            continue
        
        halo = halo_list[i]
        #Putting the observer in origo
        [xo,yo,zo] = [halo.x - x, halo.y - y,  halo.z - z]
        [xo,yo,zo] = periodic_boundaries(xo,yo,zo,boxsize)
        
            
        r = sp.sqrt(xo**2+yo**2+zo**2)
        
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
    print "Observed_halos = ", observed_halos

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






##################################### Sorting of halos #############################


# This function sorts a list of halos according to their distance
def sort_halos(x,y,z,halo_list, selected_halos, boxsize):
    distances = []
    for i in selected_halos:
        halo = halo_list[i]
        [xo,yo,zo] = [halo.x - x, halo.y - y,  halo.z - z]
        [xo,yo,zo] = periodic_boundaries(xo,yo,zo,boxsize)
        r = sp.sqrt(xo**2+yo**2+zo**2)
        distances.append(r)
        
    sort_index = sp.argsort(distances)
    selected_and_sorted_halos = sp.array(selected_halos)[sort_index]
    
    return selected_and_sorted_halos
    

##################################################################################













################################# Calculation of Hubble constants ######################

     
# This function bins and calculates Hubble constants from an observer position and a
# list of observed halos        
def Hubble(x,y,z,halo_list, selected_and_sorted_halos, bindistances,boxsize):
    hubble = 100 # Distances are in Mpc/h
    b = 1
    maxd_bin = bindistances[b]
    [rvsum, r2sum] = [0,0]
    Hubbleconstants = []
    radial_distances = []
    radial_velocities = []
    xos = []
    yos = []
    zos = []
    count = 0
    for i in selected_and_sorted_halos:
        halo = halo_list[i]
        [xo,yo,zo] = [halo.x - x, halo.y - y,  halo.z - z]
        [xo,yo,zo] = periodic_boundaries(xo,yo,zo,boxsize)
        r = sp.sqrt(xo**2+yo**2+zo**2)
        
        # We are using the CMB frame of reference for the velocities
        vx = halo.vx
        vy = halo.vy
        vz = halo.vz
        
        vr = (xo*vx+yo*vy+zo*vz)/r+r*hubble
        
        # Returning some values for making plots
        radial_distances.append(r)
        radial_velocities.append(vr)
        xos.append(xo)
        yos.append(yo)
        zos.append(zo)
#        vrtot = vr+r*hubble


        
        if maxd_bin < r:
            # Calculating the Hubble constant
            rvmean = rvsum/count
            r2mean = r2sum/count
            Hloc = rvmean/r2mean 
            Hubbleconstants.append(Hloc)
            # Move to next bin
            b = b+1
            maxd_bin = bindistances[b]
            
        rvsum = rvsum + r*vr
        r2sum = r2sum + r**2
        count = count + 1
        
    return Hubbleconstants, radial_distances, radial_velocities, xos, yos, zos

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
    

