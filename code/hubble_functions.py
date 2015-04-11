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
import tarfile
import healpy as hp
import hubble_classes as hc
import gravitational_instability as gi
import copy


plot_field = 1



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


def load_CoDECS_catalogue(halocatalogue_file):
    tar = tarfile.open(halocatalogue_file,mode='r')

    for i,catalogue in enumerate(tar.getmembers()):
        f=tar.extractfile(catalogue)
        content = sp.loadtxt(f)
        if not 'halocatalogue' in locals():
            halocatalogue = content
        else:
            sp.vstack((halocatalogue,content))
    
    tar.close()
    
    return halocatalogue
    
    
def initiate_observers_CoDECSsubhalos(parameters,halos):
    groups = load_CoDECS_catalogue(parameters.CoDECShosts_file)

    massunit = 1e10 # Msun/h

    groupmasses = sp.array(groups[:,2])*massunit
    group_IDs = groups[:,0]    
    
    submasses = sp.array([halo.mass for halo in halos])
    ID_hosts = sp.array([halo.ID_host for halo in halos])

    print "Total group masses = ", sp.sum(groupmasses)
    print "Total sub masses = ", sp.sum(submasses)
    
    localgroup_indices = sp.array(range(len(halos)))[(parameters.sub_min_m < submasses) & (submasses < parameters.sub_max_m)]    
    virgo_indices = (parameters.host_min_m < groupmasses) & (groupmasses < parameters.host_max_m)

    virgo_IDs = group_IDs[virgo_indices]
    
    observer_indices = [localgroup_index for localgroup_index in localgroup_indices if ID_hosts[localgroup_index] in virgo_IDs]


    observers = [None]*len(observer_indices)
    for ob_number, ob_index in enumerate(observer_indices):
        halo = halos[ob_index]
        position = halo.position
        observers[ob_number] = hc.Observer(ob_number,position)
        
#    pdb.set_trace()
    return observers[0:parameters.number_of_observers]


def load_grid(gridfile):
    # It is not really a halocatalogue...
    halocatalogue = sp.loadtxt(gridfile)
    return halocatalogue   

    
    
# Load the halo catalogue
def load_halocatalogue(parameters,halocatalogue_file):
    
    if parameters.CoDECS:
        halocatalogue_unsorted = load_CoDECS_catalogue(halocatalogue_file)
        halocatalogue = sp.array(sorted(halocatalogue_unsorted,
                                        key=lambda halocatalogue_unsorted: halocatalogue_unsorted[3]))

    elif parameters.use_grid:
        halocatalogue_unsorted = load_grid(parameters.gridfile)
        halocatalogue = sp.array(sorted(halocatalogue_unsorted, 
                                        key=lambda halocatalogue_unsorted: halocatalogue_unsorted[3]))            

    else:
        halocatalogue_unsorted = sp.loadtxt(halocatalogue_file)
        halocatalogue = sp.array(sorted(halocatalogue_unsorted, 
                                        key=lambda halocatalogue_unsorted: halocatalogue_unsorted[2]))
                                    
    return halocatalogue

    
    
def initiate_halos(parameters, halocatalogue):
    n_halos = len(halocatalogue)
    halos = [None]*n_halos

    if parameters.CoDECS:
        
        massunit = 1e10 # Msun/h 
        for h in range(n_halos):
            position = halocatalogue[h,[9,10,11]]/1000.
            velocity = halocatalogue[h,[12,13,14]]
            mass = halocatalogue[h,3]*massunit
            ID = int(halocatalogue[h,1])
            ID_host = int(halocatalogue[h,0])
            halo = hc.Halo(position,velocity,mass,ID,ID_host,h)
            halos[h] = halo

#        print "position of last halo = ", position
#        print "changing mass of last halo from", mass, "to", mass*1e16
#        halos[h] = hc.Halo(position,velocity,mass*1e16,ID,ID_host,h)
            
#        position = sp.array([0,0,0])
#        velocity = sp.array([0,0,0])
#        mass = sp.sum([halocatalogue[h,3]*massunit for h in range(n_halos)])*1e16
#        halos[-1] = hc.Halo(position,velocity,mass,0,0,n_halos)
 
    elif parameters.use_grid:

        for h in range(n_halos):
            position = halocatalogue[h,[0,1,2]]
            velocity = halocatalogue[h,[4,5,6]]
            mass = halocatalogue[h,3]
            ID = h
            ID_host = -1
            
            halo = hc.Halo(position,velocity,mass,ID,ID_host,h)
            halos[h] = halo
            

    else:    
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
            
            if parameters.CoDECS:
                observers = initiate_observers_CoDECSsubhalos(parameters,halos)
                
            else:
                observers = initiate_observers_subhalos(parameters,halos)
    
        if parameters.observer_choice == 'random_halos':
            observers = initiate_observers_random_halos(parameters,halos)
    
        if parameters.observer_choice == 'random_positions':
            observers = initiate_observers_random_positions(parameters)
            
        if parameters.observer_choice == 'indices_from_file':
            observers = initiate_observers_indices_from_file(parameters,halos)
    
    
    return observers
    
    
def initiate_observers_from_file(parameters):
    observer_positions = sp.array(sp.loadtxt(parameters.observerfile))
    observers = [None]*len([observer_positions])
    
    for ob in range(len([observer_positions])):
        
        position = observer_positions[ob,[0,1,2]]/1000
        observers[ob] = hc.Observer(ob,position)
        
    return observers[0:parameters.number_of_observers]
    
    
def initiate_observers_random_halos(parameters,halos):
    random_indices = sp.random.random_integers(0,len(halos)-1,parameters.number_of_observers)
    observers = [None]*len(random_indices)
    
    for ob_index, ob_number in zip(random_indices,range(len(random_indices))):
        halo = halos[ob_index]
        position = halo.position
        
        observers[ob_number] = hc.Observer(ob_number,position)
        
    return observers[0:parameters.number_of_observers]
    
def initiate_observers_subhalos(parameters,halos):
    print "WARNING! This function needs to be tested before use. Se evt. den tilsvarende fkt for CoDECS."
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
        observers[ob_number] = hc.Observer(ob_number,position)
    
    return observers[0:parameters.number_of_observers]

def initiate_observers_random_positions(parameters):

    sp.random.seed(0)
    observer_positions = parameters.boxsize*sp.rand(parameters.number_of_observers,3)
#    print "shape(observer_positions)", observer_positions
    
    observers = [None]*len(observer_positions)
        
    for ob in range(len(observer_positions)):
        position = observer_positions[ob]
        observers[ob] = hc.Observer(ob,position)
    
    
    return observers[0:parameters.number_of_observers]
    
def initiate_observers_indices_from_file(parameters,halos):

    observer_indices = sp.array(sp.loadtxt(parameters.observer_indices_file))
    observers = [None]*len([observer_indices])
    
    for ob_index, ob_number in zip([observer_indices],range(len([observer_indices]))):
        ob_index = int(ob_index)
        halo = halos[ob_index]
        position = halo.position
        
        observers[ob_number] = (hc.Observer(ob_number,position))
        
    return observers[0:parameters.number_of_observers]

def periodic_coordinate(parameters,coordinate):
    
    if coordinate > parameters.boxsize/2:
        coordinate = coordinate-parameters.boxsize
    if coordinate < -parameters.boxsize/2:
        coordinate = coordinate+parameters.boxsize
        
    return coordinate
    
    
def periodic_boundaries(parameters,position_obs,position_halo):
    

    position = position_halo-position_obs    
    x,y,z = position[0],position[1],position[2]
      
    x = periodic_coordinate(parameters,x)+position_obs[0]
    y = periodic_coordinate(parameters,y)+position_obs[1]    
    z = periodic_coordinate(parameters,z)+position_obs[2]
    
    return sp.array([x,y,z])
    
    
    
def spherical_coordinates(parameters,position_obs,position_op):


    rvec = position_op-position_obs
    r = la.norm(rvec)

    # Just to prevent problems if the observer is on top of the halo
    if r == 0:
        r = 1e-15         
        
    theta = sp.arccos((position_op[2]-position_obs[2])/r)
    phi = sp.arctan2(position_op[1]-position_obs[1],position_op[0]-position_obs[0])+sp.pi

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
    
def center_of_mass(parameters,observer_position,candidates):
    masses = sp.array([halo.mass for halo in candidates])    


    positions = sp.array([halo.position for halo in candidates])

    
    positions_op = sp.empty_like(positions)    

    for i, position in enumerate(positions):

        position_op = periodic_boundaries(parameters,observer_position,position)
        positions_op[i] = position_op

    velocities = sp.array([halo.velocity for halo in candidates])

    
    position_CoM = sp.average(positions_op,axis=0,weights=masses)


    velocity_CoM = sp.average(velocities,axis=0,weights=masses)

    
    total_mass = sp.sum(masses)

    return position_CoM, velocity_CoM, total_mass
    
    
    

def determine_CoM_for_these_halos(parameters,survey_positions,survey_masses,observer_position,local_halos,candidates):

    print "len(candidates) = ", len(candidates)

    if len(candidates) == 0:
        print "No observed halo to calculate CoM for"
        return 0,0,0,0,0,0,0,0,0

    position_CoM, \
    velocity_CoM, total_mass = center_of_mass(parameters,observer_position,candidates)        
    
    

    if parameters.correct_for_peculiar_velocities:

        velocity_correction = gi.velocity_from_matterdistribution(parameters,observer_position,position_CoM,survey_positions,survey_masses)

        velocity_CoM_nc = copy.copy(velocity_CoM)
        velocity_CoM = velocity_CoM - velocity_correction
             
    r_CoM, theta_CoM, phi_CoM = spherical_coordinates(parameters,observer_position,
                                                position_CoM)

    vr_peculiar_CoM = sp.dot(position_CoM-observer_position,velocity_CoM)/r_CoM
    
    vr_CoM = vr_peculiar_CoM + r_CoM*100
 
    if plot_field:
        import matplotlib.pyplot as plt
        plt.figure()
        ax = plt.gca()
        ax.set_xlim([observer_position[0]-parameters.survey_radius,observer_position[0]+parameters.survey_radius])
        ax.set_ylim([observer_position[1]-parameters.survey_radius,observer_position[1]+parameters.survey_radius])
#        halo_positions = sp.array(survey_positions)
#        halo_indices = (observer_position[2]-50 < halo_positions[:,2]) & (halo_positions[:,2] < observer_position[2]+50) 
#        ax.plot(halo_positions[halo_indices,0],halo_positions[halo_indices,1],'r*')
        ax.plot(observer_position[0],observer_position[1],'g*',markersize=20)
        candidate_positions = sp.array([candidate.position for candidate in candidates])
        candidate_positions_op = sp.empty_like(candidate_positions)
        for i,candidate_position in enumerate(candidate_positions):
            candidate_positions_op[i] = periodic_boundaries(parameters,observer_position,candidate_position)
            
        ax.plot(candidate_positions_op[:,0],candidate_positions_op[:,1],'y*')
        
        ax.plot(position_CoM[0],position_CoM[1],'bx')
        
        scale = 100
        candidate_velocities = sp.array([candidate.velocity for candidate in candidates])
        ax.quiver(candidate_positions_op[:,0],candidate_positions_op[:,1],candidate_velocities[:,0],candidate_velocities[:,1],color='y',scale_units='inches',scale=scale)
        ax.quiver(position_CoM[0],position_CoM[1],velocity_CoM_nc[0],velocity_CoM_nc[1],color='r',scale_units='inches',scale=scale)
        ax.quiver(position_CoM[0],position_CoM[1],velocity_correction[0],velocity_correction[1],color='black',scale_units='inches',scale=scale)
        print "Velocity = ", velocity_CoM_nc
        print "Velocity correction  = ", velocity_correction

        print "total_mass = ", total_mass
        plt.show()


    return  position_CoM, \
            r_CoM, theta_CoM, phi_CoM,\
            vr_peculiar_CoM, vr_CoM, \
            total_mass   

def plot_velocities(pos,vel,ax,color,scale):

    xs = pos[:,0]
    ys = pos[:,1]
    vxs = vel[:,0]
    vys = vel[:,1]

    ax.quiver(xs,ys,vxs,vys,color=color,scale_units='inches',scale=scale)


    

def select_candidates(parameters,candidates):

    if parameters.observed_halos == 'all':
        selected_candidates = range(len(candidates))    

    if parameters.observed_halos == 'mass_weighted':
        selected_candidates = mass_weighted_selection_of_halos(parameters,candidates)

    if parameters.observed_halos == 'random':
        selected_candidates = random_selection_of_halos(parameters,candidates)
        
    if parameters.observed_halos == 'specified_mass':
        selected_candidates = specified_mass_selection_of_halos(parameters,candidates)
        
    return selected_candidates

   
def specified_mass_selection_of_halos(parameters,candidates):
    
    selected_candidates = []
    masses = sp.array([candidate.mass for candidate in candidates])
    indices_of_candidates_of_specified_mass = sp.array(range(len(candidates)))[(parameters.SN_mass_min < masses) & (masses < parameters.SN_mass_max)]
#    pdb.set_trace()
    random.seed(0)
    if len(indices_of_candidates_of_specified_mass) == 0:
        print "No halos to observe"
    else:
        while len(selected_candidates) < parameters.number_of_SNe:
            rnd_number = random.randint(0,len(indices_of_candidates_of_specified_mass)-1)
            candidate_number = indices_of_candidates_of_specified_mass[rnd_number]
            selected_candidates.append(candidate_number)
#            print candidates[candidate_number].mass
            
    return selected_candidates
            
        
    
def random_selection_of_halos(parameters,candidates):
    
    selected_candidates = []
    random.seed(0)
    for n in range(parameters.number_of_SNe):
        rnd_candidate_number = random.randint(0,len(candidates)-1) # For some reason, the endpoint is included in the interval.
#        rnd_candidate_number = candidates[rnd_int]
        selected_candidates.append(rnd_candidate_number)
        
 
    return selected_candidates

def mass_weighted_selection_of_halos(parameters,candidates):
    
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
    
    if parameters.test_isotropy:
        
        print "I am not printing these Hubbleconstants to a file.\nInstead, I will save the observers and their observations in an array"
        sp.save(parameters.path+'observers_isotropy',observers)
#        sp.save(parameters.path+'parameters',parameters)
#        f = open(parameters.hubblefile,'w')
#        f.write("085\t")
#        
#        
#        for theta, phi in zip(parameters.directions[0],parameters.directions[1]):
#            f.write("(%s,%s)\t" % (theta,phi))
#            
#        for ob_number, ob in enumerate(observers):
#            f.write("\n%s\t" % ob_number)
#            
#            for direction in range(parameters.number_of_directions):
#                f.write("%s\t" % ob.Hubbleconstants[direction])
#                
#        f.close()
            

    elif parameters.vary_skyfraction:
        
        for row, skyfraction in enumerate(parameters.skyfractions):
            
            f = open(parameters.hubblefile+'_'+parameters.observer_choice+str(skyfraction)+'.txt','w')

            for bl in parameters.bindistances:
                f.write("%s\t" % bl)
                
            for ob_number, ob in enumerate(observers):
               f.write("\n%s\t" % ob_number)
               
               for b in range(1,len(parameters.bindistances)):
                   f.write("%s\t" % ob.Hubbleconstants[row][b])
                   
               f.write("\n")
               
            f.close()
            
    elif parameters.vary_number_of_SNe:

        f = open(parameters.hubblefile+'_'+parameters.observer_choice+'.txt','w')
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
        sp.save(parameters.path+'observers_coma',observers)
    
        f = open(parameters.hubblefile+'_'+parameters.observer_choice+'.txt','w')
        
        for bl in parameters.bindistances:
            f.write("%s\t" % bl)
            
        for ob_number, ob in enumerate(observers):
           f.write("\n%s\t" % ob_number)
           
           for b in range(1,len(parameters.bindistances)):
               f.write("%s\t" % ob.Hubbleconstants[b])
           
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
        
        if r < parameters.bindistances[0]:
            print "The CoM distances is less than the minimal bindistance"
            return sp.zeros_like(parameters.bindistances)
        
        while r < rb:
            
            rvsum[b] = rvsum[b]+r*vr
            r2sum[b] = r2sum[b]+r**2
            halo_counts[b] = halo_counts[b]+1
            b = b-1
#            pdb.set_trace()
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
    
    
def calculate_Hs_for_varying_directions(parameters,observed_halos):
    
    direction_Hubbleconstants_array = sp.zeros((parameters.number_of_directions,len(parameters.bindistances)))
    cones = []
    for direction_row, theta, phi in zip(range(parameters.number_of_directions),parameters.directions[0],parameters.directions[1]):
        halos_in_cone = [oh for oh in observed_halos if hp.rotator.angdist([theta,phi],[oh.theta,oh.phi])< parameters.max_angular_distance]
        selected_halos_in_cone_numbers = select_candidates(parameters,halos_in_cone)
        selected_halos_in_cone = [halos_in_cone[i] for i in selected_halos_in_cone_numbers]
        cone = hc.Cone(theta,phi,selected_halos_in_cone)
        cones.append(cone)

        direction_Hubbleconstants_array[direction_row] = calculate_Hs_for_these_observed_halos(parameters,selected_halos_in_cone)
     
        
    return direction_Hubbleconstants_array, cones
        
def calculate_Hs_for_varying_number_of_SNe(parameters,observed_halos):
    
    number_of_SNe_Hubbleconstants_array = sp.zeros((len(parameters.numbers_of_SNe),len(parameters.bindistances)))

    
    for row, number_of_SNe in enumerate(parameters.numbers_of_SNe):
        observed_sample_of_halos = random.sample(observed_halos,number_of_SNe)
        
        number_of_SNe_Hubbleconstants_array[row] = calculate_Hs_for_these_observed_halos(parameters,observed_sample_of_halos)
        
    return number_of_SNe_Hubbleconstants_array

def load_CoDECS_catalogue(halocatalogue_file):
    tar = tarfile.open(halocatalogue_file,mode='r')

    for i,catalogue in enumerate(tar.getmembers()):
        f=tar.extractfile(catalogue)
        content = sp.loadtxt(f)
        if not 'halocatalogue' in locals():
            halocatalogue = content
        else:
            sp.vstack((halocatalogue,content))
    
    tar.close()
    
    return halocatalogue
 


def read_snapshot(parameters):

    if parameters.CoDECS:
        tar = tarfile.open(parameters.snapshot_file,mode='r')
        for i,snapshot in enumerate(tar.getmembers()):
            f=tar.extractfile(snapshot)
    else:
        f = open(parameters.snapshot_file,'r')
    
    size_unit = 4
    print "reading file:", parameters.snapshot_file
    size_header = 5*size_unit
    f.seek(size_header)

    Npart = sp.fromfile(f,sp.uint32,count=6)
    Massarr = sp.fromfile(f,sp.double,count=6)
    Time = sp.fromfile(f,sp.double,count=1)
    Redshift = sp.fromfile(f,sp.double,count=1)
    FlagSfr = sp.fromfile(f,sp.int32,count=1)
    FlagFeedback = sp.fromfile(f,sp.int32,count=1)
    Nall = sp.fromfile(f,sp.int32,count=6)
    Flagcooling = sp.fromfile(f,sp.int32,count=1)
    NumFiles = sp.fromfile(f,sp.int32,count=1)
    BoxSize = sp.fromfile(f,sp.double,count=1)
    Omega0 = sp.fromfile(f,sp.double,count=1)
    OmegaLambda = sp.fromfile(f,sp.double,count=1)
    HubbleParam = sp.fromfile(f,sp.double,count=1)
    FlagAge = sp.fromfile(f,sp.int32,count=1)
    FlagMetals = sp.fromfile(f,sp.int32,count=1)
    NallHW = sp.fromfile(f,sp.int32,count=6)
    flag_entr_ics = sp.fromfile(f,sp.int32,count=1)

    size_fortranstuff = 6*size_unit
    f.seek(size_header+256+size_fortranstuff)

    print "Npart = ", Npart
    print "Massarr = ", Massarr

    halos = [None]*sp.sum(Npart)
    for particle_type in [0,1]:
        N = int(Npart[particle_type])
        massunit = 1e10 #Msun/h
        M = sp.double(Massarr[particle_type])*massunit
    #    print M
        
        pos = sp.fromfile(f,sp.single,count=N*3)
    
        size_more_fortran_stuff = 6*size_unit
        f.seek(size_more_fortran_stuff,1)
    
        vel = sp.fromfile(f,sp.single,count=N*3)
    
    
        xindices = sp.mod(range(N*3),3) == 0
        yindices = sp.mod(range(N*3),3) == 1
        zindices = sp.mod(range(N*3),3) == 2
        
        xs = pos[xindices]/1000.
        ys = pos[yindices]/1000.
        zs = pos[zindices]/1000.
    
        vxs = vel[xindices]    
        vys = vel[yindices]
        vzs = vel[zindices]    
    
        for p in range(N):
            position = sp.array([xs[p],ys[p],zs[p]])
            velocity = sp.array([vxs[p],vys[p],vzs[p]])
            ID = p
            ID_host = -1
            
            halo = hc.Halo(position,velocity,M,ID,ID_host,p)
            halos[p] = halo
    
#    position = sp.array([0,0,0])
#    velocity = sp.array([0,0,0])
#    M = N*M*1e16    
#    halos[p] = hc.Halo(position,velocity,M,ID,ID_host,p)
        
    if parameters.CoDECS:
        tar.close()

	f.close()

    return halos

    
    
   
                    
      
