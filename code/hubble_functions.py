from __future__ import division
import scipy as sp
import random
import scipy.linalg as linalg
import scipy.ndimage as ndi
import hubble_classes as hc
from scipy import interpolate
import pdb
import tarfile
import healpy as hp
import gravitational_instability as gi
import copy
import gc
import matplotlib.pyplot as plt
import sys


plot_field = 0



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



    
def initiate_observers_CoDECSsubhalos(parameters):
    groups = load_CoDECS_catalogue(parameters.CoDECShosts_file)

    massunit = 1e10 # Msun/h

    groupmasses = sp.array(groups[:,2])*massunit
    group_IDs = groups[:,0]    
    
    submasses = sp.array([halo.mass for halo in parameters.halos])
    ID_hosts = sp.array([halo.ID_host for halo in parameters.halos])
    
    localgroup_indices = sp.array(range(len(parameters.halos)))[(parameters.sub_min_m < submasses) & (submasses < parameters.sub_max_m)]    
    virgo_indices = (parameters.host_min_m < groupmasses) & (groupmasses < parameters.host_max_m)

    virgo_IDs = group_IDs[virgo_indices]
    
    observer_indices = [localgroup_index for localgroup_index in localgroup_indices if ID_hosts[localgroup_index] in virgo_IDs]


    observers = [None]*len(observer_indices)
    for ob_number, ob_index in enumerate(observer_indices):
        halo = parameters.halos[ob_index]
        position = halo.position
        velocity = halo.velocity
        mass = halo.mass
        observers[ob_number] = hc.Observer(ob_index,position,velocity,mass)
        
    return observers[0:parameters.number_of_observers]


def load_grid(gridfile):
    # It is not really a halocatalogue...
    halocatalogue = sp.loadtxt(gridfile)
    return halocatalogue   

    
    
# Load the halo catalogue
def load_halocatalogue(parameters,halocatalogue_file):
    
    if parameters.CoDECS:
        print "Loading CoDECS halocatalogue, file:", halocatalogue_file
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


def load_CoDECS_catalogue(halocatalogue_file):
    tar = tarfile.open(halocatalogue_file,mode='r')

    for i,catalogue in enumerate(tar.getmembers()):
        f=tar.extractfile(catalogue)
        content = sp.loadtxt(f)
        if not 'halocatalogue' in locals():
            halocatalogue = content
#            pdb.set_trace()
        else:
            halocatalogue = sp.vstack((halocatalogue,content))
#            pdb.set_trace()

#    pdb.set_trace()    
    tar.close()

    print "len(halocatalogue) = ", len(halocatalogue)
    
    return halocatalogue    
    
    
def initiate_halos(parameters, halocatalogue):
    n_halos = len(halocatalogue)
    halos = [None]*n_halos

    if parameters.CoDECS:

        halos = []
        
        massunit = 1e10 # Msun/h 
        for h in range(n_halos):
            position = halocatalogue[h,[9,10,11]]/1000.
            if parameters.calculate_pairwise_velocities:
                mass = halocatalogue[h,4]*massunit
                if (mass < parameters.min_halo_mass) | (parameters.max_halo_mass < mass):
                    continue
            else:
                mass = halocatalogue[h,3]*massunit
            velocity = halocatalogue[h,[12,13,14]]
            ID = int(halocatalogue[h,1])
            ID_host = int(halocatalogue[h,0])
            halo = hc.Halo(position,velocity,mass,ID,ID_host)
            #halos[h] = halo
            halos.append(halo)

 
    elif parameters.use_grid:

        for h in range(n_halos):
            position = halocatalogue[h,[0,1,2]]
            velocity = halocatalogue[h,[4,5,6]]
            mass = halocatalogue[h,3]
            ID = h
            ID_host = -1
            
            halo = hc.Halo(position,velocity,mass,ID,ID_host)
            halos[h] = halo
            

    else: 
        halos = []
        subhalos = []
        for h in range(n_halos):
            position = halocatalogue[h,[8,9,10]]
            velocity = halocatalogue[h,[11,12,13]]
            if parameters.calculate_pairwise_velocities:
                mass = halocatalogue[h,20]
                if (mass < parameters.min_halo_mass) | (parameters.max_halo_mass < mass):
                    continue
            else:
                mass = halocatalogue[h,2]
            ID = int(halocatalogue[h,0])
            ID_host = int(halocatalogue[h,33])
            halo = hc.Halo(position,velocity,mass,ID,ID_host)
            if parameters.calculate_pairwise_velocities:
                halos.append(halo)
            else:
                if ID_host == -1:
                    halos.append(halo)
                else:
			        subhalos.append(halo)
        
        parameters.subhalos = sp.array(subhalos)
    
    parameters.halos = sp.array(halos)    
    masses = sp.array([halo.mass for halo in parameters.halos])
    print "The number of saved halos is", len(parameters.halos)
    print "min_halo_mass, max_halo_mass = ", parameters.min_halo_mass, parameters.max_halo_mass
    print "min_mass, max_mass in catalogue = ", sp.amin(masses), sp.amax(masses)  

#    sys.exit("Now you know about the halos!")
    return 1
    

def initiate_observers(parameters):
    
    if parameters.use_lightcone:
        observers = initiate_observers_from_file(parameters)
        
    else:
    
        if parameters.observer_choice == 'from_file':
            observers = initiate_observers_from_file(parameters)
    
        if parameters.observer_choice == 'subhalos':
            
            if parameters.CoDECS:
                observers = initiate_observers_CoDECSsubhalos(parameters)
                
            else:
                observers = initiate_observers_subhalos(parameters)
    
        if parameters.observer_choice == 'random_halos':
            observers = initiate_observers_random_halos(parameters)
    
        if parameters.observer_choice == 'random_positions':
            observers = initiate_observers_random_positions(parameters)
            
        if parameters.observer_choice == 'indices_from_file':
            observers = initiate_observers_indices_from_file(parameters)
        if parameters.observer_choice == 'all':
            observers = initiate_observers_all(parameters)
    
    print "from initiate_observer, observers = ", observers 
    return observers
    
def initiate_observers_all(parameters):

    observers = sp.empty(len(parameters.halos),dtype=object)
    for observer_number,halo in enumerate(parameters.halos):
        position = halo.position
        velocity = halo.velocity
        mass = halo.mass
        observers[observer_number] = hc.Observer(observer_number,position,velocity,mass)

#    sys.exit("Have initialized observers. Exiting.")    
    return observers
    
def initiate_observers_from_file(parameters):
    observer_positions = sp.array(sp.loadtxt(parameters.observerfile))
    observers = [None]*len([observer_positions])
    
    for ob in range(len([observer_positions])):
        
        position = observer_positions[ob,[0,1,2]]/1000
        velocity = sp.array([0,0,0])
        mass = 0
        observers[ob] = hc.Observer(ob,position,velocity,0)
        
    return observers[0:parameters.number_of_observers]
    
    
def initiate_observers_random_halos(parameters):
    sp.random.seed(0)
    random_indices = sp.random.random_integers(0,len(parameters.halos)-1,parameters.number_of_observers)
    
    observers = [None]*len(random_indices)
    
    for ob_index, ob_number in zip(random_indices,range(len(random_indices))):
        halo = parameters.halos[ob_index]
        position = halo.position
        velocity = halo.velocity
        mass = halo.mass
        
        observers[ob_number] = hc.Observer(ob_index,position,velocity,mass)
        
    return observers[0:parameters.number_of_observers]
    
def initiate_observers_subhalos(parameters):

    localgroup_halos = sp.array([halo for halo in parameters.subhalos\
            if (parameters.sub_min_m < halo.mass) & (halo.mass < parameters.sub_max_m)])
    localgroup_halo_numbers = sp.array([halo_number for halo,halo_number in zip(localgroup_halos,range(len(parameters.subhalos)))\
            if (parameters.sub_min_m < halo.mass) & (halo.mass < parameters.sub_max_m)])

    virgo_halos = sp.array([halo for halo in parameters.halos\
            if (parameters.host_min_m < halo.mass) & (halo.mass < parameters.host_max_m)])

    virgo_IDs = sp.array([halo.ID for halo in virgo_halos])
    observer_halos = sp.array([halo for halo in localgroup_halos\
            if halo.ID_host in virgo_IDs])

    observer_halo_numbers = sp.array([halo_number for halo,halo_number in zip(localgroup_halos,localgroup_halo_numbers)\
            if halo.ID_host in virgo_IDs])

    observers = sp.empty(len(observer_halos),dtype=object)    
    for observer_number, observer_halo, observer_halo_number in zip(range(len(observer_halos)),observer_halos,observer_halo_numbers):
        position = observer_halo.position
        velocity = observer_halo.velocity
        mass = observer_halo.mass
        observers[observer_number] = hc.Observer(observer_halo_number,position,velocity,mass)

    submasses = sp.array([halo.mass for halo in parameters.subhalos])
    hostmasses = sp.array([halo.mass for halo in parameters.halos])
    print "The number of halos is", len(parameters.halos)
    print "The number of subhalos is", len(parameters.subhalos)

    print("The mass range of the subhalos is %0.2f --> %0.2f" % (sp.amin(submasses), sp.amax(submasses)))
    print("The mass range of the hosthalos is %0.2f --> %0.2f" % (sp.amin(hostmasses), sp.amax(hostmasses)))
    print "The number of localgroup_halos is", len(localgroup_halos)
    print "The number of virgo_halos is", len(virgo_halos)
    print "The number of potential observers is", len(observers)
    print "from hubble_functions, observers = ", observers
    
    return observers[0:parameters.number_of_observers]
    


def initiate_observers_random_positions(parameters):

    sp.random.seed(0)
    observer_positions = parameters.boxsize*sp.rand(parameters.number_of_observers,3)
    
    observers = [None]*len(observer_positions)
        
    for ob in range(len(observer_positions)):
        position = observer_positions[ob]
        position = sp.array([0,0,])
        mass = 0
        observers[ob] = hc.Observer(ob,position,velocity,mass)
    
    
    return observers[0:parameters.number_of_observers]
    
def initiate_observers_indices_from_file(parameters):

    observer_indices = sp.array(sp.loadtxt(parameters.observer_indices_file))
    observers = sp.empty(len(observer_indices),dtype=object)
        
    for ob_number,ob_index in enumerate(observer_indices):
        
        ob_index = int(ob_index)
        halo = parameters.halos[ob_index]
        position = halo.position
        velocity = halo.velocity
        mass = halo.mass
        observers[ob_number] = hc.Observer(ob_number,position,velocity,mass)
        
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
    r = linalg.norm(rvec)

    # Just to prevent problems if the observer is on top of the halo
    if r == 0:
        r = 1e-15         
        
    theta = sp.arccos((position_op[2]-position_obs[2])/r)
    phi = sp.arctan2(position_op[1]-position_obs[1],position_op[0]-position_obs[0])+sp.pi

    return r,theta, phi
    
     
    
    
def distance_correction_from_perturbed_metric(parameters,xobs,yobs,zobs,xop,yop,zop):
    
    res = 10
    f = sp.linspace(0,1,res)
    
    d = sp.array([f*(xop-xobs)+xobs, f*(yop-yobs)+yobs, f*(zop-zobs)+zobs])

    Ng = len(parameters.potential)    
    d_grid = d*Ng/parameters.boxsize-1/2.
    # For the sake of the periodic boundaries
    d_grid[d_grid < 0] = d_grid[d_grid < 0]+Ng
    d_grid[d_grid > Ng] = d_grid[d_grid > Ng]-Ng
 
    psi = ndi.map_coordinates(parameters.potential,d_grid,mode='nearest')
    tck = interpolate.splrep(f,sp.sqrt(1-2*psi),s=0)
    psi_int = interpolate.splint(0,1,tck)

    
    return psi_int
    

    print "len(chosen_halos) = ", len(chosen_halos)

    position_CoM, position_CoM_op, velocity_CoM, total_mass = \
        center_of_mass(parameters,observer_position,chosen_halos)        
                 
    r_CoM, theta_CoM, phi_CoM = spherical_coordinates(parameters,observer_position,position_CoM)



        

def choose_halos(parameters,observed_halos):

    if parameters.observed_halos == 'all':
        chosen_halos = observed_halos    

    if parameters.observed_halos == 'mass_weighted':
        chosen_halos = mass_weighted_selection_of_halos(parameters,observed_halos)

    if parameters.observed_halos == 'random':
        chosen_halos = random_selection_of_halos(parameters,observed_halos)
        
    if parameters.observed_halos == 'specified_mass':
        chosen_halos = specified_mass_selection_of_halos(parameters,observed_halos)
        
    if parameters.observed_halos == 'centered_around_massive_halo':    
        chosen_halos = find_halos_around_massive_halo(parameters,observed_halos) 
        
    return chosen_halos


def find_halos_around_massive_halo(parameters,observed_halos):
    halos_around_halo = []
    central_strip_halos = [observed_halo for observed_halo in observed_halos\
                        if (observed_halo.r-parameters.bindistances[0] > parameters.min_dist)\
                        & (parameters.bindistances[-1]-observed_halo.r > parameters.min_dist)]            

    if len(central_strip_halos) == 0:
        print "No central strip halos for this observer"
        return 0

    # Since the halos are sorted according to mass, this should be the most massive halo                    
    center_halo = central_strip_halos[-1]

    for central_strip_halo in central_strip_halos:   
        if linalg.norm(center_halo.position_op-central_strip_halo.position_op)< parameters.min_dist:
            halos_around_halo.append(central_strip_halo)


    return halos_around_halo


   
def specified_mass_selection_of_halos(parameters,observed_halos):
    
    chosen_halos = []
    chosen_halos_of_specified_mass = [observed_halo for observed_halo in observed_halos\
                                    if (parameters.SN_mass_min < observed_halo.mass)
                                    & (observed_halo.mass < parameters.SN_mass_max)]
    random.seed(0)
    if len(chosen_halos_of_specified_mass) == 0:
        print "No halos to observe"
    else:
        while len(chosen_halos) < parameters.number_of_SNe:
            chosen_halo = random.choice(chosen_halos_of_specified_mass)
            chosen_halos.append(chosen_halo)
            
    return chosen_halos
            
        
    
def random_selection_of_halos(parameters,observed_halos):
    
    chosen_halos = []
    random.seed(0)
    for n in range(parameters.number_of_SNe):
        chosen_halo = random.choice(observed_halos)
        chosen_halos.append(chosen_halo)        
 
    return chosen_halos

def mass_weighted_selection_of_halos(parameters,observed_halos):
    
    total_mass = sp.sum([observed_halo.mass for observed_halo in observed_halos])
    
    cum_mass = 0
    cum_masses = []
    
    for observed_halo in observed_halos:
        mass = observed_halo.mass
        cum_mass = cum_mass+mass
        cum_masses.append(cum_mass/total_mass)
        
    
    chosen_halos = []
    random.seed(0)
    for n in range(parameters.number_of_SNe):
        rnd = random.random()
        for chosen_halo_number,chosen_halo in enumerate(chosen_halos):
            if cum_masses[chosen_halo_number] >= rnd:
                chosen_halos.append(chosen_halo)
                break
            
    return chosen_halos
        
        
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


    halos = [None]*sp.sum(Npart)
    for particle_type in [0,1]:
        N = int(Npart[particle_type])
        massunit = 1e10 #Msun/h
        M = sp.double(Massarr[particle_type])*massunit
        
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
            
            halo = hc.Halo(position,velocity,M,ID,ID_host)
            halos[p] = halo
    
    parameters.halos = sp.array(halos)    
    return 1

    
    
   
                    
      
