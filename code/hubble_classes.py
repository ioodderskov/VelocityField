from __future__ import division
import hubble_functions as hf
import powerspectrum_functions as pf
import scipy as sp
import yaml
import healpy as hp 
import gravitational_instability as gi
import resource
import pdb



class Parameters:
    def __init__(self,parameterfile):
        self.parameterfile = parameterfile

        # Loads parameters
        with open(parameterfile, 'r') as f:
            param = yaml.load(f)
        
        self.path = param["path"]
        self.halocatalogue_file = self.path+param["halocatalogue_file"]
        self.hubblefile = self.path+param["hubblefile"]
        self.CoDECS = int(param["CoDECS"])
        if self.CoDECS:
            self.CoDECShosts_file = self.path+param["CoDECShosts_file"]
        self.gridfile = self.path+param["gridfile"]
        
        self.parallel_processing = int(param["parallel_processing"])
        
        self.snapshot = int(param["snapshot"])
        self.snapshot_file = self.path+param["snapshot_file"]


        self.use_snapshot_for_background = int(param["use_snapshot_for_background"])
        self.use_grid = int(param["use_grid"])        

        self.use_CoM = int(param["use_CoM"])
        self.correct_for_peculiar_velocities = int(param["correct_for_peculiar_velocities"])
        self.survey_radius = sp.double(param["survey_radius"])
        self.min_dist = sp.double(param["min_dist"])
        
        self.use_local_velocity = int(param["use_local_velocity"])
        self.radius_local_group = sp.double(param["radius_local_group"])
        
        self.observer_choice = param["observer_choice"]
        self.observerfile = param["observerfile"]
        self.observer_indices_file = self.path+param["observer_indices_file"]
        self.number_of_observers = int(param["number_of_observers"])
        self.host_min_m = sp.double(param["host_min_m"])
        self.host_max_m = sp.double(param["host_max_m"])
        self.sub_min_m = sp.double(param["sub_min_m"])
        self.sub_max_m = sp.double(param["sub_max_m"])

        
        self.observed_halos = param["observed_halos"]
        self.SN_mass_min = sp.double(param["SN_mass_min"])
        self.SN_mass_max = sp.double(param["SN_mass_max"])
        
        self.mind = sp.double(param["mind"])
        self.maxd = sp.double(param["maxd"])
        self.width = sp.double(param["width"])
        self.number_of_SNe = int(param["number_of_SNe"])
        self.boxsize = sp.double(param["boxsize"])
        self.omegam = sp.double(param["omegam"])
        self.h = sp.double(param["h"])
        self.number_of_cones = int(param["number_of_cones"])
        self.skyfraction = sp.double(param["skyfraction"])
        self.max_angular_distance = sp.arccos(1-2*self.skyfraction)
        
        self.calculate_std_of_deviation = int(param["calculate_std_of_deviation"])
        self.calculate_hubble_constants = int(param["calculate_hubble_constants"])
        self.calculate_redshiftdistribution = int(param["calculate_redshiftdistribution"])
        self.make_hubblediagram = int(param["make_hubblediagram"])
        self.map_velocityfield = int(param["map_velocityfield"])
        self.calculate_powerspectra = int(param["calculate_powerspectra"])
        if self.calculate_powerspectra:
            self.powerspectrafile = self.path+param["powerspectrafile"]
        self.calculate_pairwise_velocities = int(param["calculate_pairwise_velocities"])

        
        self.distances_from_perturbed_metric = int(param["distances_from_perturbed_metric"])
        if self.distances_from_perturbed_metric:
            self.potential_file = self.path+param["potential_file"]
        
        self.vary_number_of_SNe = int(param["vary_number_of_SNe"])
        self.min_number_of_SNe = int(param["min_number_of_SNe"])
        self.max_number_of_SNe = int(param["max_number_of_SNe"])
        self.step_number_of_SNe = int(param["step_number_of_SNe"])
        self.numbers_of_SNe = range(self.min_number_of_SNe,self.max_number_of_SNe+1,self.step_number_of_SNe)
        
        self.nside = int(param["nside"])
        self.lmax = int(param["lmax"])
        self.smooth_map = int(param["smooth_map"])
        self.smooth_largest_hole = int(param["smooth_largest_hole"])
        self.preset_smoothinglength = int(param["preset_smoothinglength"])
        self.smoothing_fwhm = sp.double(param["smoothing_fwhm"])

        self.badval = 1e15
        self.unseen = 0
                
        
        if self.distances_from_perturbed_metric:
            potential_from_file = sp.loadtxt(self.potential_file)
            grid = int(sp.ceil(len(potential_from_file)**(1/3)))
            self.potential = sp.zeros((grid,grid,grid))
            row = 0
            for i in range(grid):
                for j in range(grid):
                    for k in range(grid):
                        self.potential[i][j][k] = potential_from_file[row][3]
                        row = row+1

            self.potential_min = self.potential.min()
            self.potential_max = self.potential.max()
            
        self.bindistances = hf.calculate_bindistances(self.mind,self.maxd,self.width)
        
        self.vary_skyfraction = int(param["vary_skyfraction"])
        self.fraction_start = sp.double(param["fraction_start"])
        self.fraction_stop = sp.double(param["fraction_stop"])
        self.fraction_step = sp.double(param["fraction_step"])
        self.skyfractions = sp.linspace(self.fraction_start,self.fraction_stop,1/self.fraction_step)

        self.use_lightcone = int(param["use_lightcone"])
        self.halocatalogue_filebase = param["halocatalogue_filebase"]
    

        self.test_isotropy = int(param["test_isotropy"])
        nside = 2
        self.number_of_directions = hp.nside2npix(nside)
        self.directions = hp.pix2ang(nside,range(self.number_of_directions))

        self.assign_to_grid = int(param["assign_to_grid"])
#        self.velocities_on_grid = int(param["velocities_on_grid"])
        self.Ng = int(param["Ng"])
        self.smoothing = int(param["smoothing"])    
        self.smoothing_radius = sp.double(param["smoothing_radius"])
        
        self.halos = []
        self.subhalos = []

        self.max_pairwise_distance = sp.double(param["max_pairwise_distance"])
        self.min_halo_mass = sp.double(param["min_halo_mass"])
        self.max_halo_mass = sp.double(param["max_halo_mass"])


class Halo:
    def __init__(self,position,velocity,mass,ID,ID_host):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.ID = ID
        self.ID_host = ID_host


class Observed_halo:
    def __init__(self,position,position_op,velocity,r,theta,phi,vr_peculiar,vr,ID,mass):
        self.position = position
        self.position_op = position_op
        self.velocity = velocity
        self.observed_velocity = []
        self.velocity_correction = []
        self.r = r
        self.theta = theta
        self.phi = phi
        self.vr = vr
        self.vr_peculiar = vr_peculiar
        self.ID = ID
        self.mass = mass
        
        
class Cone:
    def __init__(self,theta,phi,halos_in_cone):
        self.theta = theta
        self.phi = phi
        self.halos_in_cone = halos_in_cone


class Observer:
    def __init__(self,halo_number,position,velocity,mass):
        self.halo_number = halo_number
        self.position = position
        self.velocity = velocity 
        self.mass = mass
        self.local_position = []
        self.local_position_op = []
        self.local_velocity = []
        self.local_velocity_correction = []
        self.chosen_halos = []
        self.cones = []
        self.Hubbleconstants = []
        self.ls = []
        self.cls = []
        self.vrmap = []
        self.skyfraction = []
        self.radius_of_greatest_hole = []
        self.observed_radial_peculiar_velocities = []
        self.radial_distances = []
        self.pair_masses = []
        

 
        
    def observe(self, parameters,particles):

        if parameters.use_lightcone:
            halocatalogue_file = parameters.halocatalogue_filebase + '_' + str(self.observer_number)
            halocatalogue = hf.load_halocatalogue(parameters,halocatalogue_file)
            halos = hf.initiate_halos(parameters,halocatalogue)


        observed_halos = []
        local_halos = []
        
        survey_positions = []
        survey_masses = []
          
        for h in parameters.halos:
        
            position_op = hf.periodic_boundaries(parameters,self.position,h.position)

             
            r, theta, phi = hf.spherical_coordinates(parameters,self.position,
                                                position_op)
     

            if parameters.use_local_velocity:
                if r < parameters.radius_local_group:
                    local_halo = Observed_halo(h.position,position_op,h.velocity,r,theta,phi,[],[],h.ID,h.mass)
                    local_halos.append(local_halo)
                    
            if parameters.correct_for_peculiar_velocities:
                if not parameters.use_snapshot_for_background:
                    if r < parameters.survey_radius:
                        survey_positions.append(h.position)
                        survey_masses.append(h.mass)
            
            if r < parameters.mind or r > parameters.maxd:
                continue
            
            if (parameters.vary_skyfraction == 0) & (parameters.test_isotropy == 0):
                theta_max = sp.arccos(1-2*parameters.skyfraction)
                if theta > theta_max:
                    continue
            

            vr_peculiar = sp.dot(position_op-self.position,h.velocity)/r
            vr = vr_peculiar + r*100
            ID = h.ID
            mass = h.mass
            position = h.position
            velocity = h.velocity
            
            observed_halo = Observed_halo(position,position_op,velocity,r,theta,phi,vr_peculiar,vr,ID,mass)
            observed_halos.append(observed_halo)


        if len(observed_halos) == 0:
            print "No observed halos for this observer"
            return 0

        chosen_halos = hf.choose_halos(parameters,observed_halos)

        if parameters.use_CoM:
            
            position_CoM,position_CoM_op,velocity_CoM,r_CoM,theta_CoM,phi_CoM,vr_peculiar_CoM,vr_CoM,total_mass\
            = gi.determine_CoM_for_these_halos(parameters,self.position,chosen_halos)    
            
            #Overwriting the chosen halos with their center of mass motion
            chosen_halos = []
            CoM_halo = Observed_halo(position_CoM,position_CoM_op,velocity_CoM,r_CoM,theta_CoM,phi_CoM,
                                     vr_peculiar_CoM,vr_CoM,-1,total_mass)
            chosen_halos.append(CoM_halo)
    
        # If distances are to be calculated from the perturbed metric, this is only done for
        # the chosen halos.
        if parameters.distances_from_perturbed_metric:
            for chosen_halo in chosen_halos:
                position_op = chosen_halo.position
                xop,yop,zop = position_op[0],position_op[1],position_op[2]
                psi_int = hf.distance_correction_from_perturbed_metric(parameters,self.x,self.y,self.z,xop,yop,zop)
                chosen_halo.r = chosen_halo.r*psi_int







        if parameters.correct_for_peculiar_velocities:
			            
            if parameters.use_snapshot_for_background:
                survey_positions,survey_masses = self.survey(parameters,particles)
            
            for chosen_halo in chosen_halos:
                velocity_correction = gi.velocity_from_matterdistribution(parameters,\
                                        chosen_halo.position,survey_positions,survey_masses)
                chosen_halo.observed_velocity = chosen_halo.velocity
                chosen_halo.velocity_correction = velocity_correction
                chosen_halo.velocity = chosen_halo.velocity - velocity_correction

        if parameters.use_local_velocity:
            local_position_CoM, local_position_CoM_op,local_velocity_CoM,\
            dummy, dummy, dummy, dummy, dummy,local_total_mass\
            = gi.determine_CoM_for_these_halos(parameters,self.position,local_halos)

            self.local_position = local_position_CoM
            self.local_position_op = local_position_CoM_op
            self.local_velocity = local_velocity_CoM

            if parameters.correct_for_peculiar_velocities:
                
                local_velocity_correction = gi.velocity_from_matterdistribution(parameters,\
                                        self.local_position,survey_positions,survey_masses)
                
                self.local_velocity_correction = local_velocity_correction

     
	    self.chosen_halos = chosen_halos


        return 1



    def calculate_pairwise_velocities(self,parameters):
        
        for halo_number,halo in enumerate(parameters.halos):
			
            if halo_number <= self.halo_number:
                continue
            
            position_op = hf.periodic_boundaries(parameters,self.position,halo.position)             
            pair_mass = sp.array([self.mass,halo.mass])
            r, theta, phi = hf.spherical_coordinates(parameters,self.position,
                                                position_op)

            if r > parameters.max_pairwise_distance:
                continue
            
            relative_velocity = halo.velocity-self.velocity
            vr_relative = sp.dot(position_op-self.position,relative_velocity)/r
            #print "vr_relative = ", vr_relative
            self.observed_radial_peculiar_velocities.append(vr_relative)
            self.radial_distances.append(r)
            self.pair_masses.append(pair_mass)
            #print "pair_mass = ", pair_mass
            
        return 1
            
            
        



    def survey(self,parameters,particles):


        positions = []
        masses = []
          
        for p in particles:
        
            x,y,z = p.position[0],p.position[1],p.position[2]
             
            position_op = hf.periodic_boundaries(parameters,self.position,p.position)


            r, theta, phi = hf.spherical_coordinates(parameters,self.position,
                                                position_op)

  
            if r > parameters.survey_radius:
                continue
            
            position = sp.array([x,y,z])
            positions.append(position)

            
            mass = p.mass
            masses.append(mass)

        if len(masses) == 0:
            print "No particles in survey for this observer"
            return 0

        return positions, masses
            
            
            
    def do_hubble_analysis(self,parameters):
        
        if len(self.chosen_halos) == 0:
            print "No observed halos for this observer"
            self.Hubbleconstants = sp.ones_like(parameters.bindistances)*sp.nan
            return 0
                
        if parameters.test_isotropy:
            self.Hubbleconstants, self.cones \
                = hf.calculate_Hs_for_varying_directions(parameters,self.chosen_halos)
        
        elif parameters.vary_skyfraction:
            self.Hubbleconstants = hf.calculate_Hs_for_varying_skyfractions(parameters,self.chosen_halos)
    
        elif parameters.vary_number_of_SNe:
            self.Hubbleconstants = hf.calculate_Hs_for_varying_number_of_SNe(parameters,self.chosen_halos)

            
        else:
            self.Hubbleconstants = hf.calculate_Hs_for_these_observed_halos(parameters,self.chosen_halos)
        
        return 1
                    
    def calculate_powerspectra(self,parameters):
        
        thetas = [observed_halo.theta for observed_halo in self.chosen]
        phis = [observed_halo.phi for observed_halo in self.chosen]
        vrs_peculiar = [observed_halo.vr_peculiar for observed_halo in self.chosen]
        
        vrmap = pf.create_map(parameters,thetas,phis,vrs_peculiar) 
        vrmap = pf.fill_empty_entries(parameters,vrmap)
        
        if parameters.smooth_map:
            vrmap = pf.smooth_map(parameters,vrmap)
        
        self.ls, self.cls = pf.do_harmonic_analysis(parameters,vrmap)
        self.vrmap = vrmap
        



    
