from __future__ import division
import hubble_functions as hf
import powerspectrum_functions as pf
import scipy as sp
import yaml
import healpy as hp 
import gravitational_instability as gi
import resource
import pdb
import copy
import matplotlib.pyplot as plt



class Parameters:
    def __init__(self,parameterfile):
        self.parameterfile = parameterfile

        # Loads parameters
        with open(parameterfile, 'r') as f:
            param = yaml.load(f)
            
        # The code will stop if an essential parameter is missing from the parameter file.
        # If a non-essential parameter is missing, it will be loaded as "default"
        default = 0

        # Path and choices for how and what to run        
        self.path = param["path"]
        self.parallel_processing = int(param["parallel_processing"])
        self.data_type = param["data_type"]
        self.data_to_observe = param["data_to_observe"]
        
        # Identify and save galaxies using a HOD?
        self.use_HOD = int(param.get("use_HOD",default))
        self.number_of_particle_files = int(param.get("number_of_particle_files",default))
        self.particle_file_base = str(param.get("particle_file_base",default))
        # Parameters for HOD
        self.logMmin = sp.double(param.get("logMmin",default))
        self.sigma_logM = sp.double(param.get("sigma_logM",default))
        self.logM0 = sp.double(param.get("logM0",default))
        self.logM1_prime = sp.double(param.get("logM1_prime",default))
        self.alpha = sp.double(param.get("alpha",default))
        self.alpha_c = sp.double(param.get("alpha_c",default))
        self.alpha_s = sp.double(param.get("alpha_s",default))
        
        # Halo catalogues
        self.halocatalogue_file = self.path+str(param.get("halocatalogue_file",default))
        self.CoDECShosts_file = self.path+str(param.get("CoDECShosts_file",default))

        
        # Parameters for the simulation/catalogue
        self.boxsize = sp.double(param.get("boxsize",default))
        self.omegam = sp.double(param.get("omegam",default))
        self.omega_b = sp.double(param.get("omega_b",default))
        self.h = sp.double(param.get("h",default))
        self.beta_c = sp.double(param.get("beta_c",default))

        # Division in bins
        self.mind = sp.double(param["mind"])
        self.maxd = sp.double(param["maxd"])
        self.width = sp.double(param["width"])
        self.skyfraction = sp.double(param.get("skyfraction",default))
        self.bindistances = hf.calculate_bindistances(self.mind,self.maxd,self.width)


        # Use lightcones?
        self.use_lightcone = int(param.get("use_lightcone",default))
        self.halocatalogue_filebase = param.get("halocatalogue_filebase",default)
        
        # Snapshots
        self.snapshot_file = self.path+str(param.get("snapshot_file",default))

        # Grid file (both read and write)
        self.gridfile = self.path+str(param.get("gridfile",default))

        # Outputfiles
        self.hubblefile = self.path+str(param.get("hubblefile",default))
        self.powerspectrafile = self.path+str(param.get("powerspectrafile",default))

        # Correction for peculiar velocities (using first order pertubation theory)
        self.use_snapshot_for_background = int(param.get("use_snapshot_for_background",default))
        self.correct_for_peculiar_velocities = int(param.get("correct_for_peculiar_velocities",default))
        self.survey_radius = sp.double(param.get("survey_radius",default))
        self.min_dist = sp.double(param.get("min_dist",default))
        self.plot_velocity_field = int(param.get("plot_velocity_field",default))
        self.plottet_velocity_field_file = self.path+str(param.get("plottet_velocity_field_file",default))
        self.plot_velocity_field_n = 0

        # Cepheids, tracking of local motion        
        self.use_local_velocity = int(param.get("use_local_velocity",default))
        self.radius_local_group = sp.double(param.get("radius_local_group",default))
        
        # Choices for the observers        
        self.observer_choice = param["observer_choice"]
        self.observerfile = self.path+str(param.get("observerfile",default))
        self.observer_indices_file = self.path+str(param.get("observer_indices_file",default))
        self.number_of_observers = int(param["number_of_observers"])
        self.host_min_m = sp.double(param.get("host_min_m",default))
        self.host_max_m = sp.double(param.get("host_max_m",default))
        self.sub_min_m = sp.double(param.get("sub_min_m",default))
        self.sub_max_m = sp.double(param.get("sub_max_m",default))

        # Choices for the observed halos
        self.use_CoM = int(param.get("use_CoM",default))
        self.number_of_SNe = int(param.get("number_of_SNe",default)) # Not used if use_CoM = 1 (?)
        self.observed_halos = param["observed_halos"]
        self.SN_mass_min = sp.double(param.get("SN_mass_min",default))
        self.SN_mass_max = sp.double(param.get("SN_mass_max",default))


        # Assigning values to a grid
        self.assign_to_grid = int(param.get("assign_to_grid",default))
        self.velocities_on_grid = int(param.get("velocities_on_grid",default))
        self.Ng = int(param.get("Ng",default))
        self.smoothing = int(param.get("smoothing",default))    
        self.smoothing_radius = sp.double(param.get("smoothing_radius",default))
        self.reduced_box = int(param.get("reduced_box",default))
        self.reduced_boxsize = sp.double(param.get("reduced_boxsize",default))


        
        # For investigation of isotropy
        self.test_isotropy = int(param.get("test_isotropy",default))
        nside = 2
        self.number_of_directions = hp.nside2npix(nside)
        self.directions = hp.pix2ang(nside,range(self.number_of_directions))
        self.number_of_cones = int(param.get("number_of_cones",default))
        self.max_angular_distance = sp.arccos(1-2*self.skyfraction)
        self.vary_skyfraction = int(param.get("vary_skyfraction",default))
        self.fraction_start = sp.double(param.get("fraction_start",default))
        self.fraction_stop = sp.double(param.get("fraction_stop",default))
        self.fraction_step = sp.double(param.get("fraction_step",default))

        if self.vary_skyfraction:
            self.skyfractions = sp.linspace(self.fraction_start,self.fraction_stop,1/self.fraction_step)
        
        # Stuff from the Hubble2013 project
        self.calculate_hubble_constants = int(param.get("calculate_hubble_constants",default))
        self.calculate_std_of_deviation = int(param.get("calculate_std_of_deviation",default))
        self.calculate_redshiftdistribution = int(param.get("calculate_redshiftdistribution",default))
        self.make_hubblediagram = int(param.get("make_hubblediagram",default))

        # Angular powerspectra for the radial peculiar velocities
        self.map_velocityfield = int(param.get("map_velocityfield",default))
        self.masked_map = int(param.get("masked_map",default))
        self.smoothed_map = int(param.get("smoothed_map",default))
        self.beam_fwhm = sp.double(param.get("beam_fwhm",default))
        self.calculate_powerspectra = int(param.get("calculate_powerspectra",default))
        self.nside = int(param.get("nside",default))
        self.lmax = int(param.get("lmax",default))
#        self.smooth_map = int(param.get("smooth_map",default))
#        self.smooth_largest_hole = int(param.get("smooth_largest_hole",default))
#        self.preset_smoothinglength = int(param.get("preset_smoothinglength",default))
#        self.smoothing_fwhm = sp.double(param.get("smoothing_fwhm",default))
        self.badval = 1e15
        self.unseen = 1e16


        
        # For the pairwise velocity distribution
        self.calculate_pairwise_velocities = int(param.get("calculate_pairwise_velocities",default))
        self.max_pairwise_distance = sp.double(param.get("max_pairwise_distance",default))
        self.min_halo_mass = sp.double(param.get("min_halo_mass",default))
        self.max_halo_mass = sp.double(param.get("max_halo_mass",default))

        # For calculating distances from the perturbed metric
        self.distances_from_perturbed_metric = int(param.get("distances_from_perturbed_metric",default))
        self.potential_file = self.path+str(param.get("potential_file",default))
        self.potential, self.potential_min, self.potential_max = hf.potential_on_grid(self.distances_from_perturbed_metric,self.potential_file)
        
        # For varying the number of observed halos
        self.vary_number_of_SNe = int(param.get("vary_number_of_SNe",default))
        self.min_number_of_SNe = int(param.get("min_number_of_SNe",default))
        self.max_number_of_SNe = int(param.get("max_number_of_SNe",default))
        self.step_number_of_SNe = int(param.get("step_number_of_SNe",default))
        if self.vary_number_of_SNe:
            self.numbers_of_SNe = range(self.min_number_of_SNe,self.max_number_of_SNe+1,self.step_number_of_SNe)

        # Initiating the list for halos and subhalos
        self.halos = []
#        self.subhalos = []
        # Initiating the list for galaxies
        self.galaxies = sp.empty((0,1))
        self.grid = []
        # Particle files:
        self.particle_files = []




class Particle_file:
    def __init__(self,name):
        self.name = name
        self.halo_IDs_in_file = []

            


class Galaxy:
    def __init__(self,position,velocity,ID_host):
        self.position = position
        self.velocity = velocity
        self.ID_host = ID_host




class Halo:
    def __init__(self,position,velocity,mass,vrms,ID,ID_host):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.vrms = vrms
        self.ID = ID
        self.ID_host = ID_host
        self.particle_file = []
        self.number_of_particles = []
        self.central_galaxies = []
        self.satellite_galaxies = []


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
        self.chosen_halos = []
        self.local_position = []
        self.local_position_op = []
        self.local_velocity = []
        self.local_velocity_correction = []
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
        
        if parameters.data_to_observe == 'halos':
            halos = parameters.halos
        if parameters.data_to_observe == 'grid':
            halos = parameters.grid
          
        for h in halos:
        
            position_op = hf.periodic_boundaries(parameters,self.position,h.position)

             
            r, theta, phi = hf.spherical_coordinates(parameters,self.position,
                                                position_op)
                                        

            if parameters.use_local_velocity:
                if r < parameters.radius_local_group:
                    local_halo = Observed_halo(h.position,position_op,h.velocity,r,theta,phi,[],[],h.ID,h.mass)
                    local_halos.append(local_halo)
            
#            print "the mass of this observer is", self.mass
            if parameters.correct_for_peculiar_velocities:
                if not parameters.use_snapshot_for_background:
                    if r < parameters.survey_radius:
#                        print "survey. r = ", r, "mass = ", h.mass
#                        if r < 1e-10:
#                            print "this is the observer halo, included in the survey"
                        survey_positions.append(h.position)
                        survey_masses.append(h.mass)
            
            if r < parameters.mind or r > parameters.maxd:
                continue

            
            if (parameters.vary_skyfraction == 0) & (parameters.test_isotropy == 0):
                theta_max = sp.arccos(1-2*parameters.skyfraction)
                if theta > theta_max:
#                    print "theta = ", theta, "and theta_max =", theta_max
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
            print "(observer.observe): No observed halos for this observer"
            return 0

        chosen_halos = hf.choose_halos(parameters,observed_halos)
        

        if parameters.use_CoM:
                        
            position_CoM,position_CoM_op,velocity_CoM,r_CoM,theta_CoM,phi_CoM,vr_peculiar_CoM,vr_CoM,total_mass\
            = gi.determine_CoM_for_these_halos(parameters,self.position,chosen_halos) 
            
            #Overwriting the chosen halos with their center of mass motion
            halos_around_massive_halo = copy.copy(chosen_halos)
            
            chosen_halos = []
            CoM_halo = Observed_halo(position_CoM,position_CoM_op,velocity_CoM,r_CoM,theta_CoM,phi_CoM,
                                     vr_peculiar_CoM,vr_CoM,-1,total_mass)
            chosen_halos.append(CoM_halo)
            print "the distance to the observed CoM-halo is r_CoM = ", r_CoM
    
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

            if parameters.plot_velocity_field:
                parameters.plot_velocity_field_n += 1
                gi.plot_velocity_field(parameters,self,chosen_halos,
                                       halos_around_massive_halo,survey_positions,
                                       survey_masses)
                
                parameters.plot_velocity_field -= 1


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
            print "(do_hubble_analysis): No observed halos for this observer"
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
        
        def calculate_powerspectrum(halos_in_bin):
        
            thetas = [halo_in_bin.theta for halo_in_bin in halos_in_bin]
            phis = [halo_in_bin.phi for halo_in_bin in halos_in_bin]
            vrs_peculiar = [halo_in_bin.vr_peculiar for halo_in_bin in halos_in_bin]
            
            vrmap = pf.create_map(parameters,thetas,phis,vrs_peculiar) 
                    
            if parameters.masked_map:
                vrmap = pf.mask_map(parameters,vrmap)
                
            if parameters.smoothed_map:
                nside = parameters.nside
                approx_pixsize = hp.rotator.angdist(hp.pix2ang(nside,0),hp.pix2ang(nside,1))
                if parameters.beam_fwhm > approx_pixsize:
                    print "Warning: The smoothing length shouldn't be larger than the size of the pixels!"

                vrmap = hp.smoothing(vrmap,fwhm=parameters.beam_fwhm)
                
            
            ls, cls = pf.do_harmonic_analysis(parameters,vrmap)
            
            return ls, cls, vrmap
            
        for mind_bin, maxd_bin in zip(parameters.bindistances[0:-1],parameters.bindistances[1:]):

            halos_in_bin = [observed_halo for observed_halo in self.chosen_halos\
                            if (mind_bin < observed_halo.r) & (observed_halo.r < maxd_bin)]

            ls, cls, vrmap = calculate_powerspectrum(halos_in_bin)

            self.ls.append(ls)
            self.cls.append(cls)
            self.vrmap.append(vrmap)
       
        return 0
        



    
