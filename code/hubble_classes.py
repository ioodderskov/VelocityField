from __future__ import division
import hubble_functions as hf
import powerspectrum_functions as pf
import scipy as sp
import yaml
import pdb
import healpy as hp 
import gravitational_instability as gi
import scipy.linalg as linalg


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
#        if self.snapshot:
        self.snapshot_file = self.path+param["snapshot_file"]


        self.use_snapshot_for_background = int(param["use_snapshot_for_background"])
        self.use_grid = int(param["use_grid"])        

        self.use_lonely_halo = int(param["use_lonely_halo"])
        self.use_CoM = int(param["use_CoM"])
        self.correct_for_peculiar_velocities = int(param["correct_for_peculiar_velocities"])
        self.survey_radius = sp.double(param["survey_radius"])
        self.min_dist = sp.double(param["min_dist"])
        self.cone_centered_on_halo = int(param["cone_centered_on_halo"])
        
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
        self.Ng = int(param["Ng"])
        self.smoothing = int(param["smoothing"])    
        self.smoothing_radius = sp.double(param["smoothing_radius"])


class Halo:
    def __init__(self,position,velocity,mass,ID,ID_host):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.ID = ID
        self.ID_host = ID_host




class Observed_halo:
    def __init__(self,position_op,r,theta,phi,vr_peculiar,vr,ID,mass,observed_velocity,velocity_correction):
        self.position_op = position_op
        self.r = r
        self.theta = theta
        self.phi = phi
        self.vr = vr
        self.vr_peculiar = vr_peculiar
        self.ID = ID
        self.mass = mass
        self.observed_velocity = observed_velocity
        self.velocity_correction = velocity_correction
        
class Cone:
    def __init__(self,theta,phi,halos_in_cone):
        self.theta = theta
        self.phi = phi
        self.halos_in_cone = halos_in_cone




class Observer:
    def __init__(self,observer_number,position):
        self.observer_number = observer_number
        self.position = position
        self.local_velocity = []
        self.local_velocity_correction = []
        self.observed_halos = []
        self.cones = []
        self.Hubbleconstants = []
        self.ls = []
        self.cls = []
        self.vrmap = []
        self.skyfraction = []
        self.radius_of_greatest_hole = []
        self.rho = []
        

 
        
    def observe(self, parameters,halos,particles):

        if parameters.use_lightcone:
            halocatalogue_file = parameters.halocatalogue_filebase + '_' + str(self.observer_number)
            halocatalogue = hf.load_halocatalogue(parameters,halocatalogue_file)
            halos = hf.initiate_halos(parameters,halocatalogue)


        candidates = []
        rs = []
        thetas = []
        phis = []
        vrs_peculiar = []
        vrs = []
        IDs = []
        positions_op = []
        masses = []
        
        local_halos = []
        
          
        for h in halos:
        
            position_op = hf.periodic_boundaries(parameters,self.position,h.position)

             
            r, theta, phi = hf.spherical_coordinates(parameters,self.position,
                                                position_op)
     

            if parameters.use_local_velocity:
                if r < parameters.radius_local_group:
                    local_halos.append(h)
            
            if r < parameters.mind or r > parameters.maxd:
                continue
            
            if (parameters.vary_skyfraction == 0) & (parameters.test_isotropy == 0) & (parameters.cone_centered_on_halo == 0):
                theta_max = sp.arccos(1-2*parameters.skyfraction)
                if theta > theta_max:
                    continue
            
                
            candidates.append(h)
            rs.append(r)
            thetas.append(theta)
            phis.append(phi)
            
            
            [vx,vy,vz] = h.velocity[[0,1,2]]


            vr_peculiar = sp.dot(position_op-self.position,h.velocity)/r
            vrs_peculiar.append(vr_peculiar)
            
            vr = vr_peculiar + r*100
            vrs.append(vr)

            ID = h.ID
            IDs.append(ID)
            
            positions_op.append(position_op)

            
            mass = h.mass
            masses.append(mass)

        masses = sp.array(masses)
        positions_op = sp.array(positions_op)
        if len(masses) == 0:
            print "No observed halos for this observer"
            return 0
        
        lonely_candidates = []
        if parameters.use_lonely_halo:
            central_strip_indices = sp.array(range(len(masses)))\
                                    [(rs-parameters.bindistances[0] > parameters.min_dist) \
                                    & (parameters.bindistances[-1]-rs > parameters.min_dist)]
            for index,position_halo in zip(central_strip_indices,positions_op[central_strip_indices]):
                distances_to_other_candidates = sp.array([linalg.norm(position_halo-position)\
                                                for position in positions_op if linalg.norm(position_halo-position) != 0])
                if sp.amin(distances_to_other_candidates) > parameters.min_dist:
                    lonely_candidates.append(candidates[index])
                    candidates = lonely_candidates
                    break
            if len(lonely_candidates) == 0:
                return 0
                    


            
        elif parameters.cone_centered_on_halo:            
                        
            central_strip_indices = sp.array(range(len(masses)))\
                                    [(rs-parameters.bindistances[0] > parameters.min_dist) \
                                    & (parameters.bindistances[-1]-rs > parameters.min_dist)]
            central_strip_masses = masses[central_strip_indices]
            if len(central_strip_masses) == 0:
                print "No central strip halos for this observer"
                return 0
                
            index_center = central_strip_indices[central_strip_masses == max(central_strip_masses)][0]

            position_op_center = positions_op[index_center]

#            pdb.set_trace()

            candidates_around_center = []
            for candidate,position_op in zip(candidates,positions_op):   
                if linalg.norm(position_op_center-position_op) < parameters.min_dist:
                    candidates_around_center.append(candidate)
#            pdb.set_trace()        
            candidates = candidates_around_center



            
#            theta_center = thetas[index_center]
#            phi_center = phis[index_center]
#
#            candidates_around_center = []
#            for candidate,theta,phi in zip(candidates,thetas,phis):   
#                if hp.rotator.angdist([theta_center,phi_center],[theta,phi]) < theta_max:
##                    print "dir_halo = ", theta,phi
##                    print "dir_center = ", theta_center,phi_center
##                    print "ang_dist = ",hp.rotator.angdist([theta_center,phi_center],[theta,phi])
##                    print "theta_max = ", theta_max
#                    candidates_around_center.append(candidate)
#                    
#            candidates = candidates_around_center
                    
            

        if parameters.use_CoM:

            observer_position = self.position
            if parameters.correct_for_peculiar_velocities:
                if parameters.use_snapshot_for_background:
                    survey_positions,survey_masses = self.survey(parameters,particles)
                else:
                    survey_positions,survey_masses = self.survey(parameters,halos)
            position_CoM, \
            r_CoM, theta_CoM, phi_CoM,\
            vr_peculiar_CoM, vr_CoM, \
            total_mass, observed_velocity, velocity_correction,\
            local_velocity, local_velocity_correction = hf.determine_CoM_for_these_halos(parameters,survey_positions,survey_masses,observer_position,local_halos,candidates)
            self.observed_halos.append(Observed_halo(position_CoM,r_CoM,theta_CoM,phi_CoM,vr_peculiar_CoM,vr_CoM,-1,total_mass,observed_velocity,velocity_correction))
            self.local_velocity = local_velocity
            self.local_velocity_correction = local_velocity_correction
           
        else:
        
            if parameters.test_isotropy:
                halos_to_store = range(len(candidates))
            else:
                halos_to_store = hf.select_candidates(parameters,candidates)
    
    
            # If distances are to be calculated from the perturbed metric, this is only done for
            # the selected halos.
            if parameters.distances_from_perturbed_metric:
                for halo_to_store in halos_to_store:
                    xop,yop,zop = positions_op[halo_to_store,0],positions_op[halo_to_store,1],positions_op[halo_to_store,2]
                    psi_int = hf.distance_correction_from_perturbed_metric(parameters,self.x,self.y,self.z,xop,yop,zop)
                    rs[halo_to_store] = rs[halo_to_store]*psi_int 
    
            
           
            for halo_to_store in halos_to_store:
                r = rs[halo_to_store]
                theta = thetas[halo_to_store]
                phi = phis[halo_to_store]
                vr_peculiar = vrs_peculiar[halo_to_store]
                vr = vrs[halo_to_store]
                ID = IDs[halo_to_store]
                mass = masses[halo_to_store]
                xop,yop,zop = positions_op[halo_to_store,0],positions_op[halo_to_store,1],positions_op[halo_to_store,2]

            
                self.observed_halos.append(Observed_halo(position_op,r,theta,phi,vr_peculiar,vr,ID,mass,[],[]))

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
        
        if len(self.observed_halos) == 0:
            print "No observed halos for this observer"
            self.Hubbleconstants = sp.ones_like(parameters.bindistances)*sp.nan
            return 0
                
        if parameters.test_isotropy:
            self.Hubbleconstants, self.cones = hf.calculate_Hs_for_varying_directions(parameters,self.observed_halos)
        
        elif parameters.vary_skyfraction:
            self.Hubbleconstants = hf.calculate_Hs_for_varying_skyfractions(parameters,self.observed_halos)
    
        elif parameters.vary_number_of_SNe:
            self.Hubbleconstants = hf.calculate_Hs_for_varying_number_of_SNe(parameters,self.observed_halos)

            
        else:
            self.Hubbleconstants = hf.calculate_Hs_for_these_observed_halos(parameters,self.observed_halos)
        
        return 1
                    
    def calculate_powerspectra(self,parameters):
        
        thetas = [observed_halo.theta for observed_halo in self.observed_halos]
        phis = [observed_halo.phi for observed_halo in self.observed_halos]
        vrs_peculiar = [observed_halo.vr_peculiar for observed_halo in self.observed_halos]
        
        vrmap = pf.create_map(parameters,thetas,phis,vrs_peculiar) 
#        self.radius_of_largest_hole = pf.find_largest_hole(parameters,vrmap)
#        outputfile = 'outputfile.txt'
#        f = open(outputfile,'w')
#        f.write("radius of largest hole = %s" % self.radius_of_largest_hole)
#        f.close()
        vrmap = pf.fill_empty_entries(parameters,vrmap)
        
        if parameters.smooth_map:
            vrmap = pf.smooth_map(parameters,vrmap)
        
        self.ls, self.cls = pf.do_harmonic_analysis(parameters,vrmap)
        self.vrmap = vrmap


    
