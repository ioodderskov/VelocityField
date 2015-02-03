from __future__ import division
import hubble_functions as hf
import powerspectrum_functions as pf
import scipy as sp
import yaml
import pdb



class Parameters:
    def __init__(self,parameterfile):

        # Loads parameters
        with open(parameterfile, 'r') as f:
            param = yaml.load(f)
        
        self.halocatalogue_file = param["halocatalogue_file"]
        self.hubblefile = param["hubblefile"]
        self.powerspectrafile = param["powerspectrafile"]
        
        self.parallel_processing = int(param["parallel_processing"])
        
        self.snapshot = int(param["snapshot"])
        self.snapshot_file = param["snapshot_file"]
        
        self.observer_choice = param["observer_choice"]
        self.observerfile = param["observerfile"]
        self.number_of_observers = int(param["number_of_observers"])
        self.host_min_m = sp.double(param["host_max_m"])
        self.host_max_m = sp.double(param["host_min_m"])
        self.sub_min_m = sp.double(param["sub_min_m"])
        self.sub_max_m = sp.double(param["sub_max_m"])
        
        self.observed_halos = param["observed_halos"]
        
        self.mind = sp.double(param["mind"])
        self.maxd = sp.double(param["maxd"])
        self.width = sp.double(param["width"])
        self.number_of_SNe = int(param["number_of_SNe"])
        self.boxsize = sp.double(param["boxsize"])
        self.number_of_cones = int(param["number_of_cones"])
        self.skyfraction = sp.double(param["skyfraction"])
        
        self.calculate_std_of_deviation = int(param["calculate_std_of_deviation"])
        self.calculate_hubble_constants = int(param["calculate_hubble_constants"])
        self.calculate_redshiftdistribution = int(param["calculate_redshiftdistribution"])
        self.make_hubblediagram = int(param["make_hubblediagram"])
        self.map_velocityfield = int(param["map_velocityfield"])
        self.calculate_powerspectra = int(param["calculate_powerspectra"])
        
        self.distances_from_perturbed_metric = int(param["distances_from_perturbed_metric"])
        self.potential_file = param["potential_file"]
        
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


class Halo:
    def __init__(self,position,velocity,mass,ID,ID_host,index):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.ID = ID
        self.ID_host = ID_host
        self.index = index



class Observed_halo:
    def __init__(self,r,theta,phi,vr_peculiar,vr,ID):
        self.r = r
        self.theta = theta
        self.phi = phi
        self.vr = vr
        self.vr_peculiar = vr_peculiar
        self.ID = ID




class Observer:
    def __init__(self,observer_number,position):
        self.observer_number = observer_number
        self.x = position[0]
        self.y = position[1]
        self.z = position[2]
        self.observed_halos = []
        self.Hubbleconstants = []
        self.ls = []
        self.cls = []
        self.vrmap = []
        self.skyfraction = []
        self.radius_of_greatest_hole = []
        

    
    def observe(self, parameters,halos):

        if parameters.use_lightcone:
            halocatalogue_file = parameters.halocatalogue_filebase + '_' + str(self.observer_number)
            halocatalogue = hf.load_halocatalogue(halocatalogue_file)
            halos = hf.initiate_halos(parameters,halocatalogue)
#        else:
#            halocatalogue_file = parameters.halocatalogue_file
#            halocatalogue = hf.load_halocatalogue(halocatalogue_file)
#            halos = hf.initiate_halos(parameters,halocatalogue)
            


        candidates = []
        rs = []
        thetas = []
        phis = []
        vrs_peculiar = []
        vrs = []
        IDs = []
        xops = []
        yops = []
        zops = []

        
        for h in halos:
        
            x,y,z = h.position[0],h.position[1],h.position[2]
             
            [xop,yop,zop] = hf.periodic_boundaries(parameters,self.x,self.y,self.z,x,y,z)

             
            r, theta, phi = hf.spherical_coordinates(parameters,self.x, self.y, self.z,
                                                xop,yop,zop)
     
     
            if r < parameters.mind or r > parameters.maxd:
                continue
            
            if parameters.vary_skyfraction == 0:
                theta_max = sp.arccos(1-2*parameters.skyfraction)
                if theta > theta_max:
                    continue
            
                
            candidates.append(h)
            rs.append(r)
            thetas.append(theta)
            phis.append(phi)
            
            
            [vx,vy,vz] = h.velocity[[0,1,2]]

            vr_peculiar = ((xop-self.x)*vx+(yop-self.y)*vy+(zop-self.z)*vz)/r
            vrs_peculiar.append(vr_peculiar)
            
            vr = vr_peculiar + r*100
            vrs.append(vr)

            ID = h.ID
            IDs.append(ID)
            
            xops.append(xop)
            yops.append(yop)
            zops.append(zop)  
            
        

        selected_candidates = hf.select_candidates(parameters,halos,candidates)


        # If distances are to be calculated from the perturbed metric, this is only done for
        # the selected halos.
        if parameters.distances_from_perturbed_metric:
            for selected_candidate in selected_candidates:
                [xop,yop,zop] = [xops[selected_candidate],yops[selected_candidate],zops[selected_candidate]]
                psi_int = hf.distance_correction_from_perturbed_metric(parameters,self.x,self.y,self.z,xop,yop,zop)
                rs[selected_candidate] = rs[selected_candidate]*psi_int 
#                print "rper = ",rs[selected_candidate]

        
       
        for selected_candidate in selected_candidates:
            r = rs[selected_candidate]
            theta = thetas[selected_candidate]
            phi = phis[selected_candidate]
            vr_peculiar = vrs_peculiar[selected_candidate]
            vr = vrs[selected_candidate]
            ID = IDs[selected_candidate]
            
            self.observed_halos.append(Observed_halo(r,theta,phi,vr_peculiar,vr,ID))
            
            
            
    def do_hubble_analysis(self,parameters):
        
        if parameters.vary_skyfraction:
            self.Hubbleconstants = hf.calculate_Hs_for_varying_skyfractions(parameters,self.observed_halos)
    
        elif parameters.vary_number_of_SNe:
            self.Hubbleconstants = hf.calculate_Hs_for_varying_number_of_SNe(parameters,self.observed_halos)

            
        else:
            self.Hubbleconstants = hf.calculate_Hs_for_these_observed_halos(parameters,self.observed_halos)
        
    
                
        
        
    def calculate_powerspectra(self,parameters):
        
        thetas = [observed_halo.theta for observed_halo in self.observed_halos]
        phis = [observed_halo.phi for observed_halo in self.observed_halos]
        vrs_peculiar = [observed_halo.vr_peculiar for observed_halo in self.observed_halos]
#        pdb.set_trace()
        
        vrmap = pf.create_map(parameters,thetas,phis,vrs_peculiar) 
#        pdb.set_trace()
        self.radius_of_largest_hole = pf.find_largest_hole(parameters,vrmap)
        vrmap = pf.fill_empty_entries(parameters,vrmap)
        
        if parameters.smooth_map:
            vrmap = pf.smooth_map(parameters,vrmap)
        
        self.ls, self.cls = pf.do_harmonic_analysis(parameters,vrmap)
        self.vrmap = vrmap


    
