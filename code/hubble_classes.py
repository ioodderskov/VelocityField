from __future__ import division
import hubble_functions as hf
import powerspectrum_functions as pf
import scipy as sp
import yaml
import pdb
import healpy as hp 



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
        if self.snapshot:
            self.snapshot_file = self.path+param["snapshot_file"]

        self.use_grid = int(param["use_grid"])        
        
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


class Halo:
    def __init__(self,position,velocity,mass,ID,ID_host,index):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.ID = ID
        self.ID_host = ID_host
        self.index = index



class Observed_halo:
    def __init__(self,xop,yop,zop,r,theta,phi,vr_peculiar,vr,ID,mass):
        self.xop = xop
        self.yop = yop
        self.zop = zop
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
    def __init__(self,observer_number,position):
        self.observer_number = observer_number
        self.x = position[0]
        self.y = position[1]
        self.z = position[2]
        self.observed_halos = []
        self.cones = []
        self.Hubbleconstants = []
        self.ls = []
        self.cls = []
        self.vrmap = []
        self.skyfraction = []
        self.radius_of_greatest_hole = []
        self.rho = []
        

 
    
    def observe(self, parameters,halos):

        if parameters.use_lightcone:
            halocatalogue_file = parameters.halocatalogue_filebase + '_' + str(self.observer_number)
            halocatalogue = hf.load_halocatalogue(parameters,halocatalogue_file)
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
        masses = []
        
          
        for h in halos:
        
            x,y,z = h.position[0],h.position[1],h.position[2]
             
            [xop,yop,zop] = hf.periodic_boundaries(parameters,self.x,self.y,self.z,x,y,z)

             
            r, theta, phi = hf.spherical_coordinates(parameters,self.x, self.y, self.z,
                                                xop,yop,zop)
     
             
            
            if r < parameters.mind or r > parameters.maxd:
                continue
            
            if (parameters.vary_skyfraction == 0) & (parameters.test_isotropy == 0):
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
            
            mass = h.mass
            masses.append(mass)
            
        
        if parameters.test_isotropy:
            halos_to_store = range(len(candidates))
        else:
            halos_to_store = hf.select_candidates(parameters,candidates)


        # If distances are to be calculated from the perturbed metric, this is only done for
        # the selected halos.
        if parameters.distances_from_perturbed_metric:
            for halo_to_store in halos_to_store:
                [xop,yop,zop] = [xops[halo_to_store],yops[halo_to_store],zops[halo_to_store]]
                psi_int = hf.distance_correction_from_perturbed_metric(parameters,self.x,self.y,self.z,xop,yop,zop)
                rs[halo_to_store] = rs[halo_to_store]*psi_int 
#                print "rper = ",rs[selected_candidate]

        
       
        for halo_to_store in halos_to_store:
            r = rs[halo_to_store]
            theta = thetas[halo_to_store]
            phi = phis[halo_to_store]
            vr_peculiar = vrs_peculiar[halo_to_store]
            vr = vrs[halo_to_store]
            ID = IDs[halo_to_store]
            mass = masses[halo_to_store]
            xop,yop,zop = xops[halo_to_store],yops[halo_to_store],zops[halo_to_store] 
            
            self.observed_halos.append(Observed_halo(xop,yop,zop,r,theta,phi,vr_peculiar,vr,ID,mass))

            
            
    def do_hubble_analysis(self,parameters):
        
        
        if parameters.test_isotropy:
            self.Hubbleconstants, self.cones = hf.calculate_Hs_for_varying_directions(parameters,self.observed_halos)
        
        elif parameters.vary_skyfraction:
            self.Hubbleconstants = hf.calculate_Hs_for_varying_skyfractions(parameters,self.observed_halos)
    
        elif parameters.vary_number_of_SNe:
            self.Hubbleconstants = hf.calculate_Hs_for_varying_number_of_SNe(parameters,self.observed_halos)

            
        else:
            self.Hubbleconstants = hf.calculate_Hs_for_these_observed_halos(parameters,self.observed_halos)
        
    
                    
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


    
