from __future__ import division
import numpy as sp
import scipy.linalg as linalg
import pdb
import hubble_functions as hf
from functools import partial

periodic_boundaries = 1


def velocity_from_matterdistribution(parameters,halo_position,survey_positions,survey_masses):
    
    
    positions = sp.array(survey_positions)
    masses  = sp.array(survey_masses)

    meter = 1. # meter
    second = 1. # second
    kg = 1. # kg
    Mpc = 1e6*3.085678e16*meter
    Msun = 1.989e30*kg
    H0 = 100. *1e3*meter/(second*(Mpc/parameters.h)) # s**-1
    Gc = 6.67e-11 *meter**3/(kg*second**2)
    rho_critical = 3.*H0**2/(8*sp.pi*Gc) # kg/m**3
    rho_critical_astro = rho_critical / (Msun/parameters.h) * (Mpc/parameters.h)**3


    survey_volume = 4/3*sp.pi*parameters.survey_radius**3
    rho_survey = sp.sum(masses)/survey_volume
    rho_mean = rho_critical_astro*parameters.omegam

    print "rho_survey = ", rho_survey, "and rho_critical_astro*parameters.omegam = ", rho_critical_astro*parameters.omegam


    #PERIODISKE GRAENSER:
    if periodic_boundaries:
    
        positions_oh = []
        for index,position in enumerate(positions):

            xop,yop,zop = hf.periodic_boundaries(parameters,halo_position,position)
            positions_oh.append(sp.array([xop,yop,zop]))

        positions_oh = sp.array(positions_oh)
    

    individual_contributions = sp.array([masses[i] \
                                        *(positions_oh[i]-halo_position) \
                                        /linalg.norm(positions_oh[i]-halo_position)**3 \
                                        for i in range(len(masses)) \
                                        if (linalg.norm(positions_oh[i]-halo_position) > parameters.min_dist)])
        
    
    if len(individual_contributions) == 0:
        print "No halos was found between min_dist and max_dist"
        v_from_matterdistribution = sp.array([0,0,0])

    else:

        H0_astro = 100.
        f = parameters.omegam**0.55
    
        vx_from_matterdistribution = f*H0_astro/(4*sp.pi*rho_mean)*sp.sum(individual_contributions[:,0])#+f*H0_astro/3*(halo_position[0]-observer_position[0])
        vy_from_matterdistribution = f*H0_astro/(4*sp.pi*rho_mean)*sp.sum(individual_contributions[:,1])#+f*H0_astro/3*(halo_position[1]-observer_position[1])
        vz_from_matterdistribution = f*H0_astro/(4*sp.pi*rho_mean)*sp.sum(individual_contributions[:,2])#+f*H0_astro/3*(halo_position[2]-observer_position[2])
        
        v_from_matterdistribution = sp.array([vx_from_matterdistribution,vy_from_matterdistribution,vz_from_matterdistribution])
    
    return v_from_matterdistribution
    
    
def center_of_mass(parameters,observer_position,chosen_halos):
    masses = sp.array([chosen_halo.mass for chosen_halo in chosen_halos])    
    positions = sp.array([chosen_halo.position for chosen_halo in chosen_halos])
    positions_op = sp.array([chosen_halo.position_op for chosen_halo in chosen_halos])    
    velocities = sp.array([chosen_halo.velocity for chosen_halo in chosen_halos])
    position_CoM = sp.average(positions,axis=0,weights=masses)
    position_CoM_op = sp.average(positions_op,axis=0,weights=masses)
    velocity_CoM = sp.average(velocities,axis=0,weights=masses)
    
    total_mass = sp.sum(masses)

    return position_CoM, position_CoM_op, velocity_CoM, total_mass
    
    
    
def determine_CoM_for_these_halos(parameters,observer_position,chosen_halos):

    position_CoM, position_CoM_op, velocity_CoM, total_mass = \
        center_of_mass(parameters,observer_position,chosen_halos)        
                 
    r_CoM, theta_CoM, phi_CoM = hf.spherical_coordinates(parameters,observer_position,position_CoM_op)

    vr_peculiar_CoM = sp.dot(position_CoM_op-observer_position,velocity_CoM)/r_CoM
    
    vr_CoM = vr_peculiar_CoM + r_CoM*100
 

    return position_CoM,position_CoM_op,velocity_CoM,r_CoM,theta_CoM,phi_CoM,\
            vr_peculiar_CoM,vr_CoM,total_mass 



    