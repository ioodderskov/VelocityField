from __future__ import division
import numpy as sp
import scipy.linalg as linalg
import pdb
import hubble_functions as hf
from functools import partial
import matplotlib.pyplot as plt
import copy

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

#    print "the distances and masses in the survey, with respect to the observed halo, are"
#    for i in range(len(masses)):
#        print "r = ", linalg.norm(positions_oh[i]-halo_position), "m = ", masses[i]
         

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

def translate_into_box(parameters,positions):
    positions[positions < 0] = positions[positions < 0] + parameters.boxsize
    positions[positions > parameters.boxsize] = positions[positions > parameters.boxsize]\
                                                - parameters.boxsize
    return positions

        

def plot_velocity_field(parameters,observer,chosen_halos,halos_around_massive_halo,
                        survey_positions,survey_masses):

    obs_x = observer.position[0]
    obs_y = observer.position[1]

    boxplot = 25    
    scale = 2000
    linear_halo_bias = 1
    
    plt.figure()
    ax = plt.gca()
    plt.xlabel("x [Mpc/h]")
    plt.ylabel("y [Mpc/h]")
    ax.set_xlim([obs_x-boxplot,obs_x+boxplot])
    ax.set_ylim([obs_y-boxplot,obs_y+boxplot])
    ax.set_aspect('equal', adjustable='box')
    
    positions_of_chosen_halos = sp.array([chosen_halo.position_op\
            for chosen_halo in chosen_halos])
    positions_of_halos_around_massive_halo = sp.array([halo_around_massive_halo.position_op\
            for halo_around_massive_halo in halos_around_massive_halo])
    masses_of_halos_around_massive_halo = sp.array([halo_around_massive_halo.mass\
            for halo_around_massive_halo in halos_around_massive_halo])
                
#    positions_of_chosen_halos = translate_into_box(parameters,positions_of_chosen_halos)
#    positions_of_halos_around_massive_halo = translate_into_box(parameters,positions_of_halos_around_massive_halo)

    slap = parameters.min_dist


    survey_positions_op = []
    for index,survey_position in enumerate(survey_positions):

        xop,yop,zop = hf.periodic_boundaries(parameters,observer.position,survey_position)
        survey_positions_op.append(sp.array([xop,yop,zop]))

    survey_positions = sp.array(survey_positions_op)


    survey_positions_slap = sp.array([survey_position for survey_position in survey_positions\
            if sp.abs(survey_position[2]-positions_of_chosen_halos[0,2]) < slap]) 

               
    survey_masses_slap = sp.array([survey_mass\
            for survey_position, survey_mass in zip(survey_positions,survey_masses)\
            if sp.abs(survey_position[2]-positions_of_chosen_halos[0,2]) < slap])

#    mnorm = 2.5e15
#    m = masses/mnorm
#    ax.scatter(x,y,marker='o',c='b', s=10*m )

    mnorm = sp.amax(survey_masses_slap)
    m = survey_masses_slap/mnorm
    if len(survey_positions_slap) == 0:
        print "No nearby halos to plot"
    else:
        ax.scatter(survey_positions_slap[:,0],survey_positions_slap[:,1],
                   marker='.',c='k',s=500*m)
  
    ax.plot(obs_x,obs_y,'g*',markersize=15)
    thetas = sp.linspace(0,2*sp.pi,100)
    ax.plot(positions_of_chosen_halos[:,0]+sp.cos(thetas)*parameters.min_dist,
            positions_of_chosen_halos[:,1]+sp.sin(thetas)*parameters.min_dist,'b--')
    mc = masses_of_halos_around_massive_halo/mnorm
    ax.scatter(positions_of_halos_around_massive_halo[:,0],
            positions_of_halos_around_massive_halo[:,1],marker='.',c='b',s=500*mc)
#    ax.plot(positions_of_chosen_halos[:,0],positions_of_chosen_halos[:,1],'bo',markersize=10)


    velocities_of_chosen_halos = sp.array([chosen_halo.velocity\
            for chosen_halo in chosen_halos])
    velocities_of_halos_around_massive_halo = sp.array([halo_around_massive_halo.velocity\
            for halo_around_massive_halo in halos_around_massive_halo])
    velocity_corrections_for_chosen_halos = sp.array([chosen_halo.velocity_correction\
            for chosen_halo in chosen_halos])*linear_halo_bias


#    ax.quiver(positions_of_halos_around_massive_halo[:,0],
#              positions_of_halos_around_massive_halo[:,1],
#              velocities_of_halos_around_massive_halo[:,0],
#              velocities_of_halos_around_massive_halo[:,1],color='r',scale=scale)
    ax.quiver(positions_of_chosen_halos[:,0],
              positions_of_chosen_halos[:,1],
              velocities_of_chosen_halos[:,0],
              velocities_of_chosen_halos[:,1],color='b',scale=scale)
              
    ax.quiver(positions_of_chosen_halos[:,0],
              positions_of_chosen_halos[:,1],
              velocity_corrections_for_chosen_halos[:,0],
              velocity_corrections_for_chosen_halos[:,1],color='k',scale=scale)





    plt.savefig("/home/io/Dropbox/SharedStuff/Cepheids/Observation.pdf")

    
    print "ready to plot!"
    return 0



#        ax.quiver(position_CoM[0],position_CoM[1],velocity_CoM_nc[0],velocity_CoM_nc[1],color='r',scale_units='inches',scale=scale)
#        ax.quiver(position_CoM[0],position_CoM[1],velocity_correction[0],velocity_correction[1],color='black',scale_units='inches',scale=scale)
#        ax.quiver(position_local_CoM[0],position_local_CoM[1],local_velocity[0],local_velocity[1],color='g',scale_units='inches',scale=scale)
#        ax.quiver(position_local_CoM[0],position_local_CoM[1],local_velocity_correction[0],local_velocity_correction[1],color='grey',scale_unit
#        plt.savefig("/home/io/Desktop/cepheids_LCDM_%s.png" % total_mass)

    
