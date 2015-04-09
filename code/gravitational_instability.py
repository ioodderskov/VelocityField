from __future__ import division
import numpy as sp
import scipy.linalg as linalg
import pdb


periodic_boundaries = 1

h = 0.702
meter = 1. # meter
second = 1. # second
kg = 1. # kg
Mpc = 1e6*3.085678e16*meter
Msun = 1.989e30*kg
H0 = 100. *1e3*meter/second/(Mpc/h) # s**-1
Gc = 6.67e-11 *meter**3/(kg*second**2)
rho_critical = 3.*H0**2/(8*sp.pi*Gc) # kg/m**3
rho_critical_astro = rho_critical / (Msun/h) * (Mpc/h)**3


def periodic_coordinate(coordinate,boxsize):
    
    if coordinate > boxsize/2:
        coordinate = coordinate-boxsize
    if coordinate < -boxsize/2:
        coordinate = coordinate+boxsize
        
    return coordinate


def velocity_from_matterdistribution(parameters,observer_position,halo_position,halos):
    
    positions = [halo.position for halo in halos]
    masses = [halo.mass for halo in halos]

    
#    rho_mean = sp.sum(masses)/parameters.boxsize**3
    rho_mean = rho_critical_astro*parameters.omegam

#    pdb.set_trace()

    #PERIODISKE GRAENSER:
    if periodic_boundaries:
        positions = positions-observer_position
    
        for index,position in enumerate(positions):
            x,y,z = periodic_coordinate(position[0],parameters.boxsize),periodic_coordinate(position[1],parameters.boxsize),periodic_coordinate(position[2],parameters.boxsize)
            positions[index] = [x,y,z]
            
        positions = sp.array(positions)+observer_position
    
    #SURVEY CENTERED ON HALO
#    individual_contributions = sp.array([masses[i]/rho_mean \
#                                        *(positions[i]-halo_position) \
#                                        /linalg.norm(positions[i]-halo_position)**3 \
#                                        for i in range(len(masses)) \
#                                        if (linalg.norm(positions[i]-halo_position) > parameters.min_dist) \
#                                        & (linalg.norm(positions[i]-halo_position) < parameters.survey_radius) ])

     #SURVEY CENTERED ON OBSERVER
    individual_contributions = sp.array([masses[i]/rho_mean \
                                        *(positions[i]-halo_position) \
                                        /linalg.norm(positions[i]-halo_position)**3 \
                                        for i in range(len(masses)) \
                                        if (linalg.norm(positions[i]-halo_position) > parameters.min_dist) \
                                        & (linalg.norm(positions[i]-observer_position) < parameters.survey_radius) ])
                
#    pdb.set_trace()
    
    if len(individual_contributions) == 0:
        print "No halos was found between min_dist and max_dist"
        v_from_matterdistribution = sp.array([0,0,0])

    else:

        H0_astro = 100.
        f = parameters.omegam**0.6
    
        vx_from_matterdistribution = f*H0_astro/(4*sp.pi)*sp.sum(individual_contributions[:,0])
        vy_from_matterdistribution = f*H0_astro/(4*sp.pi)*sp.sum(individual_contributions[:,1])
        vz_from_matterdistribution = f*H0_astro/(4*sp.pi)*sp.sum(individual_contributions[:,2])
        
        v_from_matterdistribution = sp.array([vx_from_matterdistribution,vy_from_matterdistribution,vz_from_matterdistribution])
    
#        pdb.set_trace()
    return v_from_matterdistribution


    