from __future__ import division
import scipy as sp
#import scipy.ndimage as nd
from copy import copy
from astropy.convolution import convolve 
from astropy_make_kernel import make_kernel

periodic_boundaries = 1

def positions_and_values_from_grid(parameters,values):

    # I write the densities to a vector (and construct another vector with the corresponding coordinates)
    # While this is done, I interpolate to a finer grid        

    # When outputting the array, it is structured like this
    # The layers corresponds to increasing values of y.
    # -------> z 
    # |
    # |
    # |
    # v
    # x
    
    
    cellsize = parameters.boxsize/parameters.Ng
    grid_x = grid_y = grid_z = sp.array(range(parameters.Ng))*cellsize+cellsize/2
    
    gridpoint_values = []
    gridpoint_positions = []

    
    for x_index,x in enumerate(grid_x):
        for y_index,y in enumerate(grid_y):
            for z_index,z in enumerate(grid_z):

                cell_value = values[x_index,y_index,z_index]
                gridpoint_values.append(cell_value)
                gridpoint_positions.append(sp.array([x,y,z]))
                

    gridpoint_values = sp.array(gridpoint_values)
        
    return gridpoint_positions, gridpoint_values


def ngp(parameters,positions,values):
    
    values_ngp = sp.zeros((parameters.Ng,parameters.Ng,parameters.Ng))
    counts_ngp = sp.zeros((parameters.Ng,parameters.Ng,parameters.Ng))
    cellsize = parameters.boxsize/parameters.Ng


    for position,pvalue in zip(positions,values):

        position = sp.array(position)
        
        position_cellunits = position/cellsize

        # cell indices
        cell_indices = sp.floor(position_cellunits)
        

        if periodic_boundaries:
            cell_indices = sp.mod(cell_indices,parameters.Ng)

        index_x, index_y, index_z = cell_indices[0],cell_indices[1],cell_indices[2]


        values_ngp[index_x][index_y][index_z] += pvalue
        counts_ngp[index_x][index_y][index_z] += 1                                    

    values_ngp = sp.array(values_ngp)/sp.array(counts_ngp)
    print "Don't mind this warning. Astropy can handle nan-values"                
    

    return values_ngp     
    
    
def velocity_field(parameters,positions,velocities):
    


    vx_ngp = ngp(parameters,positions,velocities[:,0])
    vy_ngp = ngp(parameters,positions,velocities[:,1])
    vz_ngp = ngp(parameters,positions,velocities[:,2])



    if parameters.smoothing:  
        cellsize = parameters.boxsize/parameters.Ng
        smoothing_radius_cells = parameters.smoothing_radius/cellsize
        kernel_res = int(smoothing_radius_cells*2*3)
        if sp.mod(kernel_res,2) == 0:
            kernel_res = kernel_res+1
        kernel_res = 3
        print "the kernel resolution is", kernel_res
        kernel = make_kernel([kernel_res,kernel_res,kernel_res],smoothing_radius_cells)

        vx_ngp_smoothed = copy(vx_ngp)
        vy_ngp_smoothed = copy(vy_ngp)
        vz_ngp_smoothed = copy(vz_ngp)

        n = 0
        while sp.any(sp.isnan(vx_ngp_smoothed)):
            n = n+1
            vx_ngp_smoothed = convolve(vx_ngp_smoothed,kernel,boundary='wrap') 
            vy_ngp_smoothed = convolve(vy_ngp_smoothed,kernel,boundary='wrap') 
            vz_ngp_smoothed = convolve(vz_ngp_smoothed,kernel,boundary='wrap') 
        print "Went through the convolving-loop", n, "times"


        vx,vy,vz = vx_ngp_smoothed,vy_ngp_smoothed,vz_ngp_smoothed

    else:
        vx,vy,vz = vx_ngp,vy_ngp,vz_ngp

    gridpoint_positions, gridpoint_vx = positions_and_values_from_grid(parameters,vx)
    gridpoint_positions, gridpoint_vy = positions_and_values_from_grid(parameters,vy)
    gridpoint_positions, gridpoint_vz = positions_and_values_from_grid(parameters,vz)
    

        

    gridpoint_velocities = sp.array([gridpoint_vx,gridpoint_vy,gridpoint_vz]).T
    
    
    return gridpoint_positions, gridpoint_velocities    
    
    
    
    



def tsc(parameters,positions,values):
    
    values_tsc = sp.zeros((parameters.Ng,parameters.Ng,parameters.Ng))
    cellsize = parameters.boxsize/parameters.Ng


    for position,pvalue in zip(positions,values):
        
        position = sp.array(position)
        
        position_cellunits = position/cellsize



        cell_indices = sp.floor(position_cellunits)
        leftcell_indices = cell_indices - 1
        rightcell_indices = cell_indices + 1


        cell_position = cell_indices + 0.5
        leftcell_position = leftcell_indices + 0.5
        rightcell_position = rightcell_indices + 0.5
        
        particle_cell_distances = sp.absolute(position_cellunits - cell_position)
        particle_leftcell_distances = sp.absolute(position_cellunits - leftcell_position)        
        particle_rightcell_distances = sp.absolute(position_cellunits - rightcell_position)

        weights_cell = 0.75 - sp.square(particle_cell_distances)
        weights_leftcell = 0.5*sp.square(1.5 - particle_leftcell_distances)
        weights_rightcell = 0.5*sp.square(1.5 - particle_rightcell_distances)



        if periodic_boundaries:
            cell_indices = sp.mod(cell_indices,parameters.Ng)
            leftcell_indices = sp.mod(leftcell_indices,parameters.Ng)
            rightcell_indices = sp.mod(rightcell_indices,parameters.Ng)


        indices_x, weights_x = [cell_indices[0],leftcell_indices[0],rightcell_indices[0]],\
                                [weights_cell[0],weights_leftcell[0],weights_rightcell[0]]
        indices_y, weights_y = [cell_indices[1],leftcell_indices[1],rightcell_indices[1]],\
                                [weights_cell[1],weights_leftcell[1],weights_rightcell[1]]
        indices_z, weights_z = [cell_indices[2],leftcell_indices[2],rightcell_indices[2]],\
                                [weights_cell[2],weights_leftcell[2],weights_rightcell[2]]

        for index_x,weight_x in zip(indices_x, weights_x):
            for index_y,weight_y in zip(indices_y, weights_y):
                for index_z,weight_z in zip(indices_z, weights_z):
                    values_tsc[index_x][index_y][index_z] += pvalue*weight_x*weight_y*weight_z/cellsize**3                                    


                
    return values_tsc 
    
def write_grid_to_file(parameters,gridpoint_positions,gridpoint_rho,gridpoint_velocities):

    f = open(parameters.gridfile,'w')
    
    f.write("#x \t y \t z\t rho \t vx \t vy \t vz\n")
    for gridpoint_position,rho,gridpoint_velocity in zip(gridpoint_positions,gridpoint_rho,gridpoint_velocities):
        x,y,z = gridpoint_position[0],gridpoint_position[1],gridpoint_position[2]
        vx,vy,vz = gridpoint_velocity[0],gridpoint_velocity[1],gridpoint_velocity[2]        
        f.write("%.2f\t%.2f\t%.2f\t%.2e\t%.2f\t%.2f\t%.2f\n" % (x,y,z,rho,vx,vy,vz))
                
    f.close()
    
def create_density_and_velocity_grid(parameters):
    positions = sp.array([halo.position for halo in parameters.halos])
    masses = sp.array([halo.mass for halo in parameters.halos])


    rho_tsc = tsc(parameters,positions,masses)
    gridpoint_positions, gridpoint_rho = positions_and_values_from_grid(parameters,rho_tsc)

    if parameters.velocities_on_grid:
        velocities = sp.array([halo.velocity for halo in parameters.halos])
        gridpoint_positions, gridpoint_velocities = velocity_field(parameters,positions,velocities)
    else:
        gridpoint_velocities = sp.zeros_like(gridpoint_positions)

    write_grid_to_file(parameters,gridpoint_positions,gridpoint_rho,gridpoint_velocities)    

    
    