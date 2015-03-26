from __future__ import division
import scipy as sp
import scipy.ndimage as nd

periodic_boundaries = 1

def positions_and_values_from_grid(parameters,values_tsc):

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

                cell_value = values_tsc[x_index,y_index,z_index]
                gridpoint_values.append(cell_value)
                gridpoint_positions.append(sp.array([x,y,z]))
                

    gridpoint_values = sp.array(gridpoint_values)
        
    return gridpoint_positions, gridpoint_values


def velocity_field(parameters,positions,velocities):
    
    # Need to convert from velocity-density to velocity:
    cellsize = parameters.boxsize/parameters.Ng
    
    vx_tsc = tsc(parameters,positions,velocities[:,0])*cellsize**3
    vy_tsc = tsc(parameters,positions,velocities[:,1])*cellsize**3
    vz_tsc = tsc(parameters,positions,velocities[:,2])*cellsize**3

    if parameters.smoothing:
        vx_tsc_smoothed = nd.gaussian_filter(vx_tsc,parameters.smoothing_radius,mode='wrap') 
        vy_tsc_smoothed = nd.gaussian_filter(vy_tsc,parameters.smoothing_radius,mode='wrap') 
        vz_tsc_smoothed = nd.gaussian_filter(vz_tsc,parameters.smoothing_radius,mode='wrap')         
        
        vx,vy,vz = vx_tsc_smoothed,vy_tsc_smoothed,vz_tsc_smoothed
    else:
        vx,vy,vz = vx_tsc,vy_tsc,vz_tsc
        
        
        
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
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (x,y,z,rho,vx,vy,vz))
                
    f.close()
    
def create_density_and_velocity_grid(parameters,halos):
    positions = sp.array([halo.position for halo in halos])
    masses = sp.array([halo.mass for halo in halos])
    velocities = sp.array([halo.velocity for halo in halos])

    rho_tsc = tsc(parameters,positions,masses)
    gridpoint_positions, gridpoint_rho = positions_and_values_from_grid(parameters,rho_tsc)

    gridpoint_positions, gridpoint_velocities = velocity_field(parameters,positions,velocities)

    write_grid_to_file(parameters,gridpoint_positions,gridpoint_rho,gridpoint_velocities)    

    