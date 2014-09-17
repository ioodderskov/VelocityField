from __future__ import division
import scipy as sp
import sys
sys.path.insert(0,'../SNdata')
import yaml
import SN_redshiftdistribution as SN
import hubble_functions as hf
import MeanDeviation_functions as md
#from gplot import Plot 
#plt = Plot('latex_full_hubble')
#import matplotlib.pyplot as mplt
#mplt.rc('font',family = 'serif')


# There is one argument, namely the parameterfile
if len(sys.argv) != 2:
    print "Wrong number of arguments"
    
parameterfile = sys.argv[1]

# Loads parameters
with open(parameterfile, 'r') as f:
    param = yaml.load(f)


# The halo data is loaded. Lines starting with a # are ignored.
data = sp.loadtxt(param["halofile"])

# Sort the halo data from lightest to heaviest
mass_sorted_data = sp.array(sorted(data,key=lambda data: data[2]))

# The rest of the parameters are read from the parameterfile
###############################################################################
hubblefile = param["hubblefile"]

observer_choice = param["observer_choice"]
observerfile = param["observerfile"]
number_of_observers = int(param["number_of_observers"])
[host_min_m, host_max_m] = sp.double([param["host_min_m"], param["host_max_m"]])
[sub_min_m, sub_max_m] = sp.double([param["sub_min_m"], param["sub_max_m"]])

[mind, maxd, width] = sp.double([param["mind"], param["maxd"], param["width"]])
number_of_SNe = int(param["number_of_SNe"])
boxsize = sp.double(param["boxsize"])
number_of_cones = int(param["number_of_cones"])
skyfraction = sp.double(param["skyfraction"])

calculate_std_of_deviation = int(param["calculate_std_of_deviation"])
calculate_hubble_constants = int(param["calculate_hubble_constants"])
calculate_redshiftdistribution = int(param["calculate_redshiftdistribution"])
make_hubblediagram = int(param["make_hubblediagram"])

vary_number_of_SNe = int(param["vary_number_of_SNe"])
min_number_of_SNe = int(param["min_number_of_SNe"])
max_number_of_SNe = int(param["max_number_of_SNe"])
step_number_of_SNe = int(param["step_number_of_SNe"])
###############################################################################



# Read the halo data
halo_list = hf.read_halo_file(mass_sorted_data)
    
# The observer positions are read or found:
observer_list = hf.find_observers(observer_choice,number_of_observers,boxsize,observerfile,halo_list,sub_min_m,sub_max_m,host_min_m,host_max_m)
print "The total number of observers is", len(observer_list)

#Selecting some of the observers
observer_list = observer_list[0:number_of_observers]
print "The number of observers used is", len(observer_list)




if calculate_hubble_constants:
    print "Calculating Hubble constants for increasing bin distances for all observers"     
  
#    # Calculate the bin-distances
    bindistances = hf.calculate_bindistances(mind, maxd, width)
#    
#    # Calculate the Hubbleconstants for all observers and all distances (or number of SNe)
    hf.calculate_hubble_constants_for_all_observers(range(len(observer_list)),observer_list, halo_list, number_of_SNe, bindistances, boxsize, number_of_cones, skyfraction)
#    
#    # Print the results to a file
    hf.print_hubbleconstants(hubblefile, bindistances, observer_list)


if calculate_std_of_deviation:
    print "Calculating the mean of the standard deviation over the survey-range specified in SN_redshiftdistribution.py"

    [zmin, zmax, Nbins] = [0.023, 0.1, 20]
    

    if calculate_std_of_deviation == 1:
        Wz_norm, zbins, N_tot = SN.get_table_distribution(zmin, zmax, Nbins)
        
    elif calculate_std_of_deviation == 2:
        bindistances = hf.calculate_bindistances(mind, maxd, width)
        skip, radial_velocities = hf.calculate_hubble_constants_for_all_observers(range(len(observer_list)),observer_list, halo_list, number_of_SNe, bindistances, boxsize, number_of_cones, skyfraction)
        Wz_norm, zbins, N_tot = SN.get_mock_distribution(radial_velocities, zmin, zmax, Nbins)
 
    Wz = Wz_norm*N_tot    


    mean_deviation = md.calculate_mean_deviation_over_surveyrange(observer_list, Wz, zbins,  halo_list, boxsize, number_of_cones, skyfraction)
    print "The mean deviation is", mean_deviation, "%"
    


if calculate_redshiftdistribution:
    print "Calculating the redshift distribution and comparing it to the one specified in SN_redshiftdistribution.py"
 
    [zmin, zmax, Nbins] = [0.01, 0.1, 20]   
    obs = 0
    
    bindistances = hf.calculate_bindistances(mind, maxd, width)
    skip, radial_velocities = hf.calculate_hubble_constants_for_all_observers(obs,observer_list, halo_list, number_of_SNe, bindistances, boxsize, number_of_cones, skyfraction)

    SN.get_table_distribution(zmin, zmax, Nbins)
    SN.get_mock_distribution(radial_velocities, zmin, zmax, Nbins)


    
if make_hubblediagram:

    obs = 0
    print "Saving data for the Hubble diagram for observer number", obs

    bindistances = hf.calculate_bindistances(mind, maxd, width)
    radial_distances, radial_velocities = hf.calculate_hubble_constants_for_all_observers(obs,observer_list, halo_list, number_of_SNe, bindistances, boxsize, number_of_cones, skyfraction)
    sp.save('radial_distances.npy', radial_distances)
    sp.save('radial_velocities.npy', radial_velocities)
    

    


    















