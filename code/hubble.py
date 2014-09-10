from __future__ import division
import scipy as sp
import sys
sys.path.insert(0,'/home/io/Dropbox/Projekter/Hubble/VelocityField/SNdata')
import yaml
import SN_redshiftdistribution as SN
import hubble_functions as hf
import Plotting_functions as pf
import MeanDeviation_functions as md



# There is one argument, namely the parameterfile
if len(sys.argv) != 2:
    print "Wrong number of argument"
    
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
number_of_observers = 2
observer_list = observer_list[0:number_of_observers]

print "The number of observers used is", len(observer_list)

if calculate_std_of_deviation == 1:
    mean_deviation = md.calculate_mean_deviation_over_surveyrange(observer_list, halo_list, boxsize, number_of_cones, skyfraction)
    print "The mean deviation is", mean_deviation 
    



 
   



    
    
