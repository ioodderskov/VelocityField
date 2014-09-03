from __future__ import division
import scipy as sp
import sys
import yaml
import hubble_functions as hf


## Chosen plot options
#from gplot import Plot 
#plt = Plot('prezi_hubble')
#plt.rc('font',family = 'serif')
#
## This function makes the plots a bit prettier
#def prettify(ax):
#
#    for side in ['left','right','top','bottom']:
#        ax.spines[side].set_linewidth(3)
#
#    ax.tick_params('both', length=10, width=1, which='major')
#    ax.tick_params('both', length=5, width=0.5, which='minor')



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

vary_number_of_SNe = int(param["vary_number_of_SNe"])
min_number_of_SNe = int(param["min_number_of_SNe"])
max_number_of_SNe = int(param["max_number_of_SNe"])
step_number_of_SNe = int(param["step_number_of_SNe"])
###############################################################################




# Read the halo data
halo_list = hf.read_halo_file(mass_sorted_data)
    
# The observer positions are read or found:
observer_list = hf.find_observers(observer_choice,number_of_observers,boxsize,observerfile,halo_list,sub_min_m,sub_max_m,host_min_m,host_max_m)
print "The number of observers is", len(observer_list)




if vary_number_of_SNe == 0:

    # Calculate the bin-distances
    bindistances = hf.calculate_bindistances(mind, maxd, width)
    
    # Calculate the Hubbleconstants for all observers and all distances (or number of SNe)
    hf.calculate_hubble_constants_for_all_observers(observer_list, halo_list, number_of_SNe, bindistances, boxsize, number_of_cones, skyfraction)
    
    # Print the results to a file
    hf.print_hubbleconstants(hubblefile,bindistances,observer_list)







# This is actually a terrible way of doing this! I really should make another version
# of the calculate_hubbleconstants function, adding more and more halos for each observer.
if vary_number_of_SNe == 1:
    width = maxd-mind
    
    # Calculate the bin-distances
    bindistances = hf.calculate_bindistances(mind, maxd, width)
    
    # Calculate the Hubbleconstants for all observers and all distances (or number of SNe)
    for number_of_SNe in range(min_number_of_SNe,max_number_of_SNe+1,step_number_of_SNe):
        hf.calculate_hubble_constants_for_all_observers(observer_list, halo_list, number_of_SNe, bindistances, boxsize, number_of_cones, skyfraction)
        
        for observer in observer_list:
            print observer.Hubbleconstants
    #    # Print the results to a file
    #    hf.print_hubbleconstants_varSN(hubblefile,number_of_SNe,observer_list)







#radial_distances = sp.array(radial_distances)
#radial_velocities = sp.array(radial_velocities)
#zo = sp.array(zo)
#theta = sp.arccos(zo/radial_distances)    
#
#ax = plt.subplot(1,1,1)
#ax.plot(radial_distances, radial_velocities,'x')
#plt.xlabel('$r [Mpc/h]$')
#plt.ylabel('$v_r [km/s]$')
#
#prettify(ax)
#plt.finalize(custom = True)
#
#output = '/home/io/Dropbox/PHD/DelA_rapport/Images/Hubbleplot_observer0.pdf'
#plt.savefig(output)

#hp.make_3D_plot(xo,yo,zo)
    
    
