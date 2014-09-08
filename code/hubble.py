from __future__ import division
import scipy as sp
import sys
import yaml
import hubble_functions as hf
import matplotlib.pyplot as plt

# Chosen plot options
#from gplot import Plot 
#plt = Plot('latex_full_hubble')
plt.rc('font',family = 'serif')

# This function makes the plots a bit prettier
def prettify(ax):

    for side in ['left','right','top','bottom']:
        ax.spines[side].set_linewidth(3)

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



# Calculate the bin-distances
bindistances = hf.calculate_bindistances(mind, maxd, width)

# Calculate the Hubbleconstants for all observers and all distances (or number of SNe)
radial_distances, radial_velocities = hf.calculate_hubble_constants_for_all_observers(observer_list, halo_list, number_of_SNe, bindistances, boxsize, number_of_cones, skyfraction)

# Print the results to a file
#hf.print_hubbleconstants(hubblefile,bindistances,observer_list)




radial_distances = sp.array(radial_distances)
radial_velocities = sp.array(radial_velocities)

z = radial_velocities/(3e5)
H = 100
c = 3e5
zhub = H/c*radial_distances

dL = radial_distances*(1+zhub)
Mv = -19.3 # (wikipedia)
m02 = 0.2*(5*sp.log10(dL)+Mv+25) # The apparent lumonisity is proportional to this quantity
log_cz = sp.log10(radial_velocities)

m02_min = 0.2*(5*sp.log10(mind*(1+H/c*mind))+Mv+25)
m02_max = 0.2*(5*sp.log10(maxd*(1+H/c*maxd))+Mv+25)

#zo = sp.array(zo)
#theta = sp.arccos(zo/radial_distances)    

ax = plt.subplot(1,1,1)
ax.plot(m02, log_cz,'ro')
#ax.plot([m02_min,m02_min], [3.3,4.7],'r')
#ax.plot([m02_max,m02_max], [3.3,4.7],'r')

#plt.xlabel('$r [Mpc/h]$')
#plt.ylabel('$v_r [km/s]$')


plt.xlabel('$0.2m_V$ (mag)')
plt.ylabel('$\log cz$')
plt.axis([2.6,4.0,3.3,4.7])


prettify(ax)
#plt.finalize(custom = True)

output = '/home/io/Dropbox/SharedStuff/hubble2013/redshiftdistribution0.pdf'
plt.savefig(output)

#hp.make_3D_plot(xo,yo,zo)
    
    
