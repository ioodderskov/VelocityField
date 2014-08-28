from __future__ import division
import scipy as sp
import sys
import yaml
import hubble_functions as hf
import hubble_plots as hp


# Chosen plot options
from gplot import Plot 
plt = Plot('prezi_hubble')
plt.rc('font',family = 'serif')

# This function makes the plots a bit prettier
def prettify(ax):

    for side in ['left','right','top','bottom']:
        ax.spines[side].set_linewidth(3)

    ax.tick_params('both', length=10, width=1, which='major')
    ax.tick_params('both', length=5, width=0.5, which='minor')





# There is one argument, namely the parameterfile
if len(sys.argv) != 2:
    print "Wrong number of argument"
    
parameterfile = sys.argv[1]

# Loads parameters
with open(parameterfile, 'r') as f:
    param = yaml.load(f)


# The halo data is loaded. Lines starting with a # are ignored.
data = sp.loadtxt(param["halofile"])
n_halos = data.shape[0]

[host_min_m, host_max_m] = sp.double([param["host_min_m"], param["host_max_m"]])
[sub_min_m, sub_max_m] = sp.double([param["sub_min_m"], param["sub_max_m"]])
[mind, maxd, width] = sp.double([param["mind"], param["maxd"], param["width"]])
observed_halos = int(param["observed_halos"])
boxsize = sp.double(param["boxsize"])
number_of_cones = int(param["number_of_cones"])
skyfraction = sp.double(param["skyfraction"])

# Making a list for the halos. Should I allocate memory more accurately?
halo_list = [None]*n_halos

for i in range(n_halos):
    [x,y,z] = data[i,[8,9,10]]
    [vx,vy,vz] = data[i,[11,12,13]]
    mass = data[i,2]
    ID = int(data[i,0])
    ID_host = int(data[i,33])
    halo = hf.Halos(x,y,z,vx,vy,vz,mass,ID,ID_host)
    halo_list[i] = halo
    
masses = [halo.mass for halo in halo_list]


# Identifying halos with masses corresponding to the Virgo Supercluster or the Local Group
localgroup_indices = (sub_min_m < masses) & (masses < sub_max_m)
virgo_indices = (host_min_m < masses) & (masses < host_max_m)
#
ID_hosts = sp.array([halo.ID_host for halo in halo_list])
IDs = sp.array([halo.ID for halo in halo_list])
virgo_IDs = IDs[virgo_indices]

#These methods only count the number of hosts, not the number of localgroup like subhalos
#observers = set(ID_hosts[localgroup_indices]) & set(virgo_IDs)
#observers = set(virgo_IDs).intersection(ID_hosts[localgroup_indices])

observers = [i for i, elem in enumerate(ID_hosts[localgroup_indices]) if elem in virgo_IDs]
all_indices = sp.array(range(len(halo_list)))
observer_indices = sp.array(all_indices[localgroup_indices])[observers]
print "The number of observers is", len(observers)
#observer_indices = sp.array([0,1])


observer_list = [None]*len(observer_indices)
bindistances = hf.calculate_bindistances(mind, maxd, width)
maxd = bindistances[-1]
#
##for halo_index,observer_number in zip(observers,range(len(observers_indices))):
for halo_index,observer_number in zip(sp.array(observer_indices)[[0,1]],range(1)):
#for halo_index,observer_number in zip(observer_indices,range(len(observer_indices))):
    observer = halo_list[halo_index]
    [x,y,z] = [observer.x, observer.y, observer.z]
    Hubbleconstants, xo, yo, zo, radial_distances, radial_velocities, selected_halos  = hf.find_hubble_constants_for_observer(halo_index,x,y,z,halo_list,mind,maxd,observed_halos,bindistances,boxsize,number_of_cones,skyfraction)
    observer = hf.Observers(x,y,z,selected_halos,Hubbleconstants)
    #Hastighederne er forkerte lige nu, men skal alligevel ikke bruges til noget!
    observer_list[observer_number] = observer

radial_distances = sp.array(radial_distances)
radial_velocities = sp.array(radial_velocities)
zo = sp.array(zo)
theta = sp.arccos(zo/radial_distances)    

ax = plt.subplot(1,1,1)
ax.plot(radial_distances, radial_velocities,'x')
plt.xlabel('$r [Mpc/h]$')
plt.ylabel('$v_r [km/s]$')

prettify(ax)
plt.finalize(custom = True)

output = '/home/io/Dropbox/PHD/DelA_rapport/Images/Hubbleplot_observer0.pdf'
plt.savefig(output)

#hp.make_3D_plot(xo,yo,zo)



