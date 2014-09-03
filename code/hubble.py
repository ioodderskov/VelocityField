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
mass_sorted_data = sp.array(sorted(data,key=lambda data: data[2]))
n_halos = len(mass_sorted_data)

observer_choice = param["observer_choice"]
hubblefile = param["hubblefile"]
observerfile = param["observerfile"]
number_of_observers = int(param["number_of_observers"])
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
    [x,y,z] = mass_sorted_data[i,[8,9,10]]
    [vx,vy,vz] = mass_sorted_data[i,[11,12,13]]
    mass = mass_sorted_data[i,2]
    ID = int(mass_sorted_data[i,0])
    ID_host = int(mass_sorted_data[i,33])
    halo = hf.Halos(x,y,z,vx,vy,vz,mass,ID,ID_host)
    halo_list[i] = halo
    
masses = [halo.mass for halo in halo_list]

# The observer positions are read or found:
observer_list = hf.find_observers(observer_choice,number_of_observers,boxsize,observerfile,halo_list,masses,sub_min_m,sub_max_m,host_min_m,host_max_m)


print "The number of observers is", len(observer_list)
#observer_indices = sp.array([0,1])


#observer_list = [None]*len(observer_indices)
bindistances = hf.calculate_bindistances(mind, maxd, width)
maxd = bindistances[-1]
#
for observer_number in range(len(observer_list)):
    observer = observer_list[observer_number]
    [x,y,z] = [observer.x, observer.y, observer.z]
    Hubbleconstants, radial_distances, radial_velocities, selected_halos  = hf.find_hubble_constants_for_observer(x,y,z,halo_list,mind,maxd,observed_halos,bindistances,boxsize,number_of_cones,skyfraction)

    observer.selected_halos = selected_halos
    observer.Hubbleconstants = Hubbleconstants

hf.print_hubbleconstants(hubblefile,bindistances,observer_list)


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
    
    
