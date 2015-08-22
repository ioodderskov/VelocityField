from __future__ import division
from subprocess import call
import scipy as sp
import matplotlib.pyplot as plt
import theoretical_powerspectrum as tp
import measured_powerspectrum as mp

plot_measured_powerspectrum = 0

camb_matterpowerspectrum = 'test/test_z0_matterpower.dat'
MPS_file = camb_matterpowerspectrum
do_nonlinear = 0

h = 0.678
ombh2 = 0.022068
omch2 = 0.12029
scalar_amp = 2.21381e-9

omb = ombh2/h**2
omc = omch2/h**2

d_shell = 148/h # Mpc/h
omega_m = omb+omc

# Note from CAMB:

#"The params.ini file specifies the parameters used to run the program. 
#Comments in params.ini should make this fairly self-explanatory. 
#To produce the matter power spectrum in addition to CMB Cl set get_transfer = T; 
#the do_nonlinear input parameter determines whether this is output as the linear power spectrum 
#or includes non-linear corrections from the Halofit model. "


# Edit the input parameters in camb

# Read in the old parameter file:
with open('params.ini', 'r') as file:
    # read a list of lines into data
    data = file.readlines()

# Set the new matter density and amplitude of matter density fluctuations:
data[17] = 'do_nonlinear = %d\n' % do_nonlinear

data[35] = 'ombh2          = %0.6f\n' % (omb*h**2) 
data[36] = 'omch2          = %0.6f\n' % (omc*h**2)
data[86] = 'scalar_amp(1)             = %0.2e\n' % (scalar_amp)

# and write everything back
with open('params.ini', 'w') as file:
    file.writelines( data )


# Call CAMB to generate matterpowerspectrum
call(["./camb", "params.ini"])

# Load the matter power spectrum
data = sp.loadtxt("test/test_z0_matterpower.dat")
k_over_h = data[:,0]
power_in_Mpc3_over_h3 = data[:,1]
plt.figure()
plt.xscale('log')
plt.yscale('log')
plt.plot(k_over_h,power_in_Mpc3_over_h3)


plt.figure()
tp.velocitypowerspectrum(MPS_file,0,h,omega_m,d_shell)
#velocitypowerspectrum(0.03) #Mpc. Roughly the diameter of the Milky Way
#tp.velocitypowerspectrum(MPS_file,4/h,h,omega_m,d_shell) #Mpc. 8 neighbours (from healpy.get_all_neighbours) and eq. 8 in paper 2
#tp.velocitypowerspectrum(MPS_file,150/h*0.13,h,omega_m,d_shell) # Roughly the diameter of the largest hole

# Smoothing with the scale of the Gaussian beam in healpy
R_rad = 0.05
R_phys = R_rad*d_shell
tp.velocitypowerspectrum(MPS_file,R_phys,h,omega_m,d_shell)

# Smoothing on the scale of the pixel size
R_rad = 0.14
R_phys = R_rad*d_shell
tp.velocitypowerspectrum(MPS_file,R_phys,h,omega_m,d_shell)

plt.show()
if plot_measured_powerspectrum:
    case = 'Planck512'
    skyfraction = '1'
    path = '/home/io/Dropbox/Projekter/Hubble/VelocityField/cases/'+case+'/'
    parameterfile = path+'parameters.save'
    observerfile = path+'observers.npy'
    ls, Cls = mp.get_ls_and_Cls_from_observers(parameterfile,observerfile)
    
    #Only testing with one observer
    Cls = Cls[0]
    scaled_Cls = mp.scale(ls,Cls)
    plt.plot(ls,scaled_Cls,'--')
    plt.show()
    #plt.xscale('log')
    #plt.yscale('log')