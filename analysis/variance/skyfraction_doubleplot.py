# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 12:04:29 2014

@author: io
"""

from __future__ import division
import scipy as sp
import ios_constants as ic
import matplotlib as m
#import matplotlib.pyplot as plt
import hubvar_functions as hf
import sys
from gplot import Plot
import matplotlib.colorbar as mpl



def get_data(inputdir):

    # Loading data (into a dictonary)
    fil_fullsky = '/home/io/Desktop/PHD/Hubble/output/Planck512/Hubbleconstants.txt'
    N = 600
    data = {}
    
    for f in range(1,10):
        fil = inputdir+'Hubbleconstants0.'+str(f)+'.txt'
        rmax, H, N=hf.load_Hubbleconstants(fil,1,N)
        data[f] = H
    
    
    rmax, H, N=hf.load_Hubbleconstants(fil_fullsky,1,N)
    data[10] = H
    
    # Calculating skyfractions
    #patch = sp.zeros(10)
    #frak = sp.linspace(0.1,1,10)
    #for i in range(0,10):
    #    q = frak[i]
    #    patch[i] = 2*sp.pi*(1-sp.cos(q*sp.pi))/(4*sp.pi)*100

    
    number_of_bins = len(H[1,:]);
    z68 = sp.zeros((10*number_of_bins,1));
    
    count = 0     
    for i in range(1,11):
        H = data[i]
         
        for b in range(0,number_of_bins):
            sl68, su68 = hf.confidenslimit(H[:,b],0.683)
            z68[count] = sl68+su68
            count = count+1
        
            
            
    z68 = sp.reshape(z68,(10,number_of_bins))
    z68 = z68.T

    return(rmax,z68)

# Making plots
plt = Plot('latex_full_hubble','h2c')
plt.rc('font',family = 'serif')

inputdir = '/home/io/Desktop/PHD/Hubble/output'
#output='/home/io/Dropbox/PHD/Python/tests/Planck512_skyfraction.pdf'
output = '/home/io/Dropbox/hubble2013/Planck512_skyfraction.pdf'

cdict = {
  'blue'  :  ( (0.0, 0, 0), (0.02,0.3, 0.3), (0.3,1, 1), (1.,1, 1)),
  'green':  ( (0.0, 0, 0.0), (0.3,0, 0), (0.7,1, 1), (1.,1, 1)),
  'red' :  ( (0.0, 0, 0), (0.7,0, 0), (1.,1, 1))
}

cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 100)

patch = sp.linspace(0.1,1,10)    

input_1cone = inputdir + '/skyfraction_Planck512_1cone/'
input_2cones = inputdir + '/skyfraction_Planck512_2cones/'

rmax, z68_1cone = get_data(input_1cone)
rmax, z68_2cones = get_data(input_2cones)

    
x,y = sp.meshgrid(patch,rmax)
#plt.close('all')

#L = sp.linspace(0.01,0.09,9,endpoint=True)
#L = sp.linspace(0.004,0.09,9,endpoint=True)
L = sp.array([0.004,0.01,0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09])
l = sp.log10(L)


level=sp.linspace(l[0],l[-1],100)

ax1 = plt.subplot(1,3,1)
plt.contourf(x,y,sp.log10(z68_1cone),level,shading='flat',cmap=cm)
#plt.contourf(x,y,sp.log10(z68_1cone),level,shading='flat')
plt.xlabel('Fraction of sky observed')
plt.ylabel('Distance [Mpc/h]')

ax2 = plt.subplot(1,3,2)
#im = plt.contourf(x,y,sp.log10(z68_2cones),level,shading='flat')
im = plt.contourf(x,y,sp.log10(z68_2cones),level,shading='flat',cmap=cm)
plt.xlabel('Fraction of sky observed')
#plt.ylabel('Distance [Mpc/h]')


#plt.figure.subplots_adjust(right=0.8)
cax = plt.subplot(1,3,3)
#print  plt.axes
cbar = plt.colorbar(im, cax = cax)
#print  plt.axes
cbar.set_ticks(l)
cbar.set_ticklabels(L)
cbar.ax.tick_params(labelsize=7) 

plt.suptitle('Width of 68.3% confidence interval for $H_{loc}/H_0$',x = 0.425)


plt.change_size(152.4,75)
plt.finalize()

#
plt.savefig(output)
plt.show()
#print("Here!!!!")
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as mplt
#fig = mplt.figure()
##fig.clf()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(x,y,sp.log10(z68_2cones))
#mplt.show()


