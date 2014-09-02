# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 12:04:29 2014

@author: io
"""

from __future__ import division
import scipy as sp
import ios_constants as ic
import matplotlib as m
import matplotlib.pyplot as plt
import hubvar_functions as hf
import sys

print 'Creating contourplots'

# Loading data (into a dictonary)
inputdir = sys.argv[1]#'/home/io/Desktop/PHD/Hubble/output/skyfraction_Planck512_2cones/'
fil = sys.argv[2]#'/home/io/Desktop/PHD/Hubble/output/Planck512/Hubbleconstants.txt'
N=int(sys.argv[3])
output = sys.argv[4]
write_to_tabular = int(sys.argv[5])
data = {}


for f in range(1,10):
    fil = inputdir+'Hubbleconstants0.'+str(f)+'.txt'
    rmax, H, N=hf.load_Hubbleconstants(fil,1,N)
    data[f] = H


rmax, H, N=hf.load_Hubbleconstants(fil,1,N)
data[10] = H

# Calculating skyfractions
#patch = sp.zeros(10)
#frak = sp.linspace(0.1,1,10)
#for i in range(0,10):
#    q = frak[i]
#    patch[i] = 2*sp.pi*(1-sp.cos(q*sp.pi))/(4*sp.pi)*100
patch = sp.linspace(0.1,1,10)    

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

cdict = {
  'blue'  :  ( (0.0, 0, 0), (0.02,0.3, 0.3), (0.3,1, 1), (1.,1, 1)),
  'green':  ( (0.0, 0, 0.0), (0.3,0, 0), (0.7,1, 1), (1.,1, 1)),
  'red' :  ( (0.0, 0, 0), (0.7,0, 0), (1.,1, 1))
}

cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 100)


    
x,y = sp.meshgrid(patch,rmax)
plt.close('all')

L = sp.linspace(0.01,0.09,9,endpoint=True)
l = sp.log10(L)


level=sp.linspace(l[0],l[-1],100)
plt.contourf(x,y,sp.log10(z68),level,shading='flat',cmap=cm)
cbar = plt.colorbar()

cbar.set_ticks(l)
cbar.set_ticklabels(L)
plt.xlabel('Observed percentage of the sky',fontsize=15)
plt.ylabel('Distance [Mpc/h]',fontsize=15)
plt.title('Width of 68.3% confidence interval for $H_{loc}/H_0$',fontsize=15)

plt.savefig(output)

if write_to_tabular == 1:
    mu67_pc, sigma67_pc = hf.mu_and_sigma(rmax,67,data[1]);
    mu150_pc, sigma150_pc = hf.mu_and_sigma(rmax,150,data[1]);
    mu256_pc, sigma256_pc = hf.mu_and_sigma(rmax,256,data[1]);

    tabel = open('/home/io/Dropbox/PHD/Python/tabel.txt','a')
    print >> tabel, 'skyfraction', '&', mu67_pc, '&', mu150_pc, '&', mu256_pc, '&', sigma67_pc, '&', sigma150_pc, '&', sigma256_pc, '\\\\'
    tabel.close()




