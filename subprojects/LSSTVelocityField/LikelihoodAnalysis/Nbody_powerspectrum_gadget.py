# -*- coding: utf-8 -*-
from __future__ import division
import scipy as sp
import matplotlib.pyplot as mplt
from gplot import Plot
plt = Plot('io_latex_full','full')
import sys
sys.path.insert(0, '/home/io/Dropbox/Projekter/Hubble/VelocityField/code')
import hubble_classes as hc
import hubble_functions as hf
import powerspectrum_functions as ps
import astropysics.constants as ac # link: https://pythonhosted.org/Astropysics/coremods/constants.html
# You can do, for example: ac.WMAP7Cosmology.omegaC

matterpowerspectrum_file = 'power512_011'

h = 0.678
box = 512
res = 1

data = sp.loadtxt(matterpowerspectrum_file)
ks = data[:,1]
power1 = data[:,2]
power2 = data[:,3]
#power = power1*box
#power = power1/power2**(1/3)*(2*sp.pi)**3*sp.pi**2
power = power2
plt.plot(ks,power)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('k [h/Mpc]')
plt.ylabel('?')
plt.plot([h/res,h/res],[1e-1,1e3],label='Nbody resolution')
plt.plot([h/box,h/box],[1e-1,1e3],label='Box')
plt.legend(frameon=0,loc='upper center')
mplt.show()