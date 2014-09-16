# -*- coding: utf-8 -*-
"""
Created on Fri May  2 11:58:18 2014

@author: io
"""

from __future__ import division
import scipy as sp
import hubvar_functions as hf
import sys
from gplot import Plot


# Arguments:
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
print 'output-folder=', sys.argv[1]
print 'Plottype=',sys.argv[2]

output = sys.argv[1]+'/Planck512_lightcone.pdf'
Nobs = 589
Plottype = sys.argv[2]
plt = Plot('latex_full_hubble',Plottype)
plt.rc('font',family = 'serif')


skip = list([63,91,226,246,332,367,403,406,493,545,576])


H = sp.array([])
for i in range(0,600):
    if i in skip:
        continue;
    fil = '../../cases/Planck512_lightcone/Hubbleconstants.'+str(i)+'.txt'
    rmax, Hi=hf.load_Hubbleconstants(fil,1,1);
    Hi = sp.reshape(Hi,len(rmax))
    H = sp.hstack((H,Hi))

number_of_bins = len(rmax)
H = sp.reshape(H,(Nobs,number_of_bins))

# Calculate the confidence limits, each containing the lower and upper bound  for every distance. 
sigma68, sigma95, sigma99 = hf.calc_confidence_intervals(H);


# The result is compared with the results from the standard analysis

comparisonfil = '../../cases/Planck512/Hubbleconstants.txt' 
rmax_sml, H_sml= hf.load_Hubbleconstants(comparisonfil,1,Nobs)
sigma68_sml, sigma95_sml, sigma99_sml = hf.calc_confidence_intervals(H_sml)

# Things are plotted and saved
hf.plot_patch(rmax,H,sigma68,sigma95,sigma99)


# The axes are formatted
plt.hold(True)
#
#plt.rcParams['axes.labelsize'] = 20
#

plt.xlabel('$r_{max}$ [Mpc/h]')
plt.axis([67.3, 256, 0.90, 1.10])
plt.xticks(sp.arange(80, 260, 20))

plt.ylabel('$H_{loc}/H_0$');
plt.yticks(sp.arange(0.90, 1.10+0.001, 0.025))

hf.plot_lines(rmax_sml,H_sml,sigma68_sml,sigma95_sml,sigma99_sml)

plt.change_size(152.4,130)
plt.finalize()
plt.savefig(output)


mu67_pc, sigma67_pc = hf.mu_and_sigma(rmax,67,H)
mu150_pc, sigma150_pc = hf.mu_and_sigma(rmax,150,H);
mu256_pc, sigma256_pc = hf.mu_and_sigma(rmax,256,H);

tabel = open('tabel.txt','a')
print >> tabel, 'Lightcone', '&', mu67_pc, '&', mu150_pc, '&', mu256_pc, '&', sigma67_pc, '&', sigma150_pc, '&', sigma256_pc, '\\\\'
tabel.close()


