# -*- coding: utf-8 -*-
"""
Created on Sat Apr 19 14:21:19 2014

@author: io
"""

from __future__ import division
import scipy as sp
#import matplotlib.pyplot as plt
from gplot import Plot
import hubvar_functions as hf
import sys



# Arguments:
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
print 'runname=', sys.argv[1]
print 'a=', sys.argv[2]
print 'Nobservers=', sys.argv[3]
print 'x-axis=', sys.argv[4]
print 'output-folder=', sys.argv[5]
print 'Compare=', sys.argv[6]
print 'Write to tabular=', sys.argv[7]
print 'Plus Planck-band=', sys.argv[8]
print 'Plottype=',sys.argv[9]


runname = sys.argv[1]
a = float(sys.argv[2])
N_observers = int(sys.argv[3])
xaxis = sys.argv[4]
output = sys.argv[5]+ '/'+runname+'.pdf'
compare = int(sys.argv[6])
write_to_tabular = int(sys.argv[7])
planck_band = int(sys.argv[8])
Plottype = sys.argv[9]

plt = Plot('latex_full_hubble',Plottype)
plt.rc('font',family = 'serif')


# Inputfile:
fil = '../../cases/'+runname+'/Hubbleconstants.txt'

# Load distances and Hubbleconstants, and count the maximal number of observers
lightcone = 0;
rmax, H=hf.load_Hubbleconstants(fil,a,N_observers);
print "Number of observers = ", len(H)

# Calculate the confidence limits, each containing the lower and upper bound  for every distance. 
sigma68, sigma95, sigma99 = hf.calc_confidence_intervals(H);
error = hf.mean_sys_error(H,rmax)
print "The mean systematical error is", error*100, "%" 

# If activated, the result is compared with the results from the standard analysis
if compare == 1:
    comparisonfil = '../../cases/Planck512/Hubbleconstants.txt' 
    rmax_sml, H_sml= hf.load_Hubbleconstants(comparisonfil,1,N_observers)
    sigma68_sml, sigma95_sml, sigma99_sml = hf.calc_confidence_intervals(H_sml)

# Things are plotted and saved
hf.plot_patch(rmax,H,sigma68,sigma95,sigma99)


# The axes are formatted
#plt.hold(True)
#
#plt.rcParams['axes.labelsize'] = 20
#

xmin = 67.3
xmax = 256
#xmin = 4
#xmax = 8
if planck_band == 1:
    ymin = 0.85
    ymax = 1.15
else:
    ymin = 0.9
    ymax = 1.1
#    ymin = 0
#    ymax = 2
    
if xaxis == 'rmax':
    plt.xlabel('$r_{max}$ [Mpc/h]')
    plt.axis([xmin, xmax, ymin, ymax])
    plt.xticks(sp.arange(80, 260, 20))
if xaxis == 'numSN':
    plt.xlabel('Number of observed halos')
    plt.axis([2, 200, ymin, ymax])
    plt.xticks(sp.arange(20,200,20))

plt.ylabel('$H_{loc}/H_0$');
plt.yticks(sp.arange(ymin, ymax+0.001, 0.025))

if compare == 1:
        hf.plot_lines(rmax_sml,H_sml,sigma68_sml,sigma95_sml,sigma99_sml)
        
if planck_band == 1:
    HRiess = 73.8/67.8*sp.ones(2)
    sigma = 2.4/67.8
    plt.plot([xmin,xmax],HRiess,color='green')
    plt.fill_between([xmin,xmax],HRiess-sigma, HRiess+sigma,facecolor='green',alpha=0.3)

plt.finalize()
plt.change_size(152.4,130)
plt.savefig(output)


if write_to_tabular == 1:
    mu67_pc, sigma67_pc = hf.mu_and_sigma(rmax,67,H)
    mu150_pc, sigma150_pc = hf.mu_and_sigma(rmax,150,H);
    mu256_pc, sigma256_pc = hf.mu_and_sigma(rmax,256,H);

    tabel = open('tabel.txt','a')
    print >> tabel, runname, '&', mu67_pc, '&', mu150_pc, '&', mu256_pc, '&', sigma67_pc, '&', sigma150_pc, '&', sigma256_pc, '\\\\'
    tabel.close()



#a = gca()
#a.set_xticklabels(a.get_xticks(), fontProperties)
#a.set_yticklabels(a.get_yticks(), fontProperties)

#

#
#% 
#%  HRiess = 73.8/67.8*ones(size(rmax));
#% sigmaRiess = [HRiess+2.4/67.8,fliplr(HRiess-2.4/67.8)];
#% Riess = patch([rmax fliplr(rmax)],sigmaRiess,'k');
#% plot(rmax,HRiess,'Color','k','LineWidth',1.3);
#%  set(Riess,'EdgeColor','None');
#%  alpha(Riess,0.05)
#
#
# 
# axis([64 256 0.90 1.10])
#%   axis([64 256 0.90 1.10])
# %axis([64 256 0 50])
#%  axis([2 200 0.90 1.10])
#%  set(gca,'YTick',0.90:0.025:1.10)
#  set(gca,'YTick',0.90:0.025:1.10)
# 
# lx = xlabel('$r_{max} [Mpc/h]$','FontSize',14);
#%  lx = xlabel('Number of observed supernovae','FontSize',14);
# ly = ylabel('$H_{loc}/H_0$','FontSize',14);
# %ly = ylabel('$\Bigl(\sum_i^N(\vec(H)_i\Bigr)^2/N$');
# set(lx,'Interpreter','Latex');
# set(ly,'Interpreter','Latex');
# 
# box on;
# set(gca,'Layer','top')
# size(mean(H))
# 
# % Kvadratrodsfit til antallet af supernovaer
#
#%  su68 = sigma68(:,1)'+mean(H);
#%  size(su68)
#%  size(rmax)
#
# 
# % Hvis varSN
#  %NSN = rmax;
# 
#%  % Hvis alle de andre
#%   Nrmax = 1:length(rmax);
#%   NSN = Nrmax*10;
#% 
#%  fejl =  @(p) sum((p(1)./sqrt(NSN)+p(2)-su68).^2);
#%  pmin = fminsearch(fejl,[1,1]);
#%  plot(rmax,pmin(1)./sqrt(NSN)+pmin(2))
#%  fejl(pmin)
# 
#%  print(gcf,'-depsc2',['~/Dropbox/hubble2013/',runname,'.eps']) 
# % ~/Dropbox/hubble2013/
# hold off 
