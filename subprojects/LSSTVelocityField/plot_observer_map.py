from __future__ import division
import scipy as sp
import matplotlib.pyplot as mplt
from matplotlib.patches import Polygon
import pdb
import cPickle
sys.path.insert(0, '/home/io/Dropbox/Projekter/Hubble/VelocityField/code')
import hubble_classes as hc
import healpy as hp
sys.path.insert(0, '/home/io/Dropbox/SharedStuff/code/colormap')
from option_d import get_colormap
from gplot import Plot
plt = Plot('io_latex_full','full')

def plot_powerspectrum(ls,cls,line):
    plt.plot(ls, sp.sqrt(ls*(ls+1) * cls),line)
    plt.xlabel('l') 
    plt.ylabel('sqrt(l(l+1)cl) [km/s]') 
    plt.grid()
    return 0

path = '/home/io/Dropbox/Projekter/Hubble/VelocityField/cases/Planck512/'
name = '_random_halos'
parameters = sp.load(path+'parameters'+name+'.save')
#parameters = cPickle.load(path+'parameters.save')
observers = sp.load(path+'observers'+name+'.npy')
#
#if skyfraction == '1.0':
#    title = "Full sky"
#else:
#    title = "Cut sky"
observer = observers[0]
binnumber=3
bindistance = parameters.bindistances[binnumber]
title = "Bin distance = %.0f" % bindistance    
cmap = get_colormap()
cmap.set_under("w")
hp.mollview(observer.vrmap[binnumber],min=-800,max=800,cmap=cmap, title=title)
plt.savefig('vrmap'+name+'_'+"%.0f" %bindistance +'.pdf')

plt.figure()
ls = observer.ls[binnumber]
cls = observer.cls[binnumber]
plot_powerspectrum(ls,cls,'b')
plt.savefig('powerspectrum'+name+'_'+"%.0f" %bindistance +'.pdf')
