from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import pdb
import cPickle
import sys
sys.path.insert(0, '/home/io/Dropbox/Projekter/Hubble/VelocityField/code')
import hubble_classes as hc


def get_ls_and_Cls_from_observers(parameterfile,observerfile):
    print "This function has not been tested yet"
    f = open(parameterfile,'r')    
    parameters = cPickle.load(f)
    f.close()
#    parameters = sp.load(path + 'parameters.save')
    observers = sp.load(observerfile)



#    bin_number = 0
    bin_number = 2
#    bin_number = 11
    bn = bin_number
    for mind_bin, maxd_bin in zip(parameters.bindistances[bn:bn+1],parameters.bindistances[bn+1:bn+2]):
        d = sp.mean((mind_bin,maxd_bin))
        ls = observers[0].ls[bin_number]
        Cls = sp.array([observer.cls[bin_number] for observer in observers])
    #    cls = observers[0].cls[bin_number]
        bin_number = bin_number+1

    return ls, Cls
    
def scale(ls,Cls):
    return sp.array(sp.sqrt(Cls*ls*(ls+1))) 