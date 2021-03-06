from __future__ import division
import healpy as hp
import scipy as sp
import pdb
import copy
from collections import Counter
from sympy.physics.quantum.cg import Wigner3j
from scipy import linalg


def create_map(parameters,thetas,phis,vrs):
    
    
    pix = hp.ang2pix(parameters.nside,thetas,phis)
    number_of_pixels = hp.nside2npix(parameters.nside)
    vrmap = sp.ones(number_of_pixels)*parameters.badval
    
#    pdb.set_trace()
    vrs = sp.array(vrs)
    vrs_mean_of_repeated_pixels = copy.copy(vrs)
    for p in set(pix):
        vrs_mean_of_repeated_pixels[pix == p] = sp.mean(vrs[pix == p])
    
    vrmap[pix] = vrs_mean_of_repeated_pixels

#    pdb.set_trace()    
    theta_max = sp.arccos(1-2*parameters.skyfraction)
    pix_all = sp.array(range(number_of_pixels))
    pix_unseen = pix_all[hp.pix2ang(parameters.nside,pix_all)[0]>theta_max]
    vrmap[pix_unseen] = parameters.unseen

    empty_pixels = number_of_pixels-len(set(pix))
    print "The number of empty pixels inside the survey is", empty_pixels
    print "The corresponds to a fraction of", empty_pixels/number_of_pixels
    print "The number of pixels outside survey is", len(pix_unseen)
#    pdb.set_trace()
    return vrmap
    
    
    
    
def pixelsize_in_radians(parameters):
#    example = hp.rotator.angdist(hp.pix2ang(parameters.nside,0),
#                                 hp.pix2ang(parameters.nside,1))
    return 2*sp.arccos(1-4*sp.pi/hp.nside2npix(parameters.nside)/(2*sp.pi))

def find_largest_hole(parameters,ar):
    
    minimal_distances = []
    all_pixels = sp.array(range(len(ar)))
    empty_pixels = all_pixels[(ar[all_pixels] == parameters.badval)]
    if len(empty_pixels) == 0:
        print "no empty pixels"
        return 2*sp.pi/(6*parameters.nside)
        
        
    print "the number of empty pixels is", len(empty_pixels)
    nonempty_pixels = all_pixels[(ar[all_pixels] != parameters.badval)\
                                & (ar[all_pixels] != parameters.unseen)]
        
    for p in empty_pixels:

        minimal_distance = 3.14
        theta,phi = hp.pix2ang(parameters.nside,p)
        
        for p_i in nonempty_pixels:
            
            theta_i, phi_i = hp.pix2ang(parameters.nside,p_i)
            angular_distance = hp.rotator.angdist([theta,phi],[theta_i,phi_i])
            minimal_distance = sp.minimum(minimal_distance,angular_distance)
            
        minimal_distances.append(minimal_distance)
        
         
    radius_of_largest_hole = sp.amax(minimal_distances)

    print "The angular radius of largest hole = ", radius_of_largest_hole, "rad."
    
    return radius_of_largest_hole



def apply_mask(parameters,ar,mask_badval,mask_unseen):
    mask = sp.zeros_like(ar)
    if mask_badval:
        mask[ar == parameters.badval] = 1
    if mask_unseen:
        mask[ar == parameters.unseen] = 1

    ar_masked = hp.ma(ar)
    ar_masked.mask = mask
    
    return ar_masked, mask

def apply_window(parameters,ar,window_badval,window_unseen,smoothing_radius):

        
    window = sp.ones_like(ar)
    if window_badval:
        window[ar == parameters.badval] = 0
    if window_unseen:
        window[ar == parameters.unseen] = 0            
            
        
    window = hp.smoothing(window,fwhm=smoothing_radius,verbose=False)
    
    ar_windowed = window*ar
    
    
    return ar_windowed, window

def get_pseudo_powerspectrum(ar,lmax):
    ls = sp.arange(0,lmax+1)
    print "I am removing the mono- and dipole"
    ar = hp.pixelfunc.remove_dipole(ar)    
    pseudo_cls = hp.anafast(ar, lmax=lmax,use_weights=True)
    return ls,pseudo_cls

def get_MASTER_corrected_powerspectrum(pseudo_cls,window,lmax):
    ls = sp.arange(0,lmax+1)
    Wls = hp.anafast(window,lmax=lmax,use_weights=True)
    M = sp.matrix(sp.zeros((lmax+1,lmax+1)))
    for l1 in ls:
        for l2 in ls:
            entry = (2*l2+1)/(4*sp.pi)*sp.sum(sp.array([\
                    Wl*(2*l3+1)*Wigner3j(l1,0,l2,0,l3,0).doit()**2\
                    for l3,Wl in zip(ls,Wls)]))
            M[l1,l2] = entry        
    print "The determinant of the MASTER matrix is", linalg.det(M)    
    cls_master = linalg.solve(M, pseudo_cls)
    print "Setting the master mono- and dipole to zero"
    cls_master[0] = 0
    cls_master[1] = 0
    return ls, cls_master

    



        

def do_harmonic_analysis(parameters,vrmap,bincenter):


    radius_largest_hole = find_largest_hole(parameters,vrmap)
    smoothing_fwhm_Mpch = parameters.smoothing_radius_fwhm # Mpc/h
    # bincenter is currently calculated as the mean of the distances to the halos in the bin
    smoothing_fwhm = smoothing_fwhm_Mpch/bincenter 
    

    
    
    lmax = int(sp.floor(sp.pi/smoothing_fwhm))
    print "The mean distance to the halos in the bin is", bincenter,"Mpc/h"
    print "The radius of the largest hole is", radius_largest_hole*bincenter, "Mpc/h"     
    print "The smoothing fwhm is", smoothing_fwhm_Mpch    
    print "The smoothing fwhm is", smoothing_fwhm, "rad."
    print "This corresponds to l = pi/%s = %s" \
        %(smoothing_fwhm,sp.pi/smoothing_fwhm)
    print "The pixel size is", pixelsize_in_radians(parameters), "rad."
    
            
    empty_pixels_total = len(vrmap[(vrmap == parameters.unseen)\
                                 | (vrmap == parameters.badval)])
    empty_pixels_fraction = empty_pixels_total/hp.nside2npix(parameters.nside)
    print "empty_pixels_fraction = ", empty_pixels_fraction
    
    
    ar_dummy, window = apply_window(parameters,vrmap,1,1,smoothing_fwhm)
    ar_masked, mask = apply_mask(parameters,vrmap,1,1)
    ar_final = hp.smoothing(ar_masked*window,fwhm=smoothing_fwhm,verbose=False)
    
    pixels_total = hp.nside2npix(parameters.nside)
    factor = pixels_total/(pixels_total-empty_pixels_total)
    print "Note: I am correcting the pseudo and the master coefficients by factor."
    print "factor = ", factor    
    
    ls, pseudo_cls = get_pseudo_powerspectrum(ar_final,lmax)    
    ls, master_cls = get_MASTER_corrected_powerspectrum(pseudo_cls,window,lmax)
        
    return ls,master_cls
    

    
    
    
def print_powerspectra_to_file(parameters,observers):
    
    f = open(parameters.powerspectrafile,'w')
    
    f.write("085\t")
    
    for l in range(parameters.lmax+1):
        f.write("%s\t" % l)
        
    for ob_number, ob in enumerate(observers):
        f.write("\n%s\t" % ob_number)

        for l in range(parameters.lmax+1):
            f.write("%s\t" % ob.cls[l])
        f.write("\n")
        
    f.close()
    
