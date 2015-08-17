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

def apply_window(parameters,ar,smoothing_radius):

        
    window = sp.ones_like(ar)
    window[ar == parameters.unseen] = -1            
            
        
    window = hp.smoothing(window,fwhm=smoothing_radius)
    window[ar == parameters.unseen] = 0
        
    window[window > 1] = 1
    thetas = [hp.pix2ang(parameters.nside,pix)[0] for pix in range(len(window))]

    theta_set = sorted(list(set(thetas)))        
    for i,theta in zip(range(1,len(theta_set)),theta_set[1:]):
        window_value_on_previous_ring = window[thetas == theta_set[i-1]][0]
        if window_value_on_previous_ring == 1:
            window[thetas == theta] = 1
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

    return ls, cls_master



def fill_holes_and_smooth(parameters,ar,smoothing_radius):

    while parameters.badval in ar:
        ar_masked, mask = apply_mask(parameters,ar,1,1)
        ar_smoothed = hp.smoothing(ar_masked,smoothing_radius)

    
        empty_pixels = [p for p in range(len(ar)) if ar[p] == parameters.badval]
        for empty_pixel in empty_pixels:
            
            theta,phi = hp.pix2ang(parameters.nside,empty_pixel)
            neighbours = hp.get_all_neighbours(parameters.nside,theta, phi)
            neighbours = neighbours[mask[neighbours] != 1]
            if len(neighbours) == 0:
#                print "Found a pixel with no usable neighbours"
                continue
            
            neighbours = sp.arange(0,len(ar))[neighbours]
            weights = [hp.rotator.angdist(hp.pix2ang(parameters.nside,empty_pixel),
                                          hp.pix2ang(parameters.nside,neighbour))\
                                          for neighbour in neighbours]    
            neighbours_vals = sp.reshape(ar_smoothed[neighbours],
                                         sp.shape(weights))

            ar[empty_pixel] = sp.average(neighbours_vals,weights=weights)

    # Now only unseen need to be masked       
    ar_masked, mask = apply_mask(parameters,ar,0,1)       
    ar = hp.smoothing(ar_masked,smoothing_radius)
    return ar
    



        

def do_harmonic_analysis(parameters,vrmap):


    smoothing_radius = find_largest_hole(parameters,vrmap)*2
    
    lmax = int(sp.floor(sp.pi/smoothing_radius))
    print "The smoothing length is", smoothing_radius, "rad."
    print "This corresponds to l = pi/%s = %s" \
        %(smoothing_radius,sp.pi/smoothing_radius)
    print "The pixel size is", pixelsize_in_radians(parameters), "rad."
    
            
    empty_pixels_total = len(vrmap[(vrmap == parameters.unseen)\
                                 | (vrmap == parameters.badval)])
    empty_pixels_fraction = empty_pixels_total/hp.nside2npix(parameters.nside)
    print "empty_pixels_fraction = ", empty_pixels_fraction
    
    
    # Get the window for the skycover
    ar_dummy, window = apply_window(parameters,vrmap,smoothing_radius)
    ar_filled_and_smoothed = fill_holes_and_smooth(parameters,vrmap,smoothing_radius)
    ar_final = window*ar_filled_and_smoothed
    
    
    factor = hp.nside2npix(parameters.nside)/sp.sum(window)
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
    
