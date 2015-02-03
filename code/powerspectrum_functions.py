from __future__ import division
import healpy as hp
import scipy as sp
import pdb
import copy
from collections import Counter

def find_largest_hole(parameters,ar):
    
    minimal_distances = []
    all_pixels = sp.array(range(len(ar)))
    nonempty_pixels = all_pixels[(ar[all_pixels] != parameters.badval) & (ar[all_pixels] != parameters.unseen)]

        
    for p in nonempty_pixels:

        minimal_distance = 3.14
        theta,phi = hp.pix2ang(parameters.nside,p)
        
        for p_i in nonempty_pixels:

            if p_i == p:
                continue
            
            theta_i, phi_i = hp.pix2ang(parameters.nside,p_i)
            angular_distance = hp.rotator.angdist([theta,phi],[theta_i,phi_i])
            minimal_distance = sp.minimum(minimal_distance,angular_distance)
            
        minimal_distances.append(minimal_distance)
        
    radius_of_largest_hole = max(minimal_distances)

    print "radius of largest hole = ", radius_of_largest_hole  
    return radius_of_largest_hole
            


def fill_empty_entries(parameters,ar):
    
#    pixels_without_neighbours = []
    
    while parameters.badval in ar:
    
        ar_new = copy.copy(ar)
     
        for index, x in enumerate(ar):
            if x != parameters.badval:
                continue
            
            theta,phi = hp.pix2ang(parameters.nside,index)
            neighbours = hp.get_all_neighbours(parameters.nside,theta, phi=phi)
            neighbours = neighbours[ar[neighbours] != parameters.badval]
            neighbours = neighbours[ar[neighbours] != parameters.unseen]
            
#            if len(neighbours) == 0:
#                pixels_without_neighbours.append(index)
#                continue
    
                
            x_new = sp.mean(ar[neighbours])
            ar_new[index] = x_new
            
        ar = ar_new
        
#    if len(pixels_without_neighbours) != 0:
#        print "There are some pixels with no neighbours"
#        pixel_most_common = Counter(pixels_without_neighbours).most_common(1)[0][0]
#        ar[pixel_most_common] = 1000
##        pdb.set_trace()
#        recurrence_most_common = Counter(pixels_without_neighbours).most_common(1)[0][1]
#        print "recurrence_most_common = ", recurrence_most_common
           
    return ar

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
#    pdb.set_trace()
    return vrmap



def Gaussian_kernel(angular_distance,b):
    K = sp.exp(-(angular_distance)**2/(2*b**2))
    return K
    

def Ios_smoothing(parameters,ar,b):
    
    ar_new = sp.zeros_like(ar)
    
    for index, x in enumerate(ar):
        numerator = 0
        denominator = 0
        
        theta, phi = hp.pix2ang(parameters.nside,index)
        
        
        for index_i, xi in enumerate(ar):
            
            if xi == parameters.badval:
                continue
            
            theta_i, phi_i = hp.pix2ang(parameters.nside,index_i)
            angular_distance = hp.rotator.angdist([theta, phi],[theta_i, phi_i])
            kernel = Gaussian_kernel(angular_distance,b)
            numerator = numerator + kernel*xi
            denominator = denominator + kernel
        
        x_new = numerator/denominator            
        ar_new[index] = x_new
        
    
    return ar_new


def smooth_map(parameters,vrmap):
    if not (parameters.preset_smoothinglength | parameters.smooth_largest_hole):
        print "If you want to smooth the map, you need to make a choice for the smoothinglength"

    if parameters.preset_smoothinglength:
        smoothing_fwhm = parameters.smoothing_fwhm
    
    if parameters.smooth_largest_hole:
        smoothing_fwhm = find_largest_hole(parameters,vrmap)
        
    #vrmap = Ios_smoothing(parameters,vrmap,smoothing_fwhm)
    vrmap = hp.smoothing(vrmap,fwhm=smoothing_fwhm)
        

    return vrmap
        

def do_harmonic_analysis(parameters,vrmap):

    mask = [vrmap == parameters.unseen]
    masked_vrmap = hp.ma(vrmap)
    masked_vrmap.mask = mask

    cls = hp.anafast(masked_vrmap,lmax=parameters.lmax)
    ls = sp.arange(parameters.lmax+1)
        
    return ls,cls
    

    
    
    
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
    
