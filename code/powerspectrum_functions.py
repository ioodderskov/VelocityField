from __future__ import division
import healpy as hp
import scipy as sp
import pdb
import copy

def fill_empty_entries(parameters,ar):
    
    while parameters.badval in ar:
    
        ar_new = copy.copy(ar)
     
        for index, x in enumerate(ar):
            if x != parameters.badval:
                continue
            
            theta,phi = hp.pix2ang(parameters.nside,index)
            neighbours = hp.get_all_neighbours(parameters.nside,theta, phi=phi)
            neighbours = neighbours[ar[neighbours] != parameters.badval]
            neighbours = neighbours[ar[neighbours] != parameters.unseen]
            
            if len(neighbours) == 0:
#                print "Found a pixel with no usable neighbours"
                continue
    
                
            x_new = sp.mean(ar[neighbours])
            ar_new[index] = x_new 
            
        ar = ar_new
           
    return ar

def create_map(parameters,thetas,phis,vrs):
    
    
    pix = hp.ang2pix(parameters.nside,thetas,phis)
    number_of_pixels = hp.nside2npix(parameters.nside)
    vrmap = sp.ones(number_of_pixels)*parameters.badval
    
    vrs = sp.array(vrs)
    vrs_mean_of_repeated_pixels = copy.copy(vrs)
    for p in set(pix):
        vrs_mean_of_repeated_pixels[pix == p] = sp.mean(vrs[pix == p])
    
    vrmap[pix] = vrs_mean_of_repeated_pixels
    
    theta_max = sp.arccos(1-2*parameters.skyfraction)
    pix_all = sp.array(range(number_of_pixels))
#    pdb.set_trace()
    pix_unseen = pix_all[hp.pix2ang(parameters.nside,pix_all)[0]>theta_max]
#    pdb.set_trace()
    vrmap[pix_unseen] = parameters.unseen

    filled_vrmap = fill_empty_entries(parameters,vrmap)
    
    return filled_vrmap


def smooth_map(parameters,vrmap):
    if not (parameters.preset_smoothinglength | parameters.smooth_largest_hole):
        print "If you want to smooth the map, you need to make a choice for the smoothinglength"

    if parameters.preset_smoothinglength:
        smoothing_fwhm = parameters.smoothing_fwhm
    
    if parameters.smooth_largest_hole:
        smoothing_fwhm = find_largest_hole(parameters,vrmap)
        
    vrmap = hp.smoothing(vrmap,fwhm=smoothing_fwhm)

    return vrmap
        

def do_harmonic_analysis(parameters,vrmap):

    mask = [vrmap == parameters.unseen]
    masked_vrmap = hp.ma(vrmap)
    masked_vrmap.mask = mask

    cls = hp.anafast(masked_vrmap,lmax=parameters.lmax)
    ls = sp.arange(parameters.lmax+1)
        
    return ls,cls
    
def find_largest_hole(parameters,vrmap):
    # Maybe one could use the healpy-function "Get_all_neighbours" to find the largest hole?
    print "This functions is still empty"
    
    
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
    
