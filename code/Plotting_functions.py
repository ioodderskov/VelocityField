from __future__ import division
import scipy as sp
import sys
sys.path.insert(0,'/home/io/Dropbox/Projekter/Hubble/VelocityField/SNdata')





# This function makes the plots a bit prettier
def prettify(ax):

    for side in ['left','right','top','bottom']:
        ax.spines[side].set_linewidth(3)

#    ax.tick_params('both', length=10, width=1, which='major')
#    ax.tick_params('both', length=5, width=0.5, which='minor')




def make_hubbleplot(radial_distances,radial_velocities,mind,maxd,boxsize):

#    from gplot import Plot 
#    plt = Plot('latex_full_hubble','wojtak_half2')

    
    # Chosen plot options
#    from gplot import Plot 
#    plt = Plot('latex_full_hubble','wojtak_half')
#    plt.rc('font',family = 'serif')
    
    
    radial_distances = sp.array(radial_distances)
    radial_velocities = sp.array(radial_velocities)

#    H = 100
    c = 3e5
#    zhub = H/c*radial_distances
    
#    dL = radial_distances*(1+zhub)
#    Mv = -19.3 # (wikipedia)
#    m02 = 0.2*(5*sp.log10(dL)+Mv+25) # The apparent lumonisity is proportional to this quantity
#    log_cz = sp.log10(radial_velocities)
    
#    m02_min = 0.2*(5*sp.log10(mind*(1+H/c*mind))+Mv+25)
#    m02_max = 0.2*(5*sp.log10(maxd*(1+H/c*maxd))+Mv+25)
    
#    zo = sp.array(zo)
#    theta = sp.arccos(zo/radial_distances)    
    
#    plt.plot(m02, log_cz,'ro')
    
    # I can't figure out how to make this plot okay with gplot, so I will just use matplotlib
    import matplotlib.pyplot as plt
    plt.rc('font',family = 'serif')
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(radial_distances, radial_velocities/c,'bo') 


    plt.axis([0,300,0,0.1])
    plt.xlabel('$r [Mpc/h]$',fontsize=16)
    plt.ylabel('Redshift',fontsize=16)
    
    
#    plt.xlabel('$0.2m_V$ (mag)')
#    plt.ylabel('$\log cz$')
#    plt.axis([2.6,4.0,3.3,4.7])
    
    
#    plt.finalize(custom = True)
    
