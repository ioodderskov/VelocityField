# -*- coding: utf-8 -*-
"""
Created on Mon May 19 17:37:40 2014

@author: io
"""

from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import sys


## Arguments:
#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)
#print 'fil=', sys.argv[1]
#print 'outname=', sys.argv[2]
#print 'Number of cones =', sys.argv[3]
#print 'Lightcone? =', sys.argv[4]
#print 'Box =', sys.argv[5]
#print 'Name =', sys.argv[6]
#print 'Title =', sys.argv[7]
#
#
#
#fil = sys.argv[1]
#outname = sys.argv[2]
#cones = int(sys.argv[3])
#lightcone = int(sys.argv[4])
#box = float(sys.argv[5])
#name = sys.argv[6]
#title = sys.argv[7]
#
##fil = '/home/io/Desktop/PHD/Hubble/halofinders/Rockstar/Planck512/parents_11'
##outname = '/home/io/Dropbox/hubble2013/haloplot_Planck512.png'
##cones = 0
#	
#plt.close('all')


def haloplots(fil,outname,cones,lightcone,box,name,title):
    

#    fig = plt.figure()
    
    # Load distances and Hubbleconstants, and count the maximal number of observers
    data = sp.loadtxt(fil) 
    
    
    masses = data[0:,2]
    #mmax = sp.amax(masses)
    mnorm = 2.5e15
    m = masses/mnorm
    x1 = data[0:,8]
    y1 = data[0:,9]
    z1 = data[0:,10]
    
    if lightcone == 1:
        x1 = x1+box/2
        y1 = y1+box/2
        z1 = z1+box/2
        
    index = z1 < 10
    x = x1[index]
    y = y1[index]
    
    ax = plt.add_subplot(1,1,1)
    ax.scatter(x,y,marker='o',c='b', s=10*m )
    ax.plot(box/2,box/2,c='g',marker='*',markersize=10)
    #axis('equal')
    ax.axis([0,box,0,box])
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel('$x / [Mpc/h]$', fontsize=15)
    plt.ylabel('$y / [Mpc/h]$', fontsize=15)
    plt.title(name+'\n'+title, fontsize=15)
    
    
    # Calculating bin-distances
    mind = 30
    width = 37.3
    maxd = 256
    
    mindz = mind
    maxdz = mind+width
    z = 0
    rmax = sp.zeros(100)
    while maxdz <= maxd:
        rmax[z] = maxdz
        z = z+1
        newmax = (2*maxdz**3-mindz**3)**(1/3)
        mindz = maxdz
        maxdz = newmax
    
    rmax = sp.trim_zeros(rmax)
    
    
    # Setting angle ranges for the case of 1, 2 or no cones
    if cones == 1:
        x_cone=sp.linspace(box/2,box,100)
        an_patch = 0.2*sp.pi
        angles_c = sp.linspace(0,an_patch,10)
    elif cones == 2:
        x_cone=sp.linspace(0,box,100)
        an_patch = 0.2805*sp.pi/2
        angles_c1 = sp.linspace(0,an_patch,10)
        angles_c2 = angles_c1+sp.pi
    else:
        angles = sp.linspace(0,2*sp.pi,100)    
    
            
        
    # Plotting a circle for every 5th bin
    for i in range(0,len(rmax),5):
        r=rmax[i]
        if cones == 1:
            ax.plot(box/2+r*sp.cos(angles_c),box/2+r*sp.sin(angles_c),'g',alpha=0.4)
        elif cones == 2:
            ax.plot(box/2+r*sp.cos(angles_c1),box/2+r*sp.sin(angles_c1),'g',alpha=0.4)
            ax.plot(box/2+r*sp.cos(angles_c2),box/2+r*sp.sin(angles_c2),'g',alpha=0.4)
        else:
            ax.plot(box/2+r*sp.cos(angles),box/2+r*sp.sin(angles),'g',alpha=0.4)
            
    
    # Plotting the cones
    if cones != 0:
        a = sp.tan(an_patch)
        b = box/2*(1-a)
        y_cone=a*x_cone+b
        index = sp.sqrt((x_cone-box/2)**2+(y_cone-box/2)**2) < maxd
        xc = x_cone[index]
        yc = y_cone[index]
        ax.plot(x_cone,box/2*sp.ones(len(x_cone)),'g')
        ax.plot(xc,yc,'g')
    
    
#    plt.savefig(outname)
#    
