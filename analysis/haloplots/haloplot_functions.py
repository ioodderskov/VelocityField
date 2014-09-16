# -*- coding: utf-8 -*-
"""
Created on Mon May 19 17:37:40 2014

@author: io
"""

from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt



def haloplot(fil,cones,lightcone,box,name,title):
    

#    fig = plt.figure()
    
    # Load distances and Hubbleconstants, and count the maximal number of observers
    data = sp.loadtxt(fil) 
    
    
    masses = data[0:10,2]
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
    
#    ax = plt.add_subplot(1,1,1)
    plt.scatter(x,y,marker='.',c='b', s=10*m )
    plt.plot(box/2,box/2,c='g',marker='*',markersize=10)
    #axis('equal')
    plt.axis([0,box,0,box])
#    plt.set_aspect('equal', adjustable='box')
    
    
#    plt.xlabel('$x / [Mpc/h]$', fontsize=15)
#    plt.ylabel('$y / [Mpc/h]$', fontsize=15)
    plt.title(name+'\n'+title, fontsize=8)
    plt.xlabel('$x$ [Mpc/h]')
    plt.ylabel('$y$ [Mpc/h]')
#    plt.title(name+'\n'+title)
    
    
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
            plt.plot(box/2+r*sp.cos(angles_c),box/2+r*sp.sin(angles_c),'g',alpha=0.4)
        elif cones == 2:
            plt.plot(box/2+r*sp.cos(angles_c1),box/2+r*sp.sin(angles_c1),'g',alpha=0.4)
            plt.plot(box/2+r*sp.cos(angles_c2),box/2+r*sp.sin(angles_c2),'g',alpha=0.4)
        else:
            plt.plot(box/2+r*sp.cos(angles),box/2+r*sp.sin(angles),'g',alpha=0.4)
            
    
    # Plotting the cones
    if cones != 0:
        a = sp.tan(an_patch)
        b = box/2*(1-a)
        y_cone=a*x_cone+b
        index = sp.sqrt((x_cone-box/2)**2+(y_cone-box/2)**2) < maxd
        xc = x_cone[index]
        yc = y_cone[index]
        plt.plot(x_cone,box/2*sp.ones(len(x_cone)),'g')
        plt.plot(xc,yc,'g')
        
        return(0)
    
    
#    plt.savefig(outname)
#    
