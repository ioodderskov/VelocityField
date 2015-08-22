from __future__ import division
import matplotlib.pyplot as plt
import cPickle
import scipy as sp
from scipy.interpolate import Rbf

f = open('../cases/Planck512/parameters.save')
parameters = cPickle.load(f)
f.close()

#grid = sp.loadtxt('../cases/sim16/grid.txt')

#zs_grid = sp.array(grid[:,2])
#selection = zs_grid == 0.8

#xs_grid = sp.array(grid[:,0])[selection]
#ys_grid = sp.array(grid[:,1])[selection]
#
#rhos = sp.array(grid[:,4])[selection]
#Xs, Ys = sp.meshgrid(xs_grid,ys_grid)
#Rhos = sp.empty((len(ys_grid),len(xs_grid)))
#
#for xindex,x in xs_grid:
#    for yindex,y in ys_grid:
#        rho = rhos[]



xs = [halo.position[0] for halo in parameters.halos]
ys = [halo.position[1] for halo in parameters.halos]
vxs = [halo.velocity[0] for halo in parameters.halos]
vys = [halo.velocity[1] for halo in parameters.halos]

xs_gal = [galaxy.position[0] for galaxy in parameters.galaxies]
ys_gal = [galaxy.position[1] for galaxy in parameters.galaxies]
vxs_gal = [galaxy.velocity[0] for galaxy in parameters.galaxies]
vys_gal = [galaxy.velocity[1] for galaxy in parameters.galaxies]

plt.figure(figsize=(5,5))
ax = plt.gca()

scale = 5e4
ax.plot(xs,ys,'b.',markersize=10)
ax.quiver(xs,ys,vxs,vys,color='b',scale=scale)
ax.quiver(xs_gal,ys_gal,vxs_gal,vys_gal,color='g',scale=scale)
#ax.contour(xs_grid,ys_grid,Rhos)
