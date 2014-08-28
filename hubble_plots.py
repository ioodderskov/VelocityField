from matplotlib import pyplot
import pylab
from mpl_toolkits.mplot3d import Axes3D
import random

def make_3D_plot(xo,yo,zo):
    fig = pylab.figure()
    ax = Axes3D(fig) 
    ax.scatter(xo, yo, zo)
    ax.set_xlim3d([-256,256])
    ax.set_ylim3d([-256,256])
    ax.set_zlim3d([-256,256])
    pyplot.show()
    
    fig = pylab.figure()
    pyplot.plot(xo,zo)
    pyplot.axes([-256,256,-256,256])
    pyplot.show()
