from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt

d = sp.load("../cases/sim16/line_of_sight_vector.npy")
x = d[0,:]
psi = sp.load("../cases/sim16/psi_along_the_line_of_sight.npy")
potential = sp.load("../cases/sim16/potential_array.npy")



potential_along_x_axis = potential[:,0,0]
x_pot = sp.linspace(1,15,8)

plt.plot(x_pot,potential_along_x_axis)
plt.plot(x,psi)
plt.show()
