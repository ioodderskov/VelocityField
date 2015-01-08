import scipy as sp
import matplotlib.pyplot as plt

data = sp.loadtxt('../../cases/Planck512/powerspectra.txt')

ls = data[0,2:]
Cls = data[1:,2:]

mean_Cls = sp.mean(Cls,axis=0)
scaled_Cls = sp.array(sp.sqrt(mean_Cls*ls*(ls+1)))
plt.plot(ls,scaled_Cls)
#plt.xscale('log')
#plt.yscale('log')
plt.axis([0,20, 0, 800])
plt.show()
