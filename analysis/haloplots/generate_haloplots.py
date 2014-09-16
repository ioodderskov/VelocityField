

from __future__ import division
#import matplotlib.pyplot as plt
import haloplot_functions as hp
from gplot import Plot # Import gplot
#
plt = Plot('latex_full_hubble','halo') 
plt.rc('font',family = 'serif')

indir='/home/io/Dropbox/Projekter/Hubble/VelocityField/cases'
outdir='/home/io/Dropbox/SharedStuff/hubble2013'
#outdir='/home/io/Dropbox/PHD/Python/tests'
outname = outdir+'/haloplots.png'

#plt.figure(figsize=(8.267716535433072, 10))

adj = 'box'
#adj = 'datalim'

ax1 = plt.subplot(2,2,1)
infil = indir+'/Planck512/parents_11'
text1 = "A.0, A.1 and A.2: Standard simulation"
text2 = "Random and Local Group like halos as observers"
ax1.set_aspect('equal', adjustable=adj)
hp.haloplot(infil, 0, 0, 512, text1, text2)

print 'Done with plot 1'

ax2 = plt.subplot(2,2,2)
infil = indir+'/Planck512/parents_11'
text1 = "A.3: Standard simulation"
text2 = "One cone fraction"
ax2.set_aspect('equal', adjustable=adj)
hp.haloplot(infil, 1, 0, 512, text1, text2)

print 'Done with plot 2'


ax3 = plt.subplot(2,2,3)
infil = indir+'/Planck512/parents_11'
text1 = "A.4: Standard simulation"
text2 = "Two cones fraction"
ax3.set_aspect('equal', adjustable=adj)
hp.haloplot(infil, 2, 0, 512, text1, text2)

print 'Done with plot 3'

ax4 = plt.subplot(2,2,4)
infil = indir+'/Planck512_lightcone/parents_0'
text1 = "A.5: Standard simulation" 
text2 = "Halos from past lightcone of observer"
ax4.set_aspect('equal', adjustable=adj)
hp.haloplot(infil, 0, 1, 512, text1, text2)

print 'Done with plot 4'


plt.finalize(custom = True)
plt.savefig(outname, dpi=200)
plt.show()




