from PIL import Image
import scipy as sp
from matplotlib.image import pil_to_array
import matplotlib.pyplot as plt


# Earth image extracted from basemap:
# https://github.com/matplotlib/basemap/blob/master/lib/mpl_toolkits/basemap/data/shadedrelief.jpg
#grayscale_pil_image = Image.open("shadedrelief.jpg").convert("L")
#image_array = pil_to_array(grayscale_pil_image)

#print image_array.shape

import healpy as hp

#theta = np.linspace(0,np.pi,num=image_array.shape[0])[:,None]
#phi = np.linspace(-np.pi, np.pi, num=image_array.shape[1])

theta = sp.load('thetas.npy')
phi = sp.load('phis.npy')
vprs = sp.load('vprs.npy')

nside = 8

print "Pixel area: %.2f square degrees" %hp.nside2pixarea(nside,degrees=True)

pix = hp.ang2pix(nside,theta,phi)

healpix_map = sp.zeros(hp.nside2npix(nside),dtype=sp.double)

healpix_map[pix] = vprs

mask = healpix_map == 0 
healpix_map_masked = hp.ma(healpix_map)
healpix_map_masked.mask = mask

hp.orthview(healpix_map_masked,cmap="gray", xsize=2000, flip="geo")
plt.title("Mollweide view of the Eart")
#plt.show()




smooth_healpix_map = hp.smoothing(healpix_map, fwhm=sp.radians(10))
hp.orthview(smooth_healpix_map, cmap="gray", xsize=2000, flip="geo")


LMAX = 20
cl = hp.anafast(smooth_healpix_map, lmax=LMAX)
L = sp.arange(len(cl))

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(L, sp.sqrt(L*(L+1)*cl))
plt.xlabel('L'); plt.ylabel('$\sqrt{l(l+1)C_l}$'); plt.grid()
hp.write_cl('cl.fits', cl)
#ax.set_xscale('log')
#ax.set_yscale('log')
plt.show()



#
#hp.gnomview(healpix_map, rot=(12.5, 41.9), reso=.5, xsize=1600, cmap="gray", flip="geo")
#title("Gnomonic projection of Italy")


