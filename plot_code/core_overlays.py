import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
import paths
import pyregion
from astropy.io import fits
import aplpy
from constants import distance
from astropy import wcs
import pylab as pl

regions = pyregion.open(paths.rpath('cores.reg'))

fh = fits.open(paths.dpath("w51_te_continuum_best.fits"))

fig1 = pl.figure(1)
fig1.clf()
# FF2 = aplpy.FITSFigure(paths.dpath("w51_te_continuum_best.fits"), figure=fig1)
# FF2.show_grayscale(invert=True)
# FF2.recenter(290.9230429,14.51167269,50/3600.)
# FF2.add_scalebar((1*u.pc/distance).to(u.deg, u.dimensionless_angles()).value,)
# FF2.scalebar.set_label('1 pc')
# FF2.save(paths.fpath("continuum.png"))
# FF2.show_regions(paths.rpath('cores.reg'), layer='cores')
# FF2.hide_layer('cores_txt')
# FF2.save(paths.fpath("continuum_with_cores.png"))

mywcs = wcs.WCS(fh[0].header)
ax = pl.subplot(projection=mywcs)
cm = pl.cm.gray_r
cm.set_over('black')
cm.set_bad('white')
ax.imshow(fh[0].data, cmap=cm, origin='lower', interpolation='none',
          vmin=-0.001, vmax=0.01)
ax.axis([626, 2513, 588, 2385])
fig1.savefig(paths.fpath("continuum.png"))

ra = [reg.coord_list[0] for reg in regions if len(reg.coord_list) == 3]
dec = [reg.coord_list[1] for reg in regions if len(reg.coord_list) == 3]
size = [100 * 3600 * reg.coord_list[2] for reg in regions if len(reg.coord_list) == 3]
ax.scatter(ra, dec, s=size, transform=ax.get_transform('world'), c='none',
           edgecolors='r', )

#ax.axis([290.93636, 290.90752, 14.496992, 14.524811], transform=ax.get_transform('world'))
ax.axis([626, 2513, 588, 2385])
fig1.savefig(paths.fpath("continuum_with_cores.png"))
