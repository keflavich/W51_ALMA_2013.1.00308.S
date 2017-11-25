import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
import paths
import pyregion
from astropy.io import fits
import aplpy
from constants import distance, continuum_frequency
from astropy import wcs
import radio_beam
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

ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')

ax.axis([626, 2513, 588, 2385])
fig1.savefig(paths.fpath("continuum.png"), bbox_inches='tight')

ra = [reg.coord_list[0] for reg in regions if len(reg.coord_list) == 3]
dec = [reg.coord_list[1] for reg in regions if len(reg.coord_list) == 3]
size = [100 * 3600 * reg.coord_list[2] for reg in regions if len(reg.coord_list) == 3]
scat = ax.scatter(ra, dec, s=size, transform=ax.get_transform('world'),
                  c='none', edgecolors='r', )

#ax.axis([290.93636, 290.90752, 14.496992, 14.524811], transform=ax.get_transform('world'))
ax.axis([626, 2513, 588, 2385])
fig1.savefig(paths.fpath("continuum_with_cores.png"), bbox_inches='tight')


scat.set_visible(False)
beam = radio_beam.Beam.from_fits_header(fh[0].header)
data_K = (fh[0].data * u.Unit(fh[0].header['BUNIT']) * u.beam).to(u.K, beam.jtok_equiv(continuum_frequency))

im = ax.imshow(data_K, cmap=cm, origin='lower', interpolation='none', vmin=0,
               vmax=252)
cont = ax.contour(data_K, colors=[(0,0,0,0), 'g','b','m','c','y','w','r', (0,1,0,1)], levels=[0, 50, 100,125,150,175,200,225,250])
cb = fig1.colorbar(mappable=im)
fig1.colorbar(mappable=cont, cax=cb.ax)
cb.set_label("$T_B$ [K]")

# e2 big
ax.axis((816,921,1335,1436))
fig1.savefig(paths.fpath('continuum_e2_zoom_TB.png'), bbox_inches='tight')

# e2 tight
ax.axis((846,881,1365,1396))
fig1.savefig(paths.fpath('continuum_e2_zoom_tight_TB.png'), bbox_inches='tight')

# e8
ax.axis((866,890,1226,1270))
fig1.savefig(paths.fpath('continuum_e8_zoom_TB.png'), bbox_inches='tight')

# north tight
ax.axis((1980,2017,1983,2013))
fig1.savefig(paths.fpath('continuum_north_zoom_tight_TB.png'), bbox_inches='tight')

# north broad
ax.axis((1961,2105,1966,2016))
fig1.savefig(paths.fpath('continuum_north_zoom_TB.png'), bbox_inches='tight')
