import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
import paths
import pyregion
from astropy.io import fits
import aplpy
from constants import distance
import pylab as pl

regions = pyregion.open(paths.rpath('cores.reg'))

fig1 = pl.figure(1)
fig1.clf()
FF2 = aplpy.FITSFigure('../FITS/w51_spw3_continuum_r0_mulstiscale.image.fits', figure=fig1)
FF2.show_grayscale()
FF2.add_scalebar((1*u.pc/distance).to(u.deg, u.dimensionless_angles()).value,)
FF2.scalebar.set_label('1 pc')
FF2.save(paths.fpath("continuum.png"))
FF2.show_regions(paths.rpath('cores.reg'), layer='cores')
FF2.hide_layer('cores_txt')
FF2.save(paths.fpath("continuum_with_cores.png"))
