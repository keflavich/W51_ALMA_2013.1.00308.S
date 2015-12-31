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

data = fits.open(paths.dpath('w51_spw3_continuum_r0_mulstiscale.image.fits'))
