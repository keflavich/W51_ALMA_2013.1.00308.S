import os
import numpy as np
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy import log
import paths
from astropy.io import fits
from astropy.table import Column, Table
import radio_beam
import astrodendro
from astropy import wcs
import masscalc
import photutils
from astropy.nddata.utils import Cutout2D
from astropy import coordinates


#contfile = fits.open(paths.dpath('selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_10mJy.image.pbcor.fits'))
#ln selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy.image.pbcor.fits W51_te_continuum_best.fits
#ln selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy.residual.fits W51_te_continuum_best_residual.fits
contfile = fits.open(paths.dpath('W51_te_continuum_best.fits'))
data = u.Quantity(contfile[0].data, unit=contfile[0].header['BUNIT'])
mywcs = wcs.WCS(contfile[0].header)
beam = radio_beam.Beam.from_fits_header(contfile[0].header)
pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
pixel_scale_as = pixel_scale.to(u.arcsec).value
ppbeam = (beam.sr/(pixel_scale**2)).decompose().value / u.beam

#pl.hist(data[np.isfinite(data)], bins=np.linspace(1e-4,0.1,100))

definitely_signal = data > 10*u.mJy/u.beam
total_signal = data[definitely_signal].sum() / ppbeam
print("Total flux: {0}".format(total_signal))
print("Total mass(20K): {0}".format(total_signal * masscalc.mass_conversion_factor()*u.M_sun/u.Jy))

dendro_merge = Table.read(paths.tpath('dendro_merge_continuum_and_line.ipac'), format='ascii.ipac')
corelike = dendro_merge['corelike'] == 'True'
print("Total protostar flux (0.2): {0}".format(dendro_merge['cont_flux0p2arcsec'].sum()))
print("Total protostar flux (peak): {0}".format(dendro_merge['peak_cont_flux'].sum()))
