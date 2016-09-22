
import os
import numpy as np
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.utils.console import ProgressBar
from astropy import constants
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
import pyregion
from astropy.nddata.utils import Cutout2D
from astropy import coordinates

import warnings
warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)
warnings.filterwarnings('ignore', category=RuntimeWarning,
                        message="invalid value encountered in greater", append=True)


#contfile = fits.open(paths.dpath('selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_10mJy.image.pbcor.fits'))
#ln selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy.image.pbcor.fits W51_te_continuum_best.fits
#ln selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy.residual.fits W51_te_continuum_best_residual.fits
contfile_e2e8 = fits.open(paths.dpath('longbaseline/W51e2cax.cont.image.pbcor.fits'))
datae2e8 = u.Quantity(contfile_e2e8[0].data.squeeze(), unit=contfile_e2e8[0].header['BUNIT'])

beam_e2e8 = radio_beam.Beam.from_fits_header(contfile_e2e8[0].header)

# this section repeats the analysis of dust_proerties, but it serves as a nice
# independent check since I did it 6 months later.... give or take...
beam_radius_e2e8 = ((beam_e2e8.sr/(2*np.pi))**0.5 * masscalc.distance).to(u.pc,
                                                                          u.dimensionless_angles())

contfile_north = fits.open(paths.dpath('longbaseline/W51ncax.cont.image.pbcor.fits'))
datanorth = u.Quantity(contfile_north[0].data.squeeze(), unit=contfile_north[0].header['BUNIT'])

beam_north = radio_beam.Beam.from_fits_header(contfile_north[0].header)

# this section repeats the analysis of dust_proerties, but it serves as a nice
# independent check since I did it 6 months later.... give or take...
beam_radius_north = ((beam_north.sr/(2*np.pi))**0.5 * masscalc.distance).to(u.pc,
                                                                          u.dimensionless_angles())


# what luminosity is produced within the optically thick region?
e2epeak = np.nanmax(datae2e8[2480:2600,2360:2460])
print("e2e peak flux: {0}".format(e2epeak))
e2etbpeak = e2epeak * beam_e2e8.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("e2e peak brightness tem: {0}".format(e2etbpeak))
e2e_blackbody_luminosity = (constants.sigma_sb * e2etbpeak**4 *
                            (4*np.pi*beam_radius_e2e8**2)).to(u.L_sun)
print("e2e blackbody luminosity: {0} = {0:e}".format(e2e_blackbody_luminosity,))

e8peak = np.nanmax(datae2e8[1200:1400,2500:2700])
print("e8 peak flux: {0}".format(e8peak))
e8tbpeak = e8peak * beam_e2e8.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("e8 peak brightness tem: {0}".format(e8tbpeak))
e8_blackbody_luminosity = (constants.sigma_sb * e8tbpeak**4 *
                           (4*np.pi*beam_radius_e2e8**2)).to(u.L_sun)
print("e8 blackbody luminosity: {0} = {0:e}".format(e8_blackbody_luminosity,))


northpeak = np.nanmax(datanorth[2400:2600,2400:2600])
print("north peak flux: {0}".format(northpeak))
northtbpeak = northpeak * beam_north.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("north peak brightness tem: {0}".format(northtbpeak))
north_blackbody_luminosity = (constants.sigma_sb * northtbpeak**4 *
                              (4*np.pi*beam_radius_north**2)).to(u.L_sun)
print("north blackbody luminosity: {0} = {0:e}".format(north_blackbody_luminosity,))
