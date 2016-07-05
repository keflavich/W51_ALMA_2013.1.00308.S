import numpy as np
import dust_emissivity
import radio_beam
from astropy.io import fits
from astropy import units as u
from astropy import constants
distance=5.4*u.kpc

im = fits.getdata('W51_te_continuum_best.fits')
hd = fits.getheader('W51_te_continuum_best.fits')
beam = radio_beam.Beam.from_fits_header(hd)

e2e_peak_flux = im[1350:1400,850:867].max()*u.Jy
e2e_peak_tb = e2e_peak_flux.to(u.K, beam.jtok_equiv(225*u.GHz))
print("e2e peak brightness: {0}".format(e2e_peak_tb))

e2e_luminosity = (constants.sigma_sb * (e2e_peak_tb)**4 * (4*np.pi*(beam.major*distance) * (beam.minor*distance)/(8*np.log(2)))).to(u.L_sun, u.dimensionless_angles())
print("e2e luminosity: {0}".format(e2e_luminosity))

e2e_dustmass = dust_emissivity.dust.massofsnu(225*u.GHz, e2e_peak_flux, distance=distance, temperature=e2e_peak_tb)
print("e2e dust mass: {0}".format(e2e_dustmass))
