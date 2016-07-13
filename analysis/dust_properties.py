import numpy as np
import pyregion
import dust_emissivity
from astropy import wcs
import paths
import radio_beam
from astropy.io import fits
from astropy import units as u
from astropy import constants
distance=5.4*u.kpc

im = fits.getdata(paths.dpath('W51_te_continuum_best.fits'))
hd = fits.getheader(paths.dpath('W51_te_continuum_best.fits'))
beam = radio_beam.Beam.from_fits_header(hd)

freq = 225*u.GHz

e2e_peak_flux = im[1350:1400,850:867].max()*u.Jy
e2e_peak_tb = e2e_peak_flux.to(u.K, beam.jtok_equiv(freq))
print("e2e peak brightness: {0}".format(e2e_peak_tb))

e2e_luminosity = (constants.sigma_sb * (e2e_peak_tb)**4 *
                  (4*np.pi*(beam.major*distance) *
                   (beam.minor*distance)/(8*np.log(2)))).to(u.L_sun,
                                                            u.dimensionless_angles())
print("e2e luminosity: {0}".format(e2e_luminosity))

e2e_dustmass = dust_emissivity.dust.massofsnu(freq, e2e_peak_flux,
                                              distance=distance,
                                              temperature=e2e_peak_tb)
print("e2e dust mass: {0}".format(e2e_dustmass))

e2e_peak_column = dust_emissivity.dust.colofsnu(freq, e2e_peak_flux,
                                                beamomega=beam,
                                                temperature=100*u.K).to(u.cm**-2,
                                                                        u.dimensionless_angles())
print("e2e peak column: {0}".format(e2e_peak_column))

regfn = 'e2_exclude_e2w.reg'
reg_noe2w = pyregion.open(paths.rpath(regfn))
regfn = 'e2_exclude_e2w_and_e2e.reg'
reg_noe2w_or_e2e = pyregion.open(paths.rpath(regfn))

mywcs = wcs.WCS(hd)
pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5
ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam

mask = reg_noe2w.get_mask(header=hd, shape=im.shape)

e2_total_flux = im[mask].sum()*u.Jy/u.beam / ppbeam
e2_median_flux = np.nanmedian(im[mask])*u.Jy
e2_median_tb = e2_median_flux.to(u.K, beam.jtok_equiv(freq))
print("Median brightness temperature: {0}".format(e2_median_tb))
e2_20k_total_mass = dust_emissivity.dust.massofsnu(freq, e2_total_flux,
                                                   distance=distance,
                                                   temperature=20*u.K)

e2_median_column = dust_emissivity.dust.colofsnu(freq, e2_median_flux,
                                                  beamomega=beam,
                                                  temperature=100*u.K).to(u.cm**-2,
                                                                          u.dimensionless_angles())
print("e2 median column: {0}".format(e2_median_column))

print("e2 total mass, excluding e2w's flux, assuming T={1}: {0}"
      .format(e2_20k_total_mass, 20*u.K))

e2_100k_total_mass = dust_emissivity.dust.massofsnu(freq, e2_total_flux,
                                                    distance=distance,
                                                    temperature=100*u.K)

print("e2 total mass, excluding e2w's flux, assuming T={1}: {0}"
      .format(e2_100k_total_mass, 100*u.K))

mask = reg_noe2w_or_e2e.get_mask(header=hd, shape=im.shape)

e2_noe2e_total_flux = im[mask].sum()*u.Jy/u.beam / ppbeam
e2_noe2e_median_flux = np.nanmedian(im[mask])*u.Jy
e2_noe2e_median_tb = e2_noe2e_median_flux.to(u.K, beam.jtok_equiv(freq))
print("Median brightness temperature (no e2e either): {0}".format(e2_noe2e_median_tb))
e2_noe2e_20k_total_mass = dust_emissivity.dust.massofsnu(freq,
                                                         e2_noe2e_total_flux,
                                                         distance=distance,
                                                         temperature=20*u.K)

print("e2 total mass, excluding e2w and e2e flux, assuming T={1}: {0}"
      .format(e2_noe2e_20k_total_mass, 20*u.K))

e2_noe2e_100k_total_mass = dust_emissivity.dust.massofsnu(freq,
                                                          e2_noe2e_total_flux,
                                                          distance=distance,
                                                          temperature=100*u.K)

print("e2 total mass, excluding e2w and e2e flux, assuming T={1}: {0}"
      .format(e2_noe2e_100k_total_mass, 100*u.K))
