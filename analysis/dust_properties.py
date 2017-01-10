# see also: total_mass_analysis
import numpy as np
import pyregion
import dust_emissivity
from astropy import wcs
import paths
import radio_beam
from astropy.io import fits
from astropy import units as u
from astropy import constants
from masscalc import distance, centerfreq as freq

im = fits.getdata(paths.dpath('W51_te_continuum_best.fits'))
hd = fits.getheader(paths.dpath('W51_te_continuum_best.fits'))
beam = radio_beam.Beam.from_fits_header(hd)

print("Dust opacity: ",dust_emissivity.dust.kappa(freq))

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
print("e2e peak column (T=100K): {0}".format(e2e_peak_column))

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
e2_max_flux = np.nanmax(im[mask])*u.Jy
e2_max_tb = e2_max_flux.to(u.K, beam.jtok_equiv(freq))
print("e2 cluster peak brightness temperature: {0}".format(e2_max_tb))
e2_20k_total_mass = dust_emissivity.dust.massofsnu(freq, e2_total_flux,
                                                   distance=distance,
                                                   temperature=20*u.K)

e2_median_column = dust_emissivity.dust.colofsnu(freq, e2_median_flux,
                                                 beamomega=beam,
                                                 temperature=100*u.K).to(u.cm**-2,
                                                                         u.dimensionless_angles())
print("e2 median column (T=100K): {0}".format(e2_median_column))

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





reg_e8 = pyregion.open(paths.rpath('e8cluster.reg'))
mask = reg_e8.get_mask(header=hd, shape=im.shape)

e8_total_flux = im[mask].sum()*u.Jy/u.beam / ppbeam
e8_median_flux = np.nanmedian(im[mask])*u.Jy
e8_median_tb = e8_median_flux.to(u.K, beam.jtok_equiv(freq))
print("e8 cluster Median brightness temperature: {0}".format(e8_median_tb))
e8_max_flux = np.nanmax(im[mask])*u.Jy
e8_max_tb = e8_max_flux.to(u.K, beam.jtok_equiv(freq))
print("e8 cluster peak brightness temperature: {0}".format(e8_max_tb))
e8_20k_total_mass = dust_emissivity.dust.massofsnu(freq, e8_total_flux,
                                                   distance=distance,
                                                   temperature=20*u.K)

e8_median_column = dust_emissivity.dust.colofsnu(freq, e8_median_flux,
                                                 beamomega=beam,
                                                 temperature=100*u.K).to(u.cm**-2,
                                                                         u.dimensionless_angles())
print("e8 cluster median column: {0}".format(e8_median_column))

print("e8 cluster total mass, assuming T={1}: {0}"
      .format(e8_20k_total_mass, 20*u.K))

e8_100k_total_mass = dust_emissivity.dust.massofsnu(freq, e8_total_flux,
                                                    distance=distance,
                                                    temperature=100*u.K)

print("e8 cluster total mass, assuming T={1}: {0}"
      .format(e8_100k_total_mass, 100*u.K))

e8_peak_flux = im[mask].max()*u.Jy
e8_peak_tb = e8_peak_flux.to(u.K, beam.jtok_equiv(freq))

e8_dustmass = dust_emissivity.dust.massofsnu(freq, e8_peak_flux,
                                             distance=distance,
                                             temperature=e8_peak_tb)
print("e8 (core) dust mass: {0}".format(e8_dustmass))






reg_north = pyregion.open(paths.rpath('north_exclude_HII.reg'))
mask = reg_north.get_mask(header=hd, shape=im.shape)

north_total_flux = im[mask].sum()*u.Jy/u.beam / ppbeam
north_median_flux = np.nanmedian(im[mask])*u.Jy
north_median_tb = north_median_flux.to(u.K, beam.jtok_equiv(freq))
print("north cluster Median brightness temperature: {0}".format(north_median_tb))
north_max_flux = np.nanmax(im[mask])*u.Jy
north_max_tb = north_max_flux.to(u.K, beam.jtok_equiv(freq))
print("north cluster peak brightness temperature: {0}".format(north_max_tb))
north_20k_total_mass = dust_emissivity.dust.massofsnu(freq, north_total_flux,
                                                      distance=distance,
                                                      temperature=20*u.K)

north_median_column = dust_emissivity.dust.colofsnu(freq, north_median_flux,
                                                    beamomega=beam,
                                                    temperature=100*u.K).to(u.cm**-2,
                                                                            u.dimensionless_angles())
print("north cluster median column: {0}".format(north_median_column))

print("north cluster total mass, assuming T={1}: {0}"
      .format(north_20k_total_mass, 20*u.K))

north_100k_total_mass = dust_emissivity.dust.massofsnu(freq, north_total_flux,
                                                       distance=distance,
                                                       temperature=100*u.K)

print("north cluster total mass, assuming T={1}: {0}"
      .format(north_100k_total_mass, 100*u.K))

north_peak_flux = im[mask].max()*u.Jy
north_peak_tb = north_peak_flux.to(u.K, beam.jtok_equiv(freq))

north_dustmass = dust_emissivity.dust.massofsnu(freq, north_peak_flux,
                                                distance=distance,
                                                temperature=north_peak_tb)
print("north (core) dust mass: {0}".format(north_dustmass))


print()

reg_d2 = pyregion.open(paths.rpath('d2.reg'))
mask = reg_d2.get_mask(header=hd, shape=im.shape)

d2_total_flux = im[mask].sum()*u.Jy/u.beam / ppbeam
d2_median_flux = np.nanmedian(im[mask])*u.Jy
d2_median_tb = d2_median_flux.to(u.K, beam.jtok_equiv(freq))
print("d2 Median brightness temperature: {0}".format(d2_median_tb))
d2_max_flux = np.nanmax(im[mask])*u.Jy
d2_max_tb = d2_max_flux.to(u.K, beam.jtok_equiv(freq))
print("d2 peak brightness temperature: {0}".format(d2_max_tb))
d2_20k_total_mass = dust_emissivity.dust.massofsnu(freq, d2_total_flux,
                                                      distance=distance,
                                                      temperature=20*u.K)

d2_median_column = dust_emissivity.dust.colofsnu(freq, d2_median_flux,
                                                    beamomega=beam,
                                                    temperature=100*u.K).to(u.cm**-2,
                                                                            u.dimensionless_angles())
print("d2 median column: {0}".format(d2_median_column))

print("d2 total mass, assuming T={1}: {0}"
      .format(d2_20k_total_mass, 20*u.K))

d2_100k_total_mass = dust_emissivity.dust.massofsnu(freq, d2_total_flux,
                                                       distance=distance,
                                                       temperature=100*u.K)

print("d2 total mass, assuming T={1}: {0}"
      .format(d2_100k_total_mass, 100*u.K))

d2_peak_flux = im[mask].max()*u.Jy
d2_peak_tb = d2_peak_flux.to(u.K, beam.jtok_equiv(freq))

d2_dustmass = dust_emissivity.dust.massofsnu(freq, d2_peak_flux,
                                                distance=distance,
                                                temperature=d2_peak_tb)
print("d2 (core) dust mass (should be same as above): {0}".format(d2_dustmass))



"""
e2e peak brightness: 228.1267575167404 K
e2e luminosity: 23918.78717599044 solLum
e2e dust mass: 18.742076274056984 solMass
e2e peak column: 6.25229288904971e+25 1 / cm2
Median brightness temperature: 7.002341281739255 K
e2 cluster peak brightness temperature: 228.1267575167404 K
e2 median column: 1.9178648670620061e+24 1 / cm2
e2 total mass, excluding e2w's flux, assuming T=20.0 K: 12127.9410524976 solMass
e2 total mass, excluding e2w's flux, assuming T=100.0 K: 1931.8484592796835 solMass
Median brightness temperature (no e2e either): 7.4959131800605325 K
e2 total mass, excluding e2w and e2e flux, assuming T=20.0 K: 8808.180625349396 solMass
e2 total mass, excluding e2w and e2e flux, assuming T=100.0 K: 1403.0468989391943 solMass
e8 cluster Median brightness temperature: 2.6771436081188913 K
e8 cluster peak brightness temperature: 208.06306723244958 K
e8 cluster median column: 7.345582708106389e+23 1 / cm2
e8 cluster total mass, assuming T=20.0 K: 10789.704064105663 solMass
e8 cluster total mass, assuming T=100.0 K: 1718.6822996807393 solMass
north cluster Median brightness temperature: 3.6198044207380513 K
north cluster peak brightness temperature: 260.4098839355622 K
north cluster median column: 9.909858853481762e+23 1 / cm2
north cluster total mass, assuming T=20.0 K: 12598.129830473512 solMass
north cluster total mass, assuming T=100.0 K: 2006.743302980072 solMass
"""
