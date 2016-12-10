# see also: dust_properties
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
contfile = fits.open(paths.dpath('W51_te_continuum_best.fits'))
data = u.Quantity(contfile[0].data, unit=contfile[0].header['BUNIT'])
mywcs = wcs.WCS(contfile[0].header)
beam = radio_beam.Beam.from_fits_header(contfile[0].header)
pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
pixel_scale_as = pixel_scale.to(u.arcsec).value
ppbeam = (beam.sr/(pixel_scale**2)).decompose().value / u.beam

#pl.hist(data[np.isfinite(data)], bins=np.linspace(1e-4,0.1,100))

# over what threshold are we including flux when measuring total masses?
threshold = 10*u.mJy/u.beam
# threshold above which 20 K is very thick
thick_threshold = (20*u.K).to(u.mJy, u.brightness_temperature(beam, masscalc.centerfreq)) / u.beam
threshold_column = (threshold * u.beam/u.Jy * masscalc.col_conversion_factor(beam)).to(u.cm**-2)
threshold_column = masscalc.dust.colofsnu(nu=masscalc.centerfreq,
                                          snu=threshold*u.beam,
                                          beamomega=beam.sr).to(u.cm**-2,
                                                                u.dimensionless_angles())
threshold_density = (masscalc.mass_conversion_factor(20) * (threshold*u.beam).to(u.mJy).value /
                     (4/3.*np.pi) /
                     (beam.sr.value*masscalc.distance**2)**(1.5) /
                     (2.8*constants.m_p)).to(1/u.cm**3)
definitely_signal = data > threshold
definitely_thick_if_20K = data > thick_threshold
total_signal = data[definitely_signal].sum() / ppbeam
print("Total pixels > 10mJy/beam: {0} = {1}; r_eff = {2}"
      .format(definitely_signal.sum(), definitely_signal.sum()/ppbeam,
              ((definitely_signal.sum()*(pixel_scale*masscalc.distance)**2)**0.5).to(u.pc,
                                                                                     u.dimensionless_angles())))
print("Total pixels > {3}: {0} = {1}; r_eff = {2}, fraction={4}, fraction of signal: {5}"
      .format(definitely_thick_if_20K.sum(), definitely_thick_if_20K.sum()/ppbeam,
              ((definitely_thick_if_20K.sum() *
                (pixel_scale*masscalc.distance)**2)**0.5).to(u.pc,
                                                             u.dimensionless_angles()),
              thick_threshold,
              definitely_thick_if_20K.sum()/np.isfinite(data).sum(),
              definitely_thick_if_20K.sum()/definitely_signal.sum(),
             )
     )
print("Total flux: {0}".format(total_signal))
print("Total mass(20K): {0}".format(total_signal * masscalc.mass_conversion_factor()*u.M_sun/u.Jy))
print("Threshold column (20K): {0:e}".format(threshold_column))
print("Threshold density (20K): {0:e}".format(threshold_density))
flux_of_filament = 132*u.Jy/u.beam/ppbeam
print("Total *filament* mass(20K): {0}".format(flux_of_filament * masscalc.mass_conversion_factor()*u.M_sun/u.Jy))
print("Total *filament* mass(100K): {0}".format(flux_of_filament * masscalc.mass_conversion_factor(TK=100)*u.M_sun/u.Jy))
print("Filament line mass: {0}".format(flux_of_filament *
                                       masscalc.mass_conversion_factor() *
                                       u.M_sun/u.Jy /
                                       (0.0027*u.deg*masscalc.distance).to(u.pc, u.dimensionless_angles())
                                      )
     )

print("Total mass(40K): {0}".format(total_signal * masscalc.mass_conversion_factor(TK=40)*u.M_sun/u.Jy))
threshold_column_40 = masscalc.dust.colofsnu(nu=masscalc.centerfreq,
                                             snu=threshold*u.beam,
                                             beamomega=beam.sr,
                                             temperature=40*u.K).to(u.cm**-2,
                                                                u.dimensionless_angles())
print("Threshold column (40K): {0:e}".format(threshold_column_40))
threshold_density_40 = (masscalc.mass_conversion_factor(TK=40) * (threshold*u.beam).to(u.mJy).value /
                        (4/3.*np.pi) /
                        (beam.sr.value*masscalc.distance**2)**(1.5) /
                        (2.8*constants.m_p)).to(1/u.cm**3)
print("Threshold density (40K): {0:e}".format(threshold_density_40))
flux_of_filament = 132*u.Jy/u.beam/ppbeam
print("Total *filament* mass(40K): {0}".format(flux_of_filament * masscalc.mass_conversion_factor(TK=40)*u.M_sun/u.Jy))
print("Total *filament* mass(100K): {0}".format(flux_of_filament * masscalc.mass_conversion_factor(TK=100)*u.M_sun/u.Jy))
print("Filament line mass: {0}".format(flux_of_filament *
                                       masscalc.mass_conversion_factor() *
                                       u.M_sun/u.Jy /
                                       (0.0027*u.deg*masscalc.distance).to(u.pc, u.dimensionless_angles())
                                      )
     )

dendro_merge = Table.read(paths.tpath('dendro_merge_continuum_and_line.ipac'), format='ascii.ipac')
corelike = dendro_merge['corelike'] == 'True'
dendro_protostars_pt2 = dendro_merge['cont_flux0p2arcsec'].sum()*u.Jy
print("Total protostar flux (0.2): {0}".format(dendro_protostars_pt2))
print("Total protostar flux (peak): {0}".format(dendro_merge['peak_cont_flux'].sum()*u.Jy))
total_minus_protostars = total_signal - dendro_protostars_pt2
print("Total recovered flux minus protostars: {0}".format(total_minus_protostars))
print("Total mass minus protostars (20K): {0}".format(total_minus_protostars * masscalc.mass_conversion_factor()*u.M_sun/u.Jy))

# determine BGPS total mass
regions = pyregion.open(paths.rpath("12m_pointings.reg"))
bgps_fh = fits.open("/Users/adam/work/w51/v2.0_ds2_l050_13pca_map20.fits")
mask = regions.get_mask(bgps_fh[0])
bgps_sum = bgps_fh[0].data[mask].sum() * u.Jy
bgps_ppbeam = bgps_fh[0].header['PPBEAM']
bgps_totalflux = bgps_sum/bgps_ppbeam
print("Total flux (BGPS, 271 GHz): {0}".format(bgps_totalflux))
bgps_scaled_225 = bgps_totalflux*(masscalc.centerfreq/(271.4*u.GHz))**3.5
bgps_scaled_225_4 = bgps_totalflux*(masscalc.centerfreq/(271.4*u.GHz))**4.0
bgps_scaled_225_3 = bgps_totalflux*(masscalc.centerfreq/(271.4*u.GHz))**3.0
print("Total flux (BGPS, 225 GHz, alpha=3.5): {0}".format(bgps_scaled_225))
bgps_totalmass = masscalc.dust.massofsnu(nu=271.4*u.GHz,
                                         snu=bgps_totalflux,
                                         distance=masscalc.distance)
print("Total mass (BGPS, 20K): {0}".format(bgps_totalmass))
print("*total* Fraction of recovered flux alpha=3.5: {0}".format(total_signal / bgps_scaled_225))
print("*total* Fraction of recovered flux alpha=3: {0}".format(total_signal / bgps_scaled_225_3))
print("*total* Fraction of recovered flux alpha=4: {0}".format(total_signal / bgps_scaled_225_4))
print("*nonprotostellar* Fraction of recovered flux alpha=3.5: {0}".format(total_minus_protostars / bgps_scaled_225))
print("*nonprotostellar* Fraction of recovered flux alpha=3: {0}".format(total_minus_protostars / bgps_scaled_225_3))
print("*nonprotostellar* Fraction of recovered flux alpha=4: {0}".format(total_minus_protostars / bgps_scaled_225_4))
print("Protostellar fraction: {0}".format(dendro_protostars_pt2 / bgps_scaled_225))

print()
print()


cores_merge = Table.read(paths.tpath('core_continuum_and_line.ipac'), format='ascii.ipac')

smallest_mass = masscalc.dust.massofsnu(nu=226*u.GHz,
                                        snu=cores_merge['peak'].min()*u.Jy,
                                        distance=masscalc.distance,
                                        beamomega=beam)
print("Faintest point source 20K mass: {0}"
      .format(smallest_mass))

smallest_density = (smallest_mass / (4/3.*np.pi) /
                    (beam.sr.value*masscalc.distance**2)**(1.5) /
                    (2.8*constants.m_p)).to(1/u.cm**3)

print("Faintest point source 20K density: {0}, log={1}"
      .format(smallest_density, np.log10(smallest_density.value)))



print()
print()

# this section repeats the analysis of dust_proerties, but it serves as a nice
# independent check since I did it 6 months later.... give or take...
beam_radius = ((beam.sr/(2*np.pi))**0.5 * masscalc.distance).to(u.pc,
                                                                u.dimensionless_angles())

# what luminosity is produced within the optically thick region?
e2epeak = np.nanmax(data[1300:1500,800:1000])
print("e2e peak flux: {0}".format(e2epeak))
e2etbpeak = e2epeak * beam.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("e2e peak brightness tem: {0}".format(e2etbpeak))
e2e_blackbody_luminosity = (constants.sigma_sb * e2etbpeak**4 *
                            (4*np.pi*beam_radius**2)).to(u.L_sun)
print("e2e blackbody luminosity: {0} = {0:e}".format(e2e_blackbody_luminosity,))

e8peak = np.nanmax(data[1200:1300,800:1000])
print("e8 peak flux: {0}".format(e8peak))
e8tbpeak = e8peak * beam.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("e8 peak brightness tem: {0}".format(e8tbpeak))
e8_blackbody_luminosity = (constants.sigma_sb * e8tbpeak**4 *
                           (4*np.pi*beam_radius**2)).to(u.L_sun)
print("e8 blackbody luminosity: {0} = {0:e}".format(e8_blackbody_luminosity,))


northpeak = np.nanmax(data[1900:2100,1900:2100])
print("north peak flux: {0}".format(northpeak))
northtbpeak = northpeak * beam.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("north peak brightness tem: {0}".format(northtbpeak))
north_blackbody_luminosity = (constants.sigma_sb * northtbpeak**4 *
                              (4*np.pi*beam_radius**2)).to(u.L_sun)
print("north blackbody luminosity: {0} = {0:e}".format(north_blackbody_luminosity,))

print()
print()



threecore_reg = pyregion.open(paths.rpath("three1ascores.reg"))
threecore_mask = threecore_reg.get_mask(contfile[0])
threecore_total = data[threecore_mask & (np.isfinite(data))].sum() / ppbeam
print("Total flux in the three 'main cores': {0}".format(threecore_total))
print("Fraction of total flux in the three 'main cores': {0}".format(threecore_total/total_signal))
print(" Fraction of recovered flux from BGPS alpha=3.5: {0}".format(threecore_total / bgps_scaled_225))
print(" Fraction of recovered flux from BGPS alpha=3: {0}".format(threecore_total / bgps_scaled_225_3))
print(" Fraction of recovered flux from BGPS alpha=4: {0}".format(threecore_total / bgps_scaled_225_4))
print("Area fraction in 'main cores': {0}".format(threecore_mask.sum()/mask.sum()))


planck_217 = fits.open('../../planckwmap/PLCKI_C290.925+14.509_217GHz.fits')
# https://wiki.cosmos.esa.int/planckpla/index.php/Effective_Beams
beam_planck = radio_beam.Beam(4.990*u.arcmin)
planck_217_flux = (planck_217[0].data*u.K).to(u.Jy, beam_planck.jtok_equiv(217*u.GHz))
pixel_area_planck = np.abs(planck_217[0].header['CDELT1'] * planck_217[0].header['CDELT2'])*u.deg**2
ppbeam_planck = (beam_planck.sr/pixel_area_planck).decompose()
planck_12mptg_mask = regions.get_mask(planck_217[0])
planck_12mptg_total = planck_217_flux[planck_12mptg_mask].sum() / ppbeam_planck
print("Total Planck flux in 12m ptg area: {0}".format(planck_12mptg_total))

whole_w51_reg = pyregion.open(paths.rpath('whole_w51_cloud.reg'))
whole_w51_mask_planck = whole_w51_reg.get_mask(planck_217[0])
planck_whole_total = planck_217_flux[whole_w51_mask_planck].sum() / ppbeam_planck
print("Total Planck flux in entire cloud: {0}".format(planck_whole_total))

whole_w51_mask_BGPS = whole_w51_reg.get_mask(bgps_fh[0])
whole_w51_bgps_sum = bgps_fh[0].data[whole_w51_mask_BGPS].sum() * u.Jy / bgps_ppbeam
whole_bgps_scaled_225 = whole_w51_bgps_sum*(masscalc.centerfreq/(271.4*u.GHz))**3.5
whole_bgps_scaled_225_4 = whole_w51_bgps_sum*(masscalc.centerfreq/(271.4*u.GHz))**4.0
whole_bgps_scaled_225_3 = whole_w51_bgps_sum*(masscalc.centerfreq/(271.4*u.GHz))**3.0
print("Whole BGPS flux total in W51: {0}".format(whole_w51_bgps_sum))
print("Whole BGPS flux total in W51, scaled with alpha=3.0 to 225: {0}".format(whole_bgps_scaled_225_3))
print("Whole BGPS flux total in W51, scaled with alpha=3.5 to 225: {0}".format(whole_bgps_scaled_225))
print("Whole BGPS flux total in W51, scaled with alpha=4.0 to 225: {0}".format(whole_bgps_scaled_225_4))


print()

kappa = masscalc.dust.kappa(masscalc.centerfreq)
mincol = (kappa * (2.8 * constants.m_p)).to(u.cm**2)**-1
minmass = ((beam.sr * masscalc.distance**2) / kappa).to(u.M_sun, u.dimensionless_angles())
print("minimum optically thick mass in a beam: {0}".format(minmass))
