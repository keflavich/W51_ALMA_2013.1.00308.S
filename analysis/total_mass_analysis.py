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
total_signal = data[definitely_signal].sum() / ppbeam
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
fh = fits.open("/Users/adam/work/w51/v2.0_ds2_l050_13pca_map20.fits")
mask = regions.get_mask(fh[0])
bgps_sum = fh[0].data[mask].sum()
bgps_ppbeam = fh[0].header['PPBEAM']
bgps_totalflux = bgps_sum/bgps_ppbeam
print("Total flux (BGPS, 271 GHz): {0}".format(bgps_totalflux))
bgps_scaled_225 = bgps_totalflux*(masscalc.centerfreq/(271.4*u.GHz))**3.5
bgps_scaled_225_4 = bgps_totalflux*(masscalc.centerfreq/(271.4*u.GHz))**4.0
bgps_scaled_225_3 = bgps_totalflux*(masscalc.centerfreq/(271.4*u.GHz))**3.0
print("Total flux (BGPS, 225 GHz, alpha=3.5): {0}".format(bgps_scaled_225))
bgps_totalmass = masscalc.dust.massofsnu(nu=271.4*u.GHz,
                                         snu=bgps_totalflux*u.Jy,
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


# what luminosity is produced within the optically thick region?
e2epeak = np.nanmax(data[1300:1500,800:1000])
print("e2e peak flux: {0}".format(e2epeak))
e2etbpeak = e2epeak * beam.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("e2e peak brightness tem: {0}".format(e2etbpeak))
e2eradius = (beam.sr**0.5 * masscalc.distance).to(u.pc, u.dimensionless_angles())
e2e_blackbody_luminosity = (constants.sigma_sb * e2etbpeak**4 *
                            (4*np.pi*e2eradius**2)).to(u.L_sun)
print("e2e blackbody luminosity: {0} = {0:e}".format(e2e_blackbody_luminosity,))

e8peak = np.nanmax(data[1200:1300,800:1000])
print("e8 peak flux: {0}".format(e8peak))
e8tbpeak = e8peak * beam.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("e8 peak brightness tem: {0}".format(e8tbpeak))
e8radius = (beam.sr**0.5 * masscalc.distance).to(u.pc, u.dimensionless_angles())
e8_blackbody_luminosity = (constants.sigma_sb * e8tbpeak**4 *
                           (4*np.pi*e8radius**2)).to(u.L_sun)
print("e8 blackbody luminosity: {0} = {0:e}".format(e8_blackbody_luminosity,))


northpeak = np.nanmax(data[1900:2100,1900:2100])
print("north peak flux: {0}".format(northpeak))
northtbpeak = northpeak * beam.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("north peak brightness tem: {0}".format(northtbpeak))
northradius = (beam.sr**0.5 * masscalc.distance).to(u.pc, u.dimensionless_angles())
north_blackbody_luminosity = (constants.sigma_sb * northtbpeak**4 *
                              (4*np.pi*northradius**2)).to(u.L_sun)
print("north blackbody luminosity: {0} = {0:e}".format(north_blackbody_luminosity,))
