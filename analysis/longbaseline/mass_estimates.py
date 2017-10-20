
import pylab as pl
import os
import numpy as np
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.utils.console import ProgressBar
from astropy import constants
from astropy import units as u
import astropy.visualization
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
import regions

import warnings
warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)
warnings.filterwarnings('ignore', category=RuntimeWarning,
                        message="invalid value encountered in greater", append=True)


lores_contfile = fits.open(paths.dpath('W51_te_continuum_best.fits'))
datalores = u.Quantity(lores_contfile[0].data.squeeze(), unit=lores_contfile[0].header['BUNIT'])

beam_lores = radio_beam.Beam.from_fits_header(lores_contfile[0].header)

#contfile = fits.open(paths.dpath('selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_10mJy.image.pbcor.fits'))
#ln selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy.image.pbcor.fits W51_te_continuum_best.fits
#ln selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy.residual.fits W51_te_continuum_best_residual.fits
#contfile_e2e8 = fits.open(paths.dpath('longbaseline/W51e2cax.cont.image.pbcor.fits'))
contfile_e2e8 = fits.open(paths.dpath('longbaseline/W51e2_cont_briggsSC_tclean.image.fits'))
datae2e8 = u.Quantity(contfile_e2e8[0].data.squeeze(), unit=contfile_e2e8[0].header['BUNIT'])

beam_e2e8 = radio_beam.Beam.from_fits_header(contfile_e2e8[0].header)

wcs_e2e8 = wcs.WCS(contfile_e2e8[0].header)

# this section repeats the analysis of dust_proerties, but it serves as a nice
# independent check since I did it 6 months later.... give or take...
beam_radius_e2e8 = ((beam_e2e8.sr/(2*np.pi))**0.5 * masscalc.distance).to(u.pc,
                                                                          u.dimensionless_angles())

#contfile_north = fits.open(paths.dpath('longbaseline/W51ncax.cont.image.pbcor.fits'))
contfile_north = fits.open(paths.dpath('longbaseline/W51n_cont_briggsSC_tclean.image.fits'))
datanorth = u.Quantity(contfile_north[0].data.squeeze(), unit=contfile_north[0].header['BUNIT'])

beam_north = radio_beam.Beam.from_fits_header(contfile_north[0].header)

wcs_north = wcs.WCS(contfile_north[0].header)

# this section repeats the analysis of dust_proerties, but it serves as a nice
# independent check since I did it 6 months later.... give or take...
beam_radius_north = ((beam_north.sr/(2*np.pi))**0.5 *
                     masscalc.distance).to(u.pc, u.dimensionless_angles())

e2 = regions.CircleSkyRegion(coordinates.SkyCoord('19:23:43.969 +14:30:34.518', frame='icrs', unit=(u.hour, u.deg)),
                             radius=0.616*u.arcsec)
e2pix = e2.to_pixel(wcs_e2e8)
e2mask = e2pix.to_mask()

e8 = regions.CircleSkyRegion(coordinates.SkyCoord('19:23:43.907 +14:30:28.267', frame='icrs', unit=(u.hour, u.deg)),
                             radius=0.464*u.arcsec)
e8pix = e8.to_pixel(wcs_e2e8)
e8mask = e8pix.to_mask()

north = regions.CircleSkyRegion(coordinates.SkyCoord('19:23:40.054 +14:31:05.513', frame='icrs', unit=(u.hour, u.deg)),
                                radius=0.412*u.arcsec)
northpix = north.to_pixel(wcs_north)
northmask = northpix.to_mask()

def masked_where(condition, data):
    data = data.copy()
    data[condition] = np.nan
    return data

cutout_e2 = masked_where(e2mask.data==0, e2mask.cutout(datae2e8))
cutout_e8 = masked_where(e8mask.data==0, e8mask.cutout(datae2e8))
#cutout_e8 = datae2e8[1200:1400,2500:2700]
#cutout_north = datanorth[2400:2600,2400:2600]
cutout_north = masked_where(northmask.data==0, northmask.cutout(datanorth))

e2mask_lores = e2.to_pixel(wcs.WCS(lores_contfile[0].header)).to_mask()
e8mask_lores = e8.to_pixel(wcs.WCS(lores_contfile[0].header)).to_mask()
northmask_lores = north.to_pixel(wcs.WCS(lores_contfile[0].header)).to_mask()

cutout_e2_lores = masked_where(e2mask_lores.data==0, e2mask_lores.cutout(datalores))
cutout_e8_lores = masked_where(e8mask_lores.data==0, e8mask_lores.cutout(datalores))
#cutout_e8 = datae2e8[1200:1400,2500:2700]
#cutout_north = datanorth[2400:2600,2400:2600]
cutout_north_lores = masked_where(northmask_lores.data==0, northmask_lores.cutout(datalores))


# what luminosity is produced within the optically thick region?
e2epeak = np.nanmax(cutout_e2)
print("e2e peak flux: {0}".format(e2epeak))
e2etbpeak = e2epeak * beam_e2e8.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("e2e peak brightness tem: {0}".format(e2etbpeak))
e2e_blackbody_luminosity = (constants.sigma_sb * e2etbpeak**4 *
                            (4*np.pi*beam_radius_e2e8**2)).to(u.L_sun)
print("e2e blackbody luminosity: {0} = {0:e} = 10^{1}".format(e2e_blackbody_luminosity, np.log10(e2e_blackbody_luminosity.value)))

e8peak = np.nanmax(cutout_e8)
print("e8 peak flux: {0}".format(e8peak))
e8tbpeak = e8peak * beam_e2e8.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("e8 peak brightness tem: {0}".format(e8tbpeak))
e8_blackbody_luminosity = (constants.sigma_sb * e8tbpeak**4 *
                           (4*np.pi*beam_radius_e2e8**2)).to(u.L_sun)
print("e8 blackbody luminosity: {0} = {0:e} = 10^{1}".format(e8_blackbody_luminosity, np.log10(e8_blackbody_luminosity.value)))


northpeak = np.nanmax(cutout_north)
print("north peak flux: {0}".format(northpeak))
northtbpeak = northpeak * beam_north.jtok(masscalc.centerfreq) * u.beam/u.Jy
print("north peak brightness tem: {0}".format(northtbpeak))
north_blackbody_luminosity = (constants.sigma_sb * northtbpeak**4 *
                              (4*np.pi*beam_radius_north**2)).to(u.L_sun)
print("north blackbody luminosity: {0} = {0:e} = 10^{1}".format(north_blackbody_luminosity, np.log10(north_blackbody_luminosity.value)))
print()


# compute integrated intensities
pixarea_north = wcs.utils.proj_plane_pixel_area(wcs.WCS(contfile_north[0].header)) * u.deg**2
pixarea_e2e8 = wcs.utils.proj_plane_pixel_area(wcs.WCS(contfile_e2e8[0].header)) * u.deg**2
ppbeam_north = (beam_north / pixarea_north).decompose().value
ppbeam_e2e8 = (beam_e2e8 / pixarea_e2e8).decompose().value
pixarea_lores = wcs.utils.proj_plane_pixel_area(wcs.WCS(lores_contfile[0].header)) * u.deg**2
ppbeam_lores = (beam_lores / pixarea_lores).decompose().value

integ_e2 = (cutout_e2[cutout_e2 > 3*u.mJy/u.beam]).sum() * u.beam / ppbeam_e2e8
integ_e8 = (cutout_e8[cutout_e8 > 3*u.mJy/u.beam]).sum() * u.beam / ppbeam_e2e8
integ_north = (cutout_north[cutout_north > 3*u.mJy/u.beam]).sum() * u.beam / ppbeam_north

integ_e2_lores = (cutout_e2_lores[cutout_e2_lores > 3*u.mJy/u.beam]).sum() * u.beam / ppbeam_lores
integ_e8_lores = (cutout_e8_lores[cutout_e8_lores > 3*u.mJy/u.beam]).sum() * u.beam / ppbeam_lores
integ_north_lores = (cutout_north_lores[cutout_north_lores > 3*u.mJy/u.beam]).sum() * u.beam / ppbeam_lores

print("Integrated intensity of    e2: {0:0.2f}  /  {1:0.2f} -> {2:0.1f}% recovered".format(integ_e2, integ_e2_lores, integ_e2/integ_e2_lores*100))
print("Integrated intensity of    e8: {0:0.2f}  /  {1:0.2f} -> {2:0.1f}% recovered".format(integ_e8, integ_e8_lores, integ_e8/integ_e8_lores*100))
print("Integrated intensity of north: {0:0.2f}  /  {1:0.2f} -> {2:0.1f}% recovered".format(integ_north, integ_north_lores, integ_north/integ_north_lores*100))
print()

integmass_e2 = integ_e2 * masscalc.mass_conversion_factor(e2etbpeak)/u.Jy
integmass_e8 = integ_e8 * masscalc.mass_conversion_factor(e8tbpeak)/u.Jy
integmass_north = integ_north * masscalc.mass_conversion_factor(northtbpeak)/u.Jy

print("Integrated mass of e2: {0}".format(integmass_e2))
print("Integrated mass of e8: {0}".format(integmass_e8))
print("Integrated mass of north: {0}".format(integmass_north))
print()


integmass_e2 = integ_e2 * masscalc.mass_conversion_factor(200*u.K)/u.Jy
integmass_e8 = integ_e8 * masscalc.mass_conversion_factor(200*u.K)/u.Jy
integmass_north = integ_north * masscalc.mass_conversion_factor(200*u.K)/u.Jy

print("Integrated mass of e2 at T=200K: {0}".format(integmass_e2))
print("Integrated mass of e8 at T=200K: {0}".format(integmass_e8))
print("Integrated mass of north at T=200K: {0}".format(integmass_north))


# use dendrograms to identify structures and measure masses
dende2 = astrodendro.Dendrogram.compute(cutout_e2.value, min_value=0.0005,
                                        min_delta=0.0005, wcs=wcs_e2e8.celestial)

wcs_e2cutout = wcs_e2e8.celestial[e2mask.bbox.iymin:, e2mask.bbox.ixmin:]

e2cutoutpix = e2.to_pixel(wcs_e2cutout)
struct = dende2.structure_at((int(e2cutoutpix.center.x), int(e2cutoutpix.center.y)))
trunk = struct.ancestor

e2inner = struct.values().sum() * u.Jy / ppbeam_e2e8
e2total = trunk.values().sum() * u.Jy / ppbeam_e2e8
e2filaments = e2total - e2inner

e2innermass = e2inner * masscalc.mass_conversion_factor(200*u.K) / u.Jy
e2totalmass = e2total * masscalc.mass_conversion_factor(200*u.K) / u.Jy
e2filamentsmass = e2filaments * masscalc.mass_conversion_factor(200*u.K) / u.Jy

e2innermass_pktb = e2inner * masscalc.mass_conversion_factor(e2etbpeak) / u.Jy
e2totalmass_pktb = e2total * masscalc.mass_conversion_factor(e2etbpeak) / u.Jy
e2filamentsmass_pktb = e2filaments * masscalc.mass_conversion_factor(e2etbpeak) / u.Jy

nbeams_inner = struct.get_npix() / ppbeam_e2e8

mass_per_beam_thick = ((beam_e2e8.sr*masscalc.distance**2).to(u.cm**2,
                                                              u.dimensionless_angles())
                       * (1/0.0083*u.cm**-2*u.g / (2.8*u.Da)).to(u.cm**-2) *
                       2.8*u.Da).to(u.M_sun)

print()
print("e2 nbeams at inner={0}, total={1}, filaments={2}".format(nbeams_inner,
                                                                trunk.get_npix() / ppbeam_e2e8,
                                                                (trunk.get_npix() - struct.get_npix())/ppbeam_e2e8))
print("optically thick mass of inner: {0}".format(nbeams_inner * mass_per_beam_thick))
print("Integrated intensity of inner e2 = {0:0.2f}, e2 total = {1:0.2f}, e2 filaments = {2:0.2f}"
      .format(e2inner, e2total, e2filaments))
print("Mass at 200K of inner e2 = {0:0.2f}, e2 total = {1:0.2f}, e2 filaments = {2:0.2f}"
      .format(e2innermass, e2totalmass, e2filamentsmass))
print("Mass at {3:0.1f} of inner e2 = {0:0.2f}, e2 total = {1:0.2f}, e2 filaments = {2:0.2f}"
      .format(e2innermass_pktb, e2totalmass_pktb, e2filamentsmass_pktb, e2etbpeak))


fig1 = pl.figure(1)
fig1.clf()
ax = pl.axes(projection=wcs_e2cutout)
fig1.add_axes(ax)
norm = astropy.visualization.ImageNormalize(cutout_e2, stretch=astropy.visualization.AsinhStretch())
ax.imshow(cutout_e2, origin='lower', interpolation='none', cmap='gray', norm=norm)
ax.contour(struct.get_mask(), levels=[0.5], colors=['orange'], origin='lower')
ax.contour(trunk.get_mask(), levels=[0.5], colors=['red'], origin='lower')
ax.set_xlabel("Right Ascension")
ax.set_ylabel("Declination")
fig1.savefig(paths.fpath("longbaseline/e2cutout_contours_for_masscalc.pdf"), bbox_inches='tight')


dende8 = astrodendro.Dendrogram.compute(cutout_e8.value, min_value=0.0005,
                                        min_delta=0.00001, wcs=wcs_e2e8.celestial)

e8filament = regions.CircleSkyRegion(coordinates.SkyCoord('19:23:43.899 +14:30:28.382',
                                                          frame='icrs', unit=(u.hour, u.deg)),
                                     radius=0.001*u.arcsec)
wcs_e8cutout = wcs_e2e8.celestial[e8mask.bbox.iymin:, e8mask.bbox.ixmin:]
e8cutoutpix = e8filament.to_pixel(wcs_e8cutout)
struct = dende8.structure_at((int(e8cutoutpix.center.y), int(e8cutoutpix.center.x)))
assert struct.idx == 63
print('e8 struct: {0}'.format(struct))
trunk = struct.ancestor

e8inner = struct.values().sum() * u.Jy / ppbeam_e2e8
e8total = trunk.values().sum() * u.Jy / ppbeam_e2e8
e8filaments = e8total - e8inner

e8innermass = e8inner * masscalc.mass_conversion_factor(200*u.K) / u.Jy
e8totalmass = e8total * masscalc.mass_conversion_factor(200*u.K) / u.Jy
e8filamentsmass = e8filaments * masscalc.mass_conversion_factor(200*u.K) / u.Jy

e8innermass_pktb = e8inner * masscalc.mass_conversion_factor(e8tbpeak) / u.Jy
e8totalmass_pktb = e8total * masscalc.mass_conversion_factor(e8tbpeak) / u.Jy
e8filamentsmass_pktb = e8filaments * masscalc.mass_conversion_factor(e8tbpeak) / u.Jy


print()
print("Integrated intensity of inner e8 = {0:0.2f}, e8 total = {1:0.2f}, e8 filaments = {2:0.2f}"
      .format(e8inner, e8total, e8filaments))
print("Mass at 200K of inner e8 = {0:0.2f}, e8 total = {1:0.2f}, e8 filaments = {2:0.2f}"
      .format(e8innermass, e8totalmass, e8filamentsmass))
print("Mass at {3:0.1f} of inner e8 = {0:0.2f}, e8 total = {1:0.2f}, e8 filaments = {2:0.2f}"
      .format(e8innermass_pktb, e8totalmass_pktb, e8filamentsmass_pktb, e8tbpeak))


fig1 = pl.figure(1)
fig1.clf()
ax = pl.axes(projection=wcs_e8cutout)
fig1.add_axes(ax)
norm = astropy.visualization.ImageNormalize(cutout_e8, stretch=astropy.visualization.AsinhStretch())
ax.imshow(cutout_e8, origin='lower', interpolation='none', cmap='gray', norm=norm)
ax.contour(struct.get_mask(), levels=[0.5], colors=['orange'], origin='lower')
ax.contour(trunk.get_mask(), levels=[0.5], colors=['red'], origin='lower')
ax.set_xlabel("Right Ascension")
ax.set_ylabel("Declination")
fig1.savefig(paths.fpath("longbaseline/e8cutout_contours_for_masscalc.pdf"), bbox_inches='tight')







dendnorth = astrodendro.Dendrogram.compute(cutout_north.value,
                                           min_value=0.0005, min_delta=0.0001,
                                           wcs=wcs_north.celestial)

northneighbor = regions.CircleSkyRegion(coordinates.SkyCoord(290.916939905898*u.deg, 14.51817652776441*u.deg,
                                                             frame='icrs'),
                                        radius=0.001*u.arcsec)
northcutoutpix = northneighbor.to_pixel(wcs_north.celestial[northmask.bbox.iymin:, northmask.bbox.ixmin:])
print('north cutout pix should be approximately x=73 y=118',northcutoutpix)
nstruct = dendnorth.structure_at((int(northcutoutpix.center.y), int(northcutoutpix.center.x)))
not_neighbor = [x for x in nstruct.parent.children if x is not nstruct][0]
struct = not_neighbor

northinner = struct.values().sum() * u.Jy / ppbeam_north
#northtotal = trunk.values().sum() * u.Jy / ppbeam_north
#northfilaments = northtotal - northinner

northinnermass = northinner * masscalc.mass_conversion_factor(200*u.K) / u.Jy
#northtotalmass = northtotal * masscalc.mass_conversion_factor(200*u.K) / u.Jy
#northfilamentsmass = northfilaments * masscalc.mass_conversion_factor(200*u.K) / u.Jy

northinnermass_pktb = northinner * masscalc.mass_conversion_factor(northtbpeak) / u.Jy
#northtotalmass_pktb = northtotal * masscalc.mass_conversion_factor(northtbpeak) / u.Jy
#northfilamentsmass_pktb = northfilaments * masscalc.mass_conversion_factor(northtbpeak) / u.Jy


print()
print("Integrated intensity of north = {0:0.2f}"
      .format(northinner, ))
print("Mass at 200K of north = {0:0.2f}"
      .format(northinnermass))
print("Mass at {1:0.1f} of north = {0:0.2f}"
      .format(northinnermass_pktb, northtbpeak))
