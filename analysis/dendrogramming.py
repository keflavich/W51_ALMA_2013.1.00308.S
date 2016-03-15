import numpy as np
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.utils.console import ProgressBar
from astropy import units as u
import paths
from astropy.io import fits
from astropy.table import Column
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
data = contfile[0].data

# estimate the noise from the local standard deviation of the residuals
residfile = fits.open(paths.dpath('W51_te_continuum_best_residual.fits'))
resid = residfile[0].data
smresid = convolve_fft(np.nan_to_num(resid), Gaussian2DKernel(30))
resid[np.isnan(resid)] = 0.01 # make the noise outside very high
noise = convolve_fft(np.abs(resid-smresid),  Gaussian2DKernel(30))
residfile[0].data = noise
residfile.writeto(paths.dpath('W51_te_continuum_best_noise.fits'), clobber=True)
# lowest reasonable noise level is 0.2 mJy/beam

mywcs = wcs.WCS(contfile[0].header)

# ~5-sigma and ~2-sigma
# the labels have not kept up with the values...
dend = astrodendro.Dendrogram.compute(data, min_value=0.001, min_delta=0.0004,
                                      min_npix=10,
                                      wcs=mywcs)

dend.save_to('dendrograms_min1mJy_diff1mJy.fits')

df = fits.open('dendrograms_min1mJy_diff1mJy.fits')
df[1].header.update(df[0].header)
df[2].header.update(df[0].header)
df.writeto('dendrograms_min1mJy_diff1mJy.fits', clobber=True)

beam = radio_beam.Beam.from_fits_header(contfile[0].header)
metadata = {}
metadata['data_unit'] = u.Jy/u.beam
pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
pixel_scale_as = pixel_scale.to(u.arcsec).value
metadata['spatial_scale'] = pixel_scale
metadata['beam_major'] = beam.major
metadata['beam_minor'] = beam.minor
metadata['wavelength'] = 218.22219*u.GHz
metadata['velocity_scale'] = u.km/u.s
metadata['wcs'] = mywcs
ppbeam = (beam.sr/(pixel_scale**2)).decompose().value

ppcat = astrodendro.pp_catalog(dend, metadata)

# add a 'noise' column to the catalog
keys = ['noise', 'is_leaf', 'peak_cont_flux', 'min_flux', 'mean_flux', 'peak_cont_mass',
        'peak_cont_col', 'beam_area']
radii = (0.2,0.4,0.6,0.8,1.0,1.5)*u.arcsec
columns = {k:[] for k in (keys)}
for ii, row in enumerate(ProgressBar(ppcat)):
    structure = dend[row['_idx']]
    assert structure.idx == row['_idx'] == ii
    dend_inds = structure.indices()
    columns['noise'].append(noise[dend_inds].mean())
    columns['is_leaf'].append(structure.is_leaf)
    peakflux = data[dend_inds].max()
    columns['peak_cont_flux'].append(peakflux)
    columns['min_cont_flux'].append(data[dend_inds].min())
    columns['mean_cont_flux'].append(data[dend_inds].mean())
    columns['peak_cont_mass'].append(masscalc.mass_conversion_factor()*peakflux)
    columns['peak_cont_col'].append(masscalc.col_conversion_factor()*peakflux)
    columns['beam_area'].append(beam.sr.value)
for k in columns:
    if k not in ppcat.keys():
        ppcat.add_column(Column(name=k, data=columns[k]))

cat_mask = (ppcat['is_leaf'] & (ppcat['peak_cont_flux']>8*ppcat['noise']) &
            (ppcat['mean_cont_flux']>5*ppcat['noise']) &
            (ppcat['min_cont']>1*ppcat['noise']))
pruned_ppcat = ppcat[cat_mask]
mask = dend.index_map.copy()
for ii in ProgressBar(list(range(len(ppcat)))):
    if ii not in pruned_ppcat['_idx']:
        mask[mask == ii] = -1
outf = fits.PrimaryHDU(data=mask, header=contfile[0].header)
outf.writeto('dendrograms_min1mJy_diff1mJy_mask_pruned.fits', clobber=True)

keys = ['flux{0}arcsec'.format(rr.value) for rr in radii]
columns = {k:[] for k in (keys)}
for ii, row in enumerate(ProgressBar(pruned_ppcat)):
    size = max(radii)*2.2
    xc,yc = row['x_cen'], row['y_cen']
    position = coordinates.SkyCoord(xc, yc, frame='fk5', unit=(u.deg,u.deg))
    cutout = Cutout2D(data, position, size, mywcs, mode='partial')
    for rr in radii:
        aperture = photutils.SkyCircularAperture(positions=position,
                                                 r=rr).to_pixel(cutout.wcs)
        aperture_data = photutils.aperture_photometry(data=cutout.data,
                                                      apertures=aperture,
                                                      method='exact')
        flux_jybeam = aperture_data[0]['aperture_sum']
        #flux_jybeam_average = flux_jybeam / aperture.area()
        #flux_jysr_average = flux_jybeam_average / beam.sr.value
        #flux_jysr = flux_jysr_average * aperture.area()
        #flux_jy = (flux_jysr * (pixel_scale**2).to(u.sr).value)
        columns['cont_flux{0}arcsec'.format(rr.value)].append(flux_jybeam/ppbeam)

for k in columns:
    if k not in pruned_ppcat.keys():
        pruned_ppcat.add_column(Column(name=k, data=columns[k]))

for k in pruned_ppcat.colnames:
    if "." in k:
        pruned_ppcat.rename_column(k, k.replace(".","p"))

pruned_ppcat.meta = {'ppbeam': ppbeam,
                     'beam_area_sr': beam.sr.value,
                     'pixel_scale_as': pixel_scale_as}

annulus_mean = ((pruned_ppcat['cont_flux0p4arcsec']-pruned_ppcat['cont_flux0p2arcsec']) /
                (np.pi*(0.4-0.2)**2/pixel_scale_as**2) * ppbeam)
core_ish = pruned_ppcat['peak_cont_flux'] > annulus_mean
pruned_ppcat.add_column(Column(data=core_ish, name='corelike'))

pruned_ppcat.write(paths.tpath("dendrogram_continuum_catalog.ipac"),
                   format='ascii.ipac')
