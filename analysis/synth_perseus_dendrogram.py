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
from astropy import coordinates
from astropy import table
from FITS_tools import hcongrid

contfile_original = fits.open(paths.pspath('perseus_250_to_w51_2.image.fits'))
contfile = fits.open(paths.pspath('perseus_250_2_model_tclean_clean.image.fits'))
data_original = contfile_original[0].data
data = contfile[0].data

header_in = wcs.WCS(contfile_original[0].header).celestial.to_header()
header_in['NAXIS1'] = contfile_original[0].data.shape[-1]
header_in['NAXIS2'] = contfile_original[0].data.shape[-2]
reproj_fits_in = hcongrid.hastrom(contfile_original[0].data.squeeze(),
                                  header_in,
                                  contfile[0].header)
data_original = reproj_fits_in
data_original[np.isnan(data)] = np.nan

# estimate the noise from the local standard deviation of the residuals
residfile = fits.open(paths.pspath('perseus_250_2_model_tclean_clean.residual.fits'))
resid = residfile[0].data
smresid = convolve_fft(np.nan_to_num(resid), Gaussian2DKernel(30))
# have *low* noise outside when adding noise to the input image
synthnoise = convolve_fft(np.abs(resid-smresid),  Gaussian2DKernel(30))
resid[np.isnan(resid)] = 0.01 # make the noise outside very high
noise = convolve_fft(np.abs(resid-smresid),  Gaussian2DKernel(30))
residfile[0].data = noise
residfile.writeto(paths.pspath('perseus_250_2_model_tclean_clean_noise.fits'), clobber=True)
# lowest reasonable noise level is 0.2 mJy/beam

mywcs = wcs.WCS(contfile[0].header)

# Add noise to the data so that there is a sensible floor
noise_level = 0.0002
data_original += np.random.randn(*data_original.shape) * noise_level

orig_dend = astrodendro.Dendrogram.compute(data_original, min_value=0.001,
                                           min_delta=0.0004, min_npix=10,
                                           wcs=mywcs)
orig_dend.save_to('orig_perseus_dendrograms_min1mJy_diff1mJy.fits')

# ~5-sigma and ~2-sigma
# the labels have not kept up with the values...
dend = astrodendro.Dendrogram.compute(data, min_value=0.001, min_delta=0.0004,
                                      min_npix=10,
                                      wcs=mywcs)

dend.save_to('perseus_dendrograms_min1mJy_diff1mJy.fits')

df = fits.open('perseus_dendrograms_min1mJy_diff1mJy.fits')
df[1].header.update(df[0].header)
df[2].header.update(df[0].header)
df.writeto('perseus_dendrograms_min1mJy_diff1mJy.fits', clobber=True)

beam = radio_beam.Beam.from_fits_header(contfile[0].header)
metadata = {}
metadata['data_unit'] = u.Jy/u.beam
metadata['spatial_scale'] = 0.05 * u.arcsec
metadata['beam_major'] = beam.major
metadata['beam_minor'] = beam.minor
metadata['wavelength'] = 218.22219*u.GHz
metadata['velocity_scale'] = u.km/u.s
metadata['wcs'] = mywcs

orig_ppcat = astrodendro.pp_catalog(orig_dend, metadata)
ppcat = astrodendro.pp_catalog(dend, metadata)

# add a 'noise' column to the catalog
keys = ['noise','is_leaf','peak_flux','min_flux','mean_flux']
columns = {k:[] for k in (keys)}
for ii, row in enumerate(ProgressBar(ppcat)):
    structure = dend[row['_idx']]
    assert structure.idx == row['_idx'] == ii
    dend_inds = structure.indices()
    columns['noise'].append(noise[dend_inds].mean())
    columns['is_leaf'].append(structure.is_leaf)
    columns['peak_flux'].append(contfile[0].data[dend_inds].max())
    columns['min_flux'].append(contfile[0].data[dend_inds].min())
    columns['mean_flux'].append(contfile[0].data[dend_inds].mean())

for k in columns:
    if k not in ppcat.keys():
        ppcat.add_column(Column(name=k, data=columns[k]))

cat_mask = (ppcat['is_leaf'] & (ppcat['peak_flux']>8*ppcat['noise']) &
            (ppcat['mean_flux']>5*ppcat['noise']) &
            (ppcat['min_flux']>1*ppcat['noise']))
pruned_ppcat = ppcat[cat_mask]
mask = dend.index_map.copy()
for ii in ProgressBar(list(range(len(ppcat)))):
    if ii not in pruned_ppcat['_idx']:
        mask[mask == ii] = -1
outf = fits.PrimaryHDU(data=mask, header=contfile[0].header)
outf.writeto('perseus_dendrograms_min1mJy_diff1mJy_mask_pruned.fits', clobber=True)

pruned_ppcat.write(paths.tpath("perseus_dendrogram_continuum_catalog.ipac"), format='ascii.ipac')
orig_ppcat.write(paths.tpath("perseus_original_dendrogram_continuum_catalog.ipac"), format='ascii.ipac')

# Now mask the original data
keys = ['noise','is_leaf','peak_flux','min_flux','mean_flux']
columns = {k:[] for k in (keys)}
for ii, row in enumerate(ProgressBar(orig_ppcat)):
    structure = orig_dend[row['_idx']]
    assert structure.idx == row['_idx'] == ii
    dend_inds = structure.indices()
    # do I want the real noise or the "measured" noise here?
    columns['noise'].append(noise_level) # noise[dend_inds].mean())
    columns['is_leaf'].append(structure.is_leaf)
    columns['peak_flux'].append(data_original[dend_inds].max())
    columns['min_flux'].append(data_original[dend_inds].min())
    columns['mean_flux'].append(data_original[dend_inds].mean())

for k in columns:
    if k not in orig_ppcat.keys():
        orig_ppcat.add_column(Column(name=k, data=columns[k]))

orig_cat_mask = (orig_ppcat['is_leaf'] &
                 (orig_ppcat['peak_flux']>8*orig_ppcat['noise']) &
                 (orig_ppcat['mean_flux']>5*orig_ppcat['noise']) &
                 (orig_ppcat['min_flux']>1*orig_ppcat['noise']))
pruned_orig_ppcat = orig_ppcat[orig_cat_mask]
orig_mask = orig_dend.index_map.copy()
for ii in ProgressBar(list(range(len(orig_ppcat)))):
    if ii not in pruned_orig_ppcat['_idx']:
        orig_mask[orig_mask == ii] = -1
outf = fits.PrimaryHDU(data=orig_mask, header=contfile[0].header)
outf.writeto('original_perseus_dendrograms_mask_pruned.fits', clobber=True)

c_obs = coordinates.SkyCoord(pruned_ppcat['x_cen'], pruned_ppcat['y_cen'], frame='fk5',
                             unit=(u.deg,u.deg))
c_inp = coordinates.SkyCoord(pruned_orig_ppcat['x_cen'], pruned_orig_ppcat['y_cen'], frame='fk5', unit=(u.deg,u.deg))
orig_to_obs = dict(zip(('inds','sep','_'), c_obs.match_to_catalog_sky(c_inp)))
obs_to_orig = dict(zip(('inds','sep','_'), c_inp.match_to_catalog_sky(c_obs)))

print("Observations match {0} out of {1} of the original leaves to 1 arcsec"
      .format((obs_to_orig['sep'] < 1*u.arcsec).sum(), len(pruned_orig_ppcat)))
print("Of the observations, {0} out of {1} have <1 arcsec matches to original sources"
      .format((obs_to_orig['sep'] < 1*u.arcsec).sum(), len(pruned_ppcat)))
print("Input has {0} out of {1} matches to the observed leaves within 1 arcsec"
      .format((orig_to_obs['sep'] < 1*u.arcsec).sum(), len(pruned_orig_ppcat)))

pruned_ppcat.add_column(Column(data=orig_to_obs['inds'],
                               name='match_inds'))
pruned_ppcat.add_column(Column(data=orig_to_obs['sep'],
                               name='match_separation'))
temporary = pruned_orig_ppcat[pruned_ppcat['match_inds']]
merged_orig_onto_obs = table.hstack([pruned_ppcat, temporary])

# ds9 -multiframe ../analysis/*perseus*.fits ../perseus_synth/perseus_250_2_model_tclean_clean_noise.fits -lock frames image -frame 2 -scale minmax -cmap sls -frame 3 -frame delete -frame 7 -frame delete -frame 6 -cmap sls -frame 8 -cmap sls -frame 4 -cmap sls -frame 5 -cmap value 8.5 0.05 -frame 9 -cmap value 12 0.03 -frame 1 -cmap value 8.5 0.05 -lock crosshairs image  &
