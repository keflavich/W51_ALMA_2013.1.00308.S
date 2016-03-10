import numpy as np
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy import log
import paths
import pyregion
from astropy.io import fits
from astropy.table import Table,Column
import masscalc
import radio_beam
import astrodendro
from astropy import wcs

#contfile = fits.open(paths.dpath('selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_10mJy.image.pbcor.fits'))
#ln selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy.image.pbcor.fits W51_te_continuum_best.fits
#ln selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy.residual.fits W51_te_continuum_best_residual.fits
contfile = fits.open(paths.dpath('W51_te_continuum_best.fits'))
data = contfile[0].data
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
metadata['spatial_scale'] =  0.05 * u.arcsec
metadata['beam_major'] =  beam.major
metadata['beam_minor'] =  beam.minor
metadata['wavelength'] =  218.22219*u.GHz
metadata['velocity_scale'] = u.km/u.s
metadata['wcs'] = mywcs

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

cat_mask = (ppcat['is_leaf'] & (ppcat['peak_flux']>7*ppcat['noise']) &
            (ppcat['mean_flux']>5*ppcat['noise']) &
            (ppcat['min_flux']>2*ppcat['noise']))
pruned_ppcat = ppcat[cat_mask]
mask = dend.index_map.copy()
for ii in ProgressBar(list(range(len(ppcat)))):
    if ii not in pruned_ppcat['_idx']:
        mask[mask == ii] = -1
outf = fits.PrimaryHDU(data=mask, header=contfile[0].header)
outf.writeto('dendrograms_min1mJy_diff1mJy_mask_pruned.fits', clobber=True)
