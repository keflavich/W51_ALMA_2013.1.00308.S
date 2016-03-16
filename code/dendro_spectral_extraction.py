import os
import numpy as np
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.table import Table
import FITS_tools
import pyregion

if os.path.exists('full_W51_7m12m_spw0_lines.fits'):
    tmplt = "full_W51_7m12m_spw{0}_lines.fits"
    suffix = "_7m12m"
else:
    tmplt = "full_W51_spw{0}_lines.fits"
    suffix = ''

pruned_ppcat = Table.read("dendrogram_continuum_catalog.ipac",
                          format='ascii.ipac')
#dendromask = fits.open('dendrograms_min1mJy_diff1mJy_mask_pruned.fits')

for spw in (0,1,2,3):
    cube = SpectralCube.read(tmplt.format(spw))
    print(cube)

    #mask_reproj = FITS_tools.hcongrid.hastrom(dendromask[0].data.astype('float'),
    #                                          dendromask[0].header,
    #                                          FITS_tools.strip_headers.flatten_header(cube.header))
    #rmask = np.round(mask_reproj.data).astype('int')

    for row in pruned_ppcat:
        name = row['_idx']
        print("Extracting {0} from {1}".format(name, spw))
        SL = pyregion.parse("fk5; circle({0},{1},0.5\")"
                            .format(row['x_cen'], row['y_cen']))

        #mask = rmask == name
        #dend_inds = np.where(mask)

        #view = (slice(None), # all spectral channels
        #        slice(dend_inds[0].min(), dend_inds[0].max()+1),
        #        slice(dend_inds[1].min(), dend_inds[1].max()+1),
        #       )
        #sc = cube[view].with_mask(mask[view[1:]])
        sc = cube.subcube_from_ds9region(SL)
        spec = sc.mean(axis=(1,2))
        spec.hdu.writeto("spectra/dendro{0:03d}_spw{1}_mean{2}.fits".format(name, spw, suffix),
                         clobber=True)

        bgSL = pyregion.parse("fk5; circle({0},{1},1.0\")"
                              .format(row['x_cen'], row['y_cen']))
        bgsc = cube.subcube_from_ds9region(bgSL)
        npix = np.count_nonzero(np.isfinite(bgsc[0,:,:]))
        bgspec = (bgsc.sum(axis=(1,2)) - sc.sum(axis=(1,2))) / npix
        bgspec.hdu.writeto("spectra/dendro{0:03d}_spw{1}_background_mean{2}.fits".format(name, spw, suffix),
                           clobber=True)
