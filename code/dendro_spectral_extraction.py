import numpy as np
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.table import Table

tmplt = "full_W51_spw{0}_lines.fits"

pruned_ppcat = Table.read("dendrogram_continuum_catalog.ipac",
                          format='ascii.ipac')
dendromask = fits.getdata('dendrograms_min1mJy_diff1mJy_mask_pruned.fits')

for spw in (0,1,2,3):
    cube = SpectralCube.read(tmplt.format(spw))
    print(cube)

    for row in pruned_ppcat:
        name = row['_idx']
        print("Extracting {0} from {1}".format(name, spw))

        mask = dendromask == name
        dend_inds = np.where(mask)

        view = (slice(None), # all spectral channels
                slice(dend_inds[0].min(), dend_inds[0].max()+1),
                slice(dend_inds[1].min(), dend_inds[1].max()+1),
               )
        sc = cube[view].with_mask(mask[view[1:]])
        spec = sc.mean(axis=(1,2))
        spec.hdu.writeto("spectra/dendro{0:03d}_spw{1}_mean.fits".format(name, spw),
                         clobber=True)
