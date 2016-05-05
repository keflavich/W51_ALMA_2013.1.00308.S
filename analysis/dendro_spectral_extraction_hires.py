import os
import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.table import Table
from astropy import coordinates
import radio_beam
import pyregion

if os.path.exists('full_W51_7m12m_spw1_hires_lines.fits'):
    tmplt = "full_W51_7m12m_spw{0}_hires_lines.fits"
    suffix = "_7m12m"
else:
    raise ValueError("File doesn't exist")
    tmplt = "full_W51_spw{0}_lines.fits"
    suffix = ''

pruned_ppcat = Table.read("dendrogram_continuum_catalog.ipac",
                          format='ascii.ipac')
#dendromask = fits.open('dendrograms_min1mJy_diff1mJy_mask_pruned.fits')

for spw in (1,3,2,0):
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
        coord = coordinates.SkyCoord(row['x_cen'], row['y_cen'], frame='fk5',
                                     unit=(u.deg, u.deg))

        #mask = rmask == name
        #dend_inds = np.where(mask)

        #view = (slice(None), # all spectral channels
        #        slice(dend_inds[0].min(), dend_inds[0].max()+1),
        #        slice(dend_inds[1].min(), dend_inds[1].max()+1),
        #       )
        #sc = cube[view].with_mask(mask[view[1:]])
        sc = cube.subcube_from_ds9region(SL)
        spec = sc.mean(axis=(1,2))
        spec.meta['beam'] = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in spec.beams]),
                                            minor=np.nanmedian([bm.minor.to(u.deg).value for bm in spec.beams]),
                                            pa=np.nanmedian([bm.pa.to(u.deg).value for bm in spec.beams]),
                                           )
        spec.hdu.writeto("spectra/dendro{0:03d}_spw{1}_hires_0.5as_mean{2}.fits".format(name, spw, suffix),
                         clobber=True)

        bgSL = pyregion.parse("fk5; circle({0},{1},1.0\")"
                              .format(row['x_cen'], row['y_cen']))
        bgsc = cube.subcube_from_ds9region(bgSL)
        npix = np.count_nonzero(np.isfinite(bgsc[0,:,:]))
        bgspec = (bgsc.sum(axis=(1,2)) - sc.sum(axis=(1,2))) / npix
        bgspec.meta['beam'] = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in spec.beams]),
                                              minor=np.nanmedian([bm.minor.to(u.deg).value for bm in spec.beams]),
                                              pa=np.nanmedian([bm.pa.to(u.deg).value for bm in spec.beams]),
                                             )
        bgspec.hdu.writeto("spectra/dendro{0:03d}_spw{1}_hires_background_1.0as_mean{2}.fits".format(name, spw, suffix),
                           clobber=True)

        closestpix = coord.to_pixel(cube.wcs.celestial)

        spec = cube[:, int(closestpix[1]), int(closestpix[0])]
        spec.meta['beam'] = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in spec.beams]),
                                            minor=np.nanmedian([bm.minor.to(u.deg).value for bm in spec.beams]),
                                            pa=np.nanmedian([bm.pa.to(u.deg).value for bm in spec.beams]),
                                           )
        spec.hdu.writeto("spectra/dendro{0:03d}_spw{1}_hires_pixelspec{2}.fits".format(name, spw, suffix),
                         clobber=True)
