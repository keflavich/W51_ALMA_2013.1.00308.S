import numpy as np
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
from astropy import wcs
import pyregion
import paths

fh = fits.open(paths.dpath('longbaseline/W51e2cax.cont.image.pbcor.fits'))
mywcs = wcs.WCS(fh[0].header).sub([wcs.WCSSUB_CELESTIAL])
pixscale = (mywcs.pixel_scale_matrix**2).sum()**0.5
reg = pyregion.open(paths.rpath('cores_longbaseline_spectralextractionregions_pix.reg'))

with open(paths.rpath('cores_longbaseline_spectralextractionregions.reg'), 'w') as rfh:
    rfh.write("global color=red\n")
    rfh.write("fk5\n")

    for rr in reg:
        if rr.name == 'circle':
            x,y = mywcs.wcs_pix2world(rr.coord_list[0], rr.coord_list[1], 1)
            r = pixscale * rr.coord_list[2]
            rfh.write("circle({0}, {1}, {2}) # text={{{3}}}\n".format(x, y, r, rr.attr[1]['text']))
        elif rr.name == 'ellipse':
            x,y = mywcs.wcs_pix2world(rr.coord_list[0], rr.coord_list[1], 1)
            maj,min = pixscale * np.array(rr.coord_list[2:4])
            pa = rr.coord_list[4]
            rfh.write("ellipse({0}, {1}, {2}, {3}, {4}) # text={{{5}}}\n".format(x, y, maj, min, pa, rr.attr[1]['text']))
