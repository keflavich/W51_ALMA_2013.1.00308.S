import numpy as np
import paths
from astropy.io import fits
from astropy import units as u
from astropy import coordinates
from astropy import wcs
import image_tools
import pylab as pl


for sourcename, region, center in (('e2','e2e8',("19:23:43.963","14:30:34.53")),
                                   ('e8','e2e8',("19:23:43.898","14:30:28.29")),
                                   ('north','north',("19:23:40.050","14:31:05.48")),
                                   ('ALMAmm14','ALMAmm14',("19:23:38.571","14:30:41.77")),
                                  ):

    hdu = fits.open(paths.dpath('12m/moments/CH3OH_{0}_cutout_temperaturemap.fits'.format(sourcename)))
    tmap = hdu[0].data
    header = hdu[0].header
    mywcs = wcs.WCS(header)
    crd = coordinates.SkyCoord(center[0],center[1],frame='fk5', unit=(u.hour,u.deg))
    pixcen = mywcs.wcs_world2pix(crd.ra.deg, crd.dec.deg, 0)

    weights = ((tmap > 0) & (tmap < 1000)).astype('float')

    nr, bins, rprof = image_tools.radialprofile.azimuthalAverage(tmap,
                                                                 binsize=1.0,
                                                                 weights=weights,
                                                                 center=pixcen,
                                                                 return_nr=True)
    mywcs = wcs.WCS(header)
    pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5
    pl.figure(4).clf()
    pl.plot(bins*pixscale*3600, rprof)
    pl.ylim(0,1000)
    pl.xlim(0,2.5)
    pl.xlabel("Radius (arcsec)")
    pl.ylabel("Average Temperature (K)")
    pl.savefig(paths.fpath("chemistry/ch3oh_temperature_radial_profile_{0}.png".format(sourcename)))
