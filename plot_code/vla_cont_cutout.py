import paths
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs

fnku = paths.dpath('evla/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits')
fitsKu = fits.open(fnku)
cutout_Ku = Cutout2D(fitsKu[0].data.squeeze(),
                     coordinates.SkyCoord('19:23:41.495','+14:30:40.48',unit=(u.hour,u.deg)),
                     2.22*u.arcmin, wcs=wcs.WCS(fitsKu[0].header).celestial)
fitsKu_cutout = fits.PrimaryHDU(data=cutout_Ku.data, header=cutout_Ku.wcs.to_header())
fitsKu_fn = paths.dpath("rgb/Kuband_e2e_cutout.fits")
fitsKu_cutout.writeto(fitsKu_fn, clobber=True)
