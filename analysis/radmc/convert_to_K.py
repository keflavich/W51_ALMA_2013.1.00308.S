import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import wcs

def convert_to_K(radmc_fits_img, distance=5400*u.pc):
    fh = fits.open(radmc_fits_img)
    mywcs = wcs.WCS(fh[0].header)
    pix_area = np.abs(mywcs.celestial.pixel_scale_matrix.diagonal().prod()) * u.deg**2
    conv = u.Jy.to(u.K, equivalencies=u.brightness_temperature(pix_area, mywcs.wcs.crval[2]*u.Hz))
    fh[0].data *= conv
    fh[0].header['BUNIT'] = 'K'
    return fh
