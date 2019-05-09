from astropy.nddata import Cutout2D
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy.wcs import utils as wcsutils
from astropy import stats
import reproject
import image_registration

cont3mm = fits.open('/Users/adam/work/w51/alma/FITS/longbaseline/w51e2_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits.gz')
cont7mm = fits.open('/Users/adam/work/w51/vla_q/FITS/W51e2w_QbandAarray_cont_spws_continuum_cal_clean_2terms_robust0_wproj_selfcal9.image.tt0.pbcor.fits')

cutout_center = coordinates.SkyCoord('19:23:43.941 +14:30:34.559',
                                     unit=(u.hour, u.deg), frame='icrs')
size = 2*u.arcsec

cutout_3mm = Cutout2D(cont3mm[0].data.squeeze(), cutout_center, size, wcs=wcs.WCS(cont3mm[0].header).celestial)
cutout_7mm = Cutout2D(cont7mm[0].data.squeeze(), cutout_center, size, wcs=wcs.WCS(cont7mm[0].header).celestial)
proj_7mmto3mm,_ = reproject.reproject_interp((cutout_7mm.data, cutout_7mm.wcs),
                                             cutout_3mm.wcs,
                                             shape_out=cutout_3mm.shape)

pixscale = wcsutils.proj_plane_pixel_area(cutout_3mm.wcs)**0.5*u.deg

errest = stats.mad_std(cutout_3mm.data)

chi2shift = image_registration.chi2_shift(proj_7mmto3mm, cutout_3mm.data, err=errest, upsample_factor=1000)
print(chi2shift)
print(chi2shift[:2] * cutout_3mm.wcs.wcs.cdelt * 3600)
print(chi2shift[:2] * cutout_3mm.wcs.wcs.cdelt)
"""
[-4.295500000000004, 3.9625000000000057, 0.0034999999999999892, 0.003500000000000003]
[0.0214775 0.0198125]
[5.96597222e-06 5.50347222e-06]
"""
ichi2shift = image_registration.chi2_shift_iterzoom(proj_7mmto3mm, cutout_3mm.data, err=errest, upsample_factor=1000)
print(ichi2shift)

from image_registration.fft_tools import shift
xoff, yoff = chi2shift[:2]
corrected_7mm = shift.shiftnd(proj_7mmto3mm, (yoff, xoff))

import pylab as pl
pl.figure(1).clf()
pl.subplot(2,2,1).imshow(proj_7mmto3mm, origin='lower')
pl.subplot(2,2,2).imshow(cutout_3mm.data, origin='lower')
pl.subplot(2,2,3).imshow(proj_7mmto3mm * 3 - cutout_3mm.data, origin='lower')
pl.subplot(2,2,4).imshow(corrected_7mm * 3 - cutout_3mm.data, origin='lower')
