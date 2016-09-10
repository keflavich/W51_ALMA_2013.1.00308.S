from astropy import units as u
from astropy.coordinates import SkyCoord
import paths
import pyregion

e2e = SkyCoord('19:23:43.9661','14:30:34.497', unit=(u.hour, u.deg))
e8fil = SkyCoord('19:23:43.800','14:30:23.01', unit=(u.hour, u.deg))

e8 = SkyCoord('19:23:43.9042','14:30:28.245', unit=(u.hour, u.deg))
e8south = SkyCoord('19:23:43.783','14:30:21.15', unit=(u.hour, u.deg))

# ds9 pix -> wcs
# ww = wcs.WCS(fits.getheader('../FITS/longbaseline/SiO_m32to55kms_north.fits'))
# coordinates.SkyCoord(*ww.celestial.wcs_pix2world(264.4,258.5,1), frame='fk5', unit=(u.deg, u.deg)).to_string(style='hmsdms')
north = SkyCoord('19:23:40.0513','14:31:05.4775', unit=(u.hour, u.deg))
between_e2e_and_e8 = SkyCoord('19:23:43.911', '+14:30:30.68', unit=(u.hour, u.deg))
d2 = SkyCoord('19:23:39.819','14:31:04.83', unit=(u.hour, u.deg))
e1 = SkyCoord(290.93263,14.50745,unit=('deg','deg'), frame='fk5')

reg = pyregion.open(paths.rpath("e2eoutflow_reference_vector.reg"))[0]
e2e_reference_vector = SkyCoord(reg.coord_list[::2], reg.coord_list[1::2],
                                unit=(u.deg, u.deg))

lacy = SkyCoord('19:23:39.759', '+14:31:05.08', unit=(u.hour, u.deg))
