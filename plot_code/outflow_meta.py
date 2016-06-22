from astropy import units as u
from astropy.coordinates import SkyCoord
import paths
import pyregion

e2e = SkyCoord('19:23:43.967','14:30:34.54', unit=(u.hour, u.deg))
e8fil = SkyCoord('19:23:43.800','14:30:23.01', unit=(u.hour, u.deg))
north = SkyCoord('19:23:40.065','14:31:05.33', unit=(u.hour, u.deg))
between_e2e_and_e8 = SkyCoord('19:23:43.911', '+14:30:30.68', unit=(u.hour, u.deg))

reg = pyregion.open(paths.rpath("e2eoutflow_reference_vector.reg"))[0]
e2e_reference_vector = SkyCoord(reg.coord_list[::2], reg.coord_list[1::2],
                                unit=(u.deg, u.deg))
