"""
Make bite-sized cutouts of full cube data around the main sources
"""
from astropy import units as u
from astropy import coordinates
from spectral_cube import SpectralCube
import pyregion
import glob

try:
    import paths
    regions = pyregion.open(paths.rpath("e2e8northcutouts.reg"))
except ImportError:
    regions = pyregion.open("e2e8northcutouts.reg")

corners = {reg.attr[1]['text']:
           {'lowerleft': coordinates.SkyCoord([reg.coord_list[:2]],
                                              frame='fk5', unit=(u.deg, u.deg)),
            'upperright': coordinates.SkyCoord([reg.coord_list[2:4]],
                                               frame='fk5', unit=(u.deg, u.deg)),}
           for reg in regions
          }

for source in ('e2','e8','north'):
    for cubefn in glob.glob("full*fits"):

        try:
            cube = SpectralCube.read(cubefn)
        except Exception as ex:
            print(cubefn, ex)
            raise ex
            continue

        outfn = "{0}cutout_{1}".format(source, cubefn)

        lowerleft, upperright = corners[source]['lowerleft'],corners[source]['upperright']
        bl_x, bl_y = cube.wcs.celestial.wcs_world2pix(lowerleft.ra, lowerleft.dec, 0)
        tr_x, tr_y = cube.wcs.celestial.wcs_world2pix(upperright.ra, upperright.dec, 0)
        assert tr_y > bl_y+2
        assert tr_x > bl_x+2

        view = slice(None), slice(bl_y, tr_y), slice(bl_x, tr_x)
        print(cube)
        print("View: ",view)

        cutout = cube[view]
        print(cutout)
        cutout.write(outfn, overwrite=True)
