"""
Make bite-sized cutouts of the long-baseline data
"""
import re
from astropy import units as u
from astropy import coordinates
from spectral_cube import SpectralCube
import pyregion
import paths

view = slice(None), slice(2400,3100), slice(2300,2850)

regions = pyregion.open(paths.rpath("e2e8northcutouts.reg"))

corners = {reg.attr[1]['text']: {'lowerleft': coordinates.SkyCoord([reg.coord_list[:2]], frame='fk5', unit=(u.deg, u.deg)),
                                 'upperright': coordinates.SkyCoord([reg.coord_list[2:4]], frame='fk5', unit=(u.deg, u.deg)),}
           for reg in regions
          }

#corners = {'e2':dict(lowerleft=coordinates.SkyCoord("19:23:44.002 +14:30:33.82", frame='fk5', unit=(u.hour, u.deg)),
#                     upperright=coordinates.SkyCoord("19:23:43.813 +14:30:37.32", frame='fk5', unit=(u.hour, u.deg))),
#           #'north':dict(lowerleft=coordinates.SkyCoord("19:23: +14:31", frame='fk5', unit=(u.hour, u.deg)),
#           #             upperright=coordinates.SkyCoord("19:23: +14:31", frame='fk5', unit=(u.hour, u.deg))),
#          }
repl = re.compile("W51[en]2?cax")

for source,cubefn in [#('e2', "W51e2cax.CH3CN_K3_nat.image.fits"),
                      #('e2', "W51e2cax.CH3CN_K3_nat_all.image.fits"),
                      # done ('e2', "W51e2cax.CH3CN_K8.image.pbcor.fits"),
                      # done ('e2', "W51e2cax.H30alpha.image.pbcor.fits"),
                      # done ('e2', "W51e2cax.SPW2_ALL.image.fits"),
                      # done ('e2', "W51e2cax.SPW4_ALL.image.fits"),
                      # done ('e2', "W51e2cax.SPW6_ALL.image.fits"),
                      # done ('e2', "W51e2cax.SPW1_ALL.image.fits"),
                      # done ('e2', "W51e2cax.SPW3_ALL.image.fits"),
                      # done ('e2', "W51e2cax.SPW5_ALL.image.fits"),
                      # done ('e2', "W51e2cax.SPW7_ALL.image.fits"),
                      # done ('e2', "W51e2cax.SPW8_ALL.image.fits"),
                      ('e2', "W51e2cax.SPW9_ALL.image.fits"),
                      #('e8', "W51e2cax.SPW1_ALL.image.fits"),
                      #('e8', "W51e2cax.SPW3_ALL.image.fits"),
                      #('e8', "W51e2cax.SPW5_ALL.image.fits"),
                      #('e8', "W51e2cax.SPW7_ALL.image.fits"),
                      #('e8', "W51e2cax.SPW8_ALL.image.fits"),
                      #('e8', "W51e2cax.SPW9_ALL.image.fits"),
                      #('e8', "W51e2cax.CH3CN_K3_nat.image.fits"),
                      #('e8', "W51e2cax.CH3CN_K3_nat_all.image.fits"),
                      # done ('e8', "W51e2cax.CH3CN_K8.image.pbcor.fits"),
                      # done ('e8', "W51e2cax.H30alpha.image.pbcor.fits"),
                      # done ('e8', "W51e2cax.SPW2_ALL.image.fits"),
                      # done ('e8', "W51e2cax.SPW4_ALL.image.fits"),
                      # done ('e8', "W51e2cax.SPW6_ALL.image.fits"),
                      # done ('north', "W51ncax.H30alpha.image.pbcor.fits"),
                      # done ('north', "W51ncax.CH3CN_K8.image.pbcor.fits"),
                      # done ('north', "W51ncax.SPW2_ALL.image.fits"),
                      # done ('north', "W51ncax.SPW4_ALL.image.fits"),
                      # done ('north', "W51ncax.SPW6_ALL.image.fits"),
                      ]:
    suffix = ".image.fits" if ".image.fits" in cubefn else ".image.pbcor.fits"
    print("Cube {0} cutout {1} beginning".format(cubefn, source))

    try:
        cube = SpectralCube.read(cubefn)
    except Exception as ex:
        print(cubefn, ex)
        continue
        raise ex

    outfn = repl.sub("W51{0}cax".format(source), cubefn)

    lowerleft, upperright = corners[source]['lowerleft'],corners[source]['upperright'],
    bl_x, bl_y = cube.wcs.celestial.wcs_world2pix(lowerleft.ra, lowerleft.dec, 0)
    tr_x, tr_y = cube.wcs.celestial.wcs_world2pix(upperright.ra, upperright.dec, 0)
    assert tr_y > bl_y+2
    assert tr_x > bl_x+2

    view = slice(None), slice(bl_y, tr_y), slice(bl_x, tr_x)
    print(cube)
    print("View: ",view)

    cutout = cube[view]
    print(cutout)
    cutout.write(outfn.replace(suffix,"_cutout.fits"), overwrite=True)
    print("Median calculation")
    med = cutout.median(axis=0)
    print("Median subtraction")
    cutout.allow_huge_operations=True
    cutoutms = cutout-med
    print("Velocity conversion")
    vcutoutms = cutoutms.with_spectral_unit(u.km/u.s,
                                            velocity_convention='radio')
    print("Writing velocity cutout",vcutoutms)
    vcutoutms.write(outfn.replace(suffix,"_medsub_cutout.fits"), overwrite=True)
