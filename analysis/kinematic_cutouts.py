from spectral_cube import SpectralCube
import pyregion
import glob
import os
from astropy import units as u

# this is meant to be run on orion

files = glob.glob("W51_b6*.image.pbcor.fits")

for region_fn,name in (('w51e2e8box_kinematics.reg', 'e2e8'),
                       ('w51northbox_kinematics.reg', 'north'),
                       ('ALMAmm14box_kinematics.reg', 'ALMAmm14'),
                      ):
    regions = pyregion.open(region_fn)

    for fn in files:
        outname = fn.replace(".fits","_{0}cutout.fits".format(name))
        if not os.path.exists(outname):
            print(outname)
            cube = SpectralCube.read(fn).minimal_subcube()
            vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')

            cutout = vcube.subcube_from_ds9region(regions)

            cutout.write(outname,
                         overwrite=True)
