import os
import glob

import paths
from spectral_cube import SpectralCube
from astropy import units as u
import radio_beam

sourcename = 'e2e8'

out_beam = radio_beam.Beam(0.35*u.arcsec, 0.35*u.arcsec)

for fn in glob.glob(paths.dpath('12m/cutouts/*.CH3OH*{0}*fits'.format(sourcename))):

    outpath = fn.replace(".fits","_0.35arcsec.fits").replace("cutouts","cutouts/commonres")
    if not os.path.exists(outpath):

        vrcube = SpectralCube.read(fn)
        cube = vrcube.convolve_to(out_beam)
        cube.write(outpath)

cubes = [SpectralCube.read(fn) for fn in glob.glob(paths.dpath('12m/cutouts/commonres/W51_b6_12M.CH3OH*.image.pbcor_{0}cutout_0.35arcsec.fits'
                                                               .format(sourcename)))]
