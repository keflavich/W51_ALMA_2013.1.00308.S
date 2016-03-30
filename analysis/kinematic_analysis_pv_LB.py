""" LB = long baseline """
import pvextractor
from astropy import units as u
from astropy import coordinates
from spectral_cube import SpectralCube
import os
import glob
import paths
import pyregion

e2ereg = pyregion.open(paths.rpath('w51e2e.reg'))[0]
w51e2e = coordinates.SkyCoord(e2ereg.coord_list[0]*u.deg, e2ereg.coord_list[1]*u.deg, frame='fk5')

files = glob.glob("/Volumes/passport/w51-alma/*.fits")


for fn in files:
    try:
        cube = SpectralCube.read(fn)
    except:
        continue
    cube.allow_huge_operations=True
    print(cube)

    out = os.path.splitext(os.path.split(fn)[1])[0]
    print(out)

    diskaxis_out = 'W51e2_PV_diskaxis_{0}.fits'.format(out)
    if not os.path.exists(diskaxis_out):
        diskycoords = "19:23:44.197,+14:30:37.34,19:23:43.960,+14:30:34.55,19:23:43.882,+14:30:32.21,19:23:43.851,+14:30:31.26".split(",")
        diskycoords = coordinates.SkyCoord(["{0} {1}".format(diskycoords[jj], diskycoords[jj+1]) for jj in (0,2,4)], unit=(u.hour, u.deg), frame='fk5')
        P = pvextractor.Path(diskycoords, 0.2*u.arcsec)
        extracted = pvextractor.extract_pv_slice(cube, P)
        extracted.writeto(diskaxis_out, clobber=True)

    outflowaxis_out = 'W51e2_PV_outflowaxis_{0}.fits'.format(out)
    if not os.path.exists(outflowaxis_out):
        outflow_coords = coordinates.SkyCoord(["19:23:44.127 +14:30:32.30", "19:23:43.822 +14:30:36.64"], unit=(u.hour, u.deg), frame='fk5')
        outflowpath = pvextractor.Path(outflow_coords, 0.2*u.arcsec)
        extracted = pvextractor.extract_pv_slice(cube, outflowpath)
        extracted.writeto(outflowaxis_out, clobber=True)
