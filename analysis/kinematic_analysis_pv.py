import pvextractor
from astropy import units as u
from astropy import coordinates
from spectral_cube import SpectralCube

fntemplate = 'full_W51e2cutout_spw{0}_lines.fits'
medsubtemplate = 'full_W51e2cutout_spw{0}_lines_medsub.fits'

for ii in range(3):
    cube = SpectralCube.read(fntemplate.format(ii))
    cube.allow_huge_operations=True
    med = cube.median(axis=0)
    medsub = cube - med
    medsub.write(medsubtemplate.format(ii), overwrite=True)

    diskycoords = "19:23:44.197,+14:30:37.34,19:23:43.960,+14:30:34.55,19:23:43.882,+14:30:32.21,19:23:43.851,+14:30:31.26".split(",")
    diskycoords = coordinates.SkyCoord(["{0} {1}".format(diskycoords[ii], diskycoords[ii+1]) for ii in (0,2,4)], unit=(u.hour, u.deg), frame='fk5')
    P = pvextractor.Path(diskycoords, 0.2*u.arcsec)
    extracted = pvextractor.extract_pv_slice(medsub, P)
    extracted.writeto('W51e2_PV_diskaxis_spw{0}.fits'.format(ii), clobber=True)

    outflow_coords = coordinates.SkyCoord(["19:23:44.127 +14:30:32.30", "19:23:43.822 +14:30:36.64"], unit=(u.hour, u.deg), frame='fk5')
    outflowpath = pvextractor.Path(outflow_coords, 0.2*u.arcsec)
    extracted = pvextractor.extract_pv_slice(medsub, outflowpath)
    extracted.writeto('W51e2_PV_outflowaxis_spw{0}.fits'.format(ii), clobber=True)
