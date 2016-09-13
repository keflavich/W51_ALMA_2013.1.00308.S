import numpy as np
from astropy import coordinates
from astropy import units as u
import radio_beam
from spectral_cube import SpectralCube
from astropy.table import Table

tmplt = "full_W51{2}_spw{0}{1}_lines.fits"

tbl = Table.read('continuum_photometry.ipac', format='ascii.ipac')

for extra1 in ("","_7m12m"):
    for extra2 in ("","_hires"):
        for spw in (0,1,2,3):
            try:
                cube = SpectralCube.read(tmplt.format(spw, extra2, extra1))
            except IOError:
                print("didn't find {0}".format(tmplt.format(spw, extra2, extra1)))
                continue
            print(cube)
            try:
                beam = radio_beam.Beam.from_fits_header(cube.header)
            except TypeError:
                if hasattr(cube, 'beams'):
                    beam = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in cube.beams]),
                                           minor=np.nanmedian([bm.minor.to(u.deg).value for bm in cube.beams]),
                                           pa=np.nanmedian([bm.pa.to(u.deg).value for bm in cube.beams]),
                                          )
                else:
                    beam = None

            for row in tbl:
                name = row['name']

                print("Extracting {0} from {1}".format(name, spw))

                coord = coordinates.SkyCoord(row['PeakRA'], row['PeakDec'],
                                             frame='fk5', unit=(u.deg, u.deg))
                xpix, ypix = cube.wcs.celestial.wcs_world2pix(coord.ra.deg,
                                                              coord.dec.deg,
                                                              0)
                spec = cube[:,int(ypix),int(xpix)]
                assert not all(np.isnan(spec))

                if beam is not None:
                    spec.meta['beam'] = beam
                spec.hdu.writeto("spectra/{0}_spw{1}{2}{3}_peak.fits".format(name, spw, extra1, extra2),
                                 clobber=True)
