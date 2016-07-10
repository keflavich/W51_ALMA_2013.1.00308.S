import numpy as np
from astropy import units as u
import pyregion
import radio_beam
from spectral_cube import SpectralCube

tmplt = "full_W51{2}_spw{0}{1}_lines.fits"

region_list = pyregion.open('outflow_spectralextraction_ellipses.reg') + pyregion.open("cores.reg")

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

            for reg in region_list:
                if 'text' not in reg.attr[1]:
                    continue
                name = reg.attr[1]['text']
                if name:
                    print("Extracting {0} from {1}".format(name, spw))
                    SL = pyregion.ShapeList([reg])
                    sc = cube.subcube_from_ds9region(SL)
                    spec = sc.mean(axis=(1,2))
                    assert not all(np.isnan(spec))

                    if beam is not None:
                        spec.meta['beam'] = beam
                    spec.hdu.writeto("spectra/{0}_spw{1}{2}{3}_mean.fits".format(name, spw, extra1, extra2),
                                     clobber=True)
