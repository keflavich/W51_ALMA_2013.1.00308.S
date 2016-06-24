import numpy as np
from astropy import units as u
import pyregion
import radio_beam
from spectral_cube import SpectralCube

tmplt = '{reg}cax.SPW{0}_ALL.image.fits'

region_list = pyregion.open("cores_longbaseline_spectralextractionregions.reg")

for region in ('W51e2', 'W51n'):
    for spw in (2,4,6):
        try:
            cube = SpectralCube.read(tmplt.format(spw, reg=region))
        except IOError:
            print("didn't find {0}".format(tmplt.format(spw, reg=region)))
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
                try:
                    sc = cube.subcube_from_ds9region(SL)
                except ValueError as ex:
                    print(ex)
                    continue
                spec = sc.mean(axis=(1,2))
                assert not all(np.isnan(spec))

                if beam is not None:
                    spec.meta['beam'] = beam
                spec.hdu.writeto("spectra/{0}_{reg}_spw{1}_mean.fits".format(name, spw, reg=region),
                                 clobber=True)
