import numpy as np
from astropy import coordinates
from astropy import units as u
import copy
import pyregion
import radio_beam
from spectral_cube import SpectralCube
import os

tmplt = "full_W51{2}_spw{0}{1}_lines.fits"

def rpath(x, basepath='/Users/adam/work/w51/alma/regions/'):
    return os.path.join(basepath,x)

#region_list = (pyregion.open(rpath('outflow_spectralextraction_ellipses.reg')) +
#               pyregion.open(rpath("cores.reg")) +
#               pyregion.open(rpath('cm_continuum_apertures.reg')))
region_list = (pyregion.open(rpath('vla_pointsource_centroids.reg')))

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
                if name and reg.name in ('circle',):
                    print("Extracting {0} from {1}".format(name, spw))
                    SL = pyregion.ShapeList([reg])
                    sc = cube.subcube_from_ds9region(SL)
                    mask = sc.mask.include().max(axis=0)
                    spec = sc.mean(axis=(1,2))
                    assert not all(np.isnan(spec))

                    # make a 'background region' that has the same area
                    bgreg = copy.copy(reg)
                    bgreg.coord_list[2] *= 2**0.5
                    SLbg = pyregion.ShapeList([bgreg])
                    scbg = cube.subcube_from_ds9region(SLbg)
                    bgspec = (scbg.sum(axis=(1,2)) - sc.sum(axis=(1,2))) / mask.sum()

                    if beam is not None:
                        spec.meta['beam'] = beam
                        bgspec.meta['beam'] = beam
                    spec.hdu.writeto("spectra/{0}_spw{1}{2}{3}_mean.fits".format(name, spw, extra1, extra2),
                                     overwrite=True)
                    bgspec.hdu.writeto("spectra/{0}_spw{1}{2}{3}_background_mean.fits".format(name, spw, extra1, extra2),
                                       overwrite=True)
                elif name and reg.name in ('point',):
                    print("Extracting {0} from {1}".format(name, spw))
                    coord = coordinates.SkyCoord(reg.coord_list[0], reg.coord_list[1],
                                                 frame='fk5', unit=(u.deg, u.deg))
                    xpix, ypix = cube.wcs.celestial.wcs_world2pix(coord.ra.deg,
                                                                  coord.dec.deg,
                                                                  0)
                    try:
                        spec = cube[:,int(np.round(ypix)),int(np.round(xpix))]
                    except IndexError:
                        print("Skipping {0} because it is outside the cube."
                              "  xpix={1} ypix={2} ra={3} dec={4}"
                              .format(name, xpix, ypix, coord.ra, coord.dec))
                        continue
                    assert not all(np.isnan(spec))

                    if beam is not None:
                        spec.meta['beam'] = beam
                    spec.hdu.writeto("spectra/{0}_spw{1}{2}{3}_singlepixelspectrum.fits"
                                     .format(name, spw, extra1, extra2),
                                     overwrite=True)
