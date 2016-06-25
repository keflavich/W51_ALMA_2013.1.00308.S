import os
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from astropy import units as u
import paths
from astropy.io import fits
from astropy import wcs
from astropy import log


def get_mom0(fn, vrange=[51,60]*u.km/u.s, exclude_vrange=[40,75]*u.km/u.s,
             percentile=30, iterate=True):
    savename = paths.dpath(os.path.join(os.path.split(fn)[0], 'moments',
                                        os.path.split(fn)[1].replace(".fits","_p30sub.fits")))
    if os.path.exists(savename):
        print("Loading {0} from disk".format(savename))
        fh = fits.open(savename)

        m0 = Projection(value=fh[0].data, header=fh[0].header,
                        wcs=wcs.WCS(fh[0].header),
                        unit=u.Unit(fh[0].header['BUNIT']),)

    else:
        print("Computing {0} from analysis of {1}".format(savename, fn))
        cube = SpectralCube.read(fn)
        #bm = cube.beams[0]
        #jtok = bm.jtok(cube.wcs.wcs.restfrq*u.Hz)

        slab = cube.spectral_slab(*vrange)
        cube.beam_threshold = 1
        #contguess = cube.spectral_slab(0*u.km/u.s, 40*u.km/u.s).percentile(50, axis=0)
        #contguess = cube.spectral_slab(70*u.km/u.s, 100*u.km/u.s).percentile(50, axis=0)
        mask = (cube.spectral_axis<exclude_vrange[0]) | (cube.spectral_axis > exclude_vrange[1])
        masked_cube = cube.with_mask(mask[:,None,None])
        try:
            log.info("Continuum measurement.")
            contguess = masked_cube.percentile(percentile, axis=0,
                                               iterate_rays=iterate,
                                               progressbar=iterate)
        except ValueError as ex:
            print("skipping {0}".format(fn))
            print(ex)
            raise
        log.info("Subtract the continuum")
        slabsub = (slab-contguess)
        slabsub.beam_threshold = 0.25
        log.info("Compute moment 0")
        m0 = slabsub.moment0()
        log.info("Save results to {0}".format(savename))
        m0.hdu.writeto(savename)

    return m0
