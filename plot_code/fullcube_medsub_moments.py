"""
Moment mapping based on full cubes
"""
from line_to_image_list import line_to_image_list
from astropy import units as u
from astropy import coordinates
from spectral_cube import SpectralCube
import pyregion
import glob
import os


try:
    import paths
    from paths import fpath,dpath
except ImportError:
    dpath = lambda x: x
    fpath = lambda x: os.path.join('moments',x)

for cubefn in glob.glob("full*fits"):
    for (linename, freq_, dv, spw) in line_to_image_list:
        freq = u.Quantity(float(freq_.strip('GHz')), unit=u.GHz)

        try:
            cube = SpectralCube.read(cubefn)
        except Exception as ex:
            print(cubefn, ex)
            raise ex
            continue

        outfn = "{0}_{1}".format(linename, cubefn)

        m0fn = paths.dpath("moments/{0}_medsub_moment0.fits".format(outfn))
        m1fn = paths.dpath("moments/{0}_medsub_moment1.fits".format(outfn))
        m2fn = paths.dpath("moments/{0}_medsub_moment2.fits".format(outfn))
        maxfn = paths.dpath("moments/{0}_medsub_max.fits".format(outfn))

        vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
                                        rest_value=freq).spectral_slab(20*u.km/u.s,
                                                                       95*u.km/u.s)

        if vcube.shape[0] < 10:
            print("Skipping {0} for line {1} because it's not in the cube".format(cubefn, linename))
            continue
        if hasattr(vcube, 'beams') and not any([bm.isfinite for bm in vcube.beams]):
            print("Skipping {0} for line {1} because of bad beams".format(cubefn, linename))
            continue

        vcube.beam_threshold = 100
        med = vcube.with_mask((((vcube.spectral_axis < 35*u.km/u.s) &
                                (vcube.spectral_axis > 20*u.km/u.s)) |
                               ((vcube.spectral_axis > 80*u.km/u.s) &
                                (vcube.spectral_axis < 95*u.km/u.s))
                              )[:,None,None]).median(axis=0)
        vcube = vcube.spectral_slab(50*u.km/u.s, 65*u.km/u.s)
        
        # I hope this isn't needed.... but it is
        if vcube.shape[0] > 200:
            raise ValueError("Warning: this is really, really big.  Danger of crash.")
        vcube.allow_huge_operations=True

        vcube_msub = vcube - med

        m0 = vcube_msub.moment0(axis=0)
        m1 = vcube_msub.moment1(axis=0)
        m2 = vcube_msub.moment2(axis=0)
        pmax = vcube_msub.max(axis=0)

        m0.hdu.writeto(m0fn, clobber=True)
        m1.hdu.writeto(m1fn, clobber=True)
        m2.hdu.writeto(m2fn, clobber=True)
        pmax.hdu.writeto(maxfn, clobber=True)
