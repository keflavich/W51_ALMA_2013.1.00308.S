import os
import glob
import collections
from astropy import units as u
from spectral_cube import SpectralCube
#from astropy.stats import mad_std
from mad import mad_std_nan
import numpy as np

import matplotlib
matplotlib.use('Agg')

try:
    from paths import fpath,dpath
    files = [
        'W51_b6_12M.H2CO303_202.image.pbcor.fits',
        'W51_b6_12M.H2CO321_220.image.pbcor.fits',
        'W51_b6_12M.H2CO322_221.image.pbcor.fits',
        'W51_b6_12M.OCS18-17.image.pbcor.fits',
        'W51_b6_12M.CH3OH422-312.image.pbcor.fits',
        'W51_b6_12M.HC3N24-23.image.pbcor.fits',
        'W51_b6_12M.HNCO10110-919.image.pbcor.fits',
        'W51_b6_12M.OCS19-18.image.pbcor.fits',
        'W51_b6_12M.SO65-54.image.pbcor.fits',
        'W51_b6_12M.HNCO1028-927.image.pbcor.fits',
        'W51_b6_12M.CH3OH423-514.image.pbcor.fits',
        'W51_b6_12M.CH3OH5m42-6m43.image.pbcor.fits',
        'W51_b6_12M.CH3OH808-716.image.pbcor.fits',
        'W51_b6_12M.CH3OH1029-936.image.pbcor.fits',
        'W51_b6_12M.CH3OH18315-17414.image.pbcor.fits',
        'W51_b6_12M.CH3OH23519-22617.image.pbcor.fits',
        'W51_b6_12M.CH3OH25322-24420.image.pbcor.fits',
        'W51_b6_12M.13CS5-4.image.pbcor.fits',
        'W51_b6_12M.NH2CHO11210-1029.image.pbcor.fits',
        'W51_b6_12M.NH2CHO1156-1055.image.pbcor.fits',
        'W51_b6_12M.HC3Nv7=124-23.image.pbcor.fits',
        'W51_b6_12M.H30alpha.image.pbcor.fits',
        'W51_b6_12M.C18O2-1.image.pbcor.fits',
        'W51_b6_12M.H2CCO11-10.image.pbcor.fits',
        'W51_b6_12M.HCOOH431-524.image.pbcor.fits',
        'W51_b6_12M.CH3OCHO17314-16313E.image.pbcor.fits',
        'W51_b6_12M.CH3CH2CN24321-23320.image.pbcor.fits',
        'W51_b6_12M.HC3Nv7=1_24-23.image.pbcor.fits',
        'W51_b6_12M.Acetone21120-20219AE.image.pbcor.fits',
        'W51_b6_12M.Acetone21120-20119EE.image.pbcor.fits',
        'W51_b6_12M.CH3CH2CN24222-23221.image.pbcor.fits',
        'W51_b6_12M.H213CO312-211.image.pbcor.fits',
        'W51_b6_12M.H2CN303-202_F3_2.image.pbcor.fits',
        'W51_b6_12M.H2CN322-221_F5_2.image.pbcor.fits',
        'W51_b6_12M.CH3OCHO17413-16412A.image.pbcor.fits',
        'W51_b6_12M.CH3CH2OH550-541.image.pbcor.fits',
        'W51_b6_12M.CH3OCH313013-12112AA.image.pbcor.fits',
        'W51_b6_12M.HNCO10010-909.image.pbcor.fits',
        'W51_b6_12M.CH3OCHO17314-16313A.image.pbcor.fits',
        'W51_b6_12M.O13CS18-17.image.pbcor.fits',
        'W51_b6_12M.CH3OCH323321-23222AA.image.pbcor.fits',
        'W51_b6_12M.CH3OCH323321-23222EE.image.pbcor.fits',
        'W51_b6_12M.N2D+_3-2.image.pbcor.fits',
        'W51_b6_12M.12CO2-1.image.pbcor.fits',
        'W51_b6_12M.HNCO1055-954.image.pbcor.fits',
        'W51_b6_12M.HNCO1046-945.image.pbcor.fits',
        'W51_b6_12M.HNCO1038-937.image.pbcor.fits',
        'W51_b6_12M.13CH3OH515-414.image.pbcor.fits',
        'W51_b6_12M.PN5-4.image.pbcor.fits',
        'W51_b6_12M.HNCO28128-29029.image.pbcor.fits',
        'W51_b6_12M.SO2_22715-23618.image.pbcor.fits',
        'W51_b6_12M.SO2_16610-17513.image.pbcor.fits',
        'W51_b6_12M.SO2v2=1_20218-19317.image.pbcor.fits',
        'W51_b6_12M.SO2v2=1_22220-22121.image.pbcor.fits',
        'W51_b6_12M.SO2v2=1_16313-16214.image.pbcor.fits',
        'W51_b6_12M.SO2v2=1_642-735.image.pbcor.fits',
        'W51_b6_12M.CH3NCO_25124-24123.image.pbcor.fits',
        'W51_b6_12M.CH3NCO_27226-26225.image.pbcor.fits',
        'W51_b6_12M.CH3SH_152-151.image.pbcor.fits',
        'W51_b6_12M.CH3SH_162-161.image.pbcor.fits',
        'W51_b6_12M.CH3SH_73-82.image.pbcor.fits',
        'W51_b6_12M.CH3SH_232-231.image.pbcor.fits',
    ]

except ImportError:
    print("Failed to import dpath")
    files = glob.glob("W51_b6*.image.pbcor.fits")
    dpath = lambda x: x
    fpath = lambda x: os.path.join('moments',x)

cont_percentiles = collections.defaultdict(lambda: 50)
# these lines almost intersect with SO, causing problems near outflows
cont_percentiles['CH3OH23519-22617'] = 10
cont_percentiles['CH3OH25322-24420'] = 10
cont_percentiles['CH3OH18315-17414'] = 10
cont_percentiles['HNCO28128-29029'] = 10

def load_projection(filename):
    from astropy import wcs
    from astropy.io import fits
    from spectral_cube.lower_dimensional_structures import Projection

    fh = fits.open(filename)

    return Projection(value=fh[0].data,
                      header=fh[0].header,
                      wcs=wcs.WCS(fh[0].header),
                      unit=u.Unit(fh[0].header['BUNIT']),
                     )


for fn in files:
    fname = os.path.splitext(os.path.basename(fn))[0]
    print(fname)

    m0fn = dpath("moments/{0}_medsub_moment0.fits".format(fname))
    m1fn = dpath("moments/{0}_medsub_moment1.fits".format(fname))
    m2fn = dpath("moments/{0}_medsub_moment2.fits".format(fname))
    madstdfn = dpath("moments/{0}_medsub_madstd.fits".format(fname))
    maxfn = dpath("moments/{0}_medsub_max.fits".format(fname))
    argmaxfn = dpath("moments/{0}_medsub_argmax.fits".format(fname))

    if os.path.exists(m0fn) and os.path.exists(madstdfn):
        m0 = load_projection(m0fn)
        m1 = load_projection(m1fn)
        m2 = load_projection(m2fn)
        pmax = load_projection(maxfn)
        madstd = load_projection(madstdfn)
#        argmaxfn = dpath("moments/{0}_medsub_argmax.fits".format(fname))
    else:
        cube = SpectralCube.read(dpath(fn)).minimal_subcube()
        vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
        vcube.beam_threshold = 100

        pct = 50
        for key in cont_percentiles:
            if key in fname:
                pct = cont_percentiles[key]

        contmasked_cube = vcube.with_mask(((vcube.spectral_axis < 35*u.km/u.s)
                                           | (vcube.spectral_axis >
                                              80*u.km/u.s))[:,None,None])

        med = contmasked_cube.percentile(pct, axis=0)
        vcube = vcube.spectral_slab(50*u.km/u.s, 65*u.km/u.s)
        # I hope this isn't needed....
        vcube.allow_huge_operations=True
        vcube_msub = vcube - med

        m0 = vcube_msub.moment0(axis=0)
        m1 = vcube_msub.moment1(axis=0)
        m2 = vcube_msub.moment2(axis=0)
        # required if mad_std can't handle NaNs
        #madstd = contmasked_cube.apply_function(mad_std, axis=0,
        #                                        projection=True,
        #                                        progressbar=True,
        #                                        unit=cube.unit,)
        madstd = contmasked_cube.apply_numpy_function(mad_std_nan,
                                                      axis=0,
                                                      projection=True,
                                                      unit=cube.unit,
                                                      fill=np.nan,
                                                      how='cube',
                                                      reduce=True,
                                                      progressbar=True)
        pmax = vcube_msub.max(axis=0)

#        argmax = vcube.spectral_axis[vcube_msub.argmax(axis=0)]


    m0.quicklook()
    m0.hdu.writeto(m0fn, clobber=True)

    m1.quicklook()
    m1.hdu.writeto(m1fn, clobber=True)

    m2.quicklook()
    m2.hdu.writeto(m2fn, clobber=True)

    pmax.quicklook()
    pmax.hdu.writeto(maxfn, clobber=True)

    madstd.quicklook()
    madstd.hdu.writeto(madstdfn, clobber=True)

#    argmax.quicklook()
#    argmax.hdu.writeto(argmaxfn, clobber=True)

    try:
        m0.FITSFigure.save(fpath("moments/{0}_moment0.png".format(fname)))
        m1.FITSFigure.show_colorscale(cmap='viridis', vmin=45, vmax=68)
        m1.FITSFigure.show_contour(m0.hdu, levels=[4], colors=['k'])
        m1.FITSFigure.save(fpath("moments/{0}_moment1.png".format(fname)))
        m2.FITSFigure.show_colorscale(cmap='viridis', vmin=0, vmax=30)
        m2.FITSFigure.show_contour(m0.hdu, levels=[4], colors=['k'])
        m2.FITSFigure.save(fpath("moments/{0}_moment2.png".format(fname)))
        pmax.FITSFigure.show_colorscale(cmap='viridis')
        pmax.FITSFigure.show_contour(m0.hdu, levels=[4], colors=['k'])
        pmax.FITSFigure.save(fpath("moments/{0}_max.png".format(fname)))
#        argmax.FITSFigure.show_colorscale(cmap='RdYlBu_r', vmin=45, vmax=68)
#        argmax.FITSFigure.save(fpath("moments/{0}_argmax.png".format(fname)))

        for ra,dec,rad,name in ((290.93253, 14.508016, 0.00311407, 'e2e8'),
                                (290.9118, 14.512366, 0.00311407, 'southeast'),
                                (290.9166, 14.518094, 0.00311407, 'north')):

            m0.FITSFigure.recenter(ra, dec, rad)
            m0.FITSFigure.save(fpath("moments/{0}_{1}zoom_moment0.png".format(fname,name)))
            m1.FITSFigure.recenter(ra, dec, rad)
            m1.FITSFigure.save(fpath("moments/{0}_{1}zoom_moment1.png".format(fname,name)))
            m2.FITSFigure.recenter(ra, dec, rad)
            m2.FITSFigure.save(fpath("moments/{0}_{1}zoom_moment2.png".format(fname,name)))
            pmax.FITSFigure.recenter(ra, dec, rad)
            pmax.FITSFigure.save(fpath("moments/{0}_{1}zoom_max.png".format(fname,name)))
#            argmax.FITSFigure.recenter(ra, dec, rad)
#            argmax.FITSFigure.save(fpath("moments/{0}_{1}zoom_argmax.png".format(fname,name)))

        m0.FITSFigure.close()
        m1.FITSFigure.close()
        m2.FITSFigure.close()
        pmax.FITSFigure.close()
#        argmax.FITSFigure.close()
    except AttributeError:
        continue
