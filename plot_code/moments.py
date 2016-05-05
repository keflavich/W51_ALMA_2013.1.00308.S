import os
import glob
from astropy import units as u
from spectral_cube import SpectralCube

try:
    from paths import fpath,dpath
    files = [
        #'12m/w51_H2CO_303_202_contsub.image.pbcor.fits',
        #'12m/w51_C18O_21_contsub.image.pbcor.fits',
        #'12m/w51_12CO_21_contsub_hires.image.pbcor.fits',
        #'12m/w51_H2CO_303_202_merge7m12m_contsub.image.pbcor.fits',
        #'12m/w51_SO_65-54_contsub.fits',
        #'12m/w51_h41alpha_contsub.image.pbcor.fits',
        #'12m/w51_H2CO_321_220_contsub.image.pbcor.fits',
        'W51_b6_12M.H2CO303_202.image.pbcor.fits',
        'W51_b6_12M.H2CO321_220.image.pbcor.fits',
        'W51_b6_12M.H2CO322_221.image.pbcor.fits',
        'W51_b6_12M.OCS18-17.image.pbcor.fits',
        'W51_b6_12M.CH3OH422-312.image.pbcor.fits',
        'W51_b6_12M.HC3N24-23.image.pbcor.fits',
        'W51_b6_7M_12M.13CS5-4.image.pbcor.fits',
        'W51_b6_7M_12M.CH3OCH3_13013-12112.image.pbcor.fits',
        'W51_b6_7M_12M.CH3OH422-312.image.pbcor.fits',
        'W51_b6_7M_12M.CH3OH423-514.image.pbcor.fits',
        'W51_b6_7M_12M.CH3OH5m42-6m43.image.pbcor.fits',
        'W51_b6_7M_12M.CH3OH808-716.image.pbcor.fits',
        'W51_b6_7M_12M.H2CO303_202.image.pbcor.fits',
        'W51_b6_7M_12M.H2CO321_220.image.pbcor.fits',
        'W51_b6_7M_12M.H2CO322_221.image.pbcor.fits',
        'W51_b6_7M_12M.HC3N24-23.image.pbcor.fits',
        'W51_b6_7M_12M.HNCO10110-919.image.pbcor.fits',
        'W51_b6_7M_12M.HNCO1028-927.image.pbcor.fits',
        'W51_b6_7M_12M.NH2CHO11210-1029.image.pbcor.fits',
        'W51_b6_7M_12M.OCS18-17.image.pbcor.fits',
        'W51_b6_7M_12M.OCS19-18.image.pbcor.fits',
        'W51_b6_7M_12M.SO65-54.image.pbcor.fits',
    ]

except ImportError:
    files = glob.glob("W51_b6*.image.pbcor.fits")
    dpath = lambda x: x
    fpath = lambda x: os.path.join('moments',x)

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

    m0fn = dpath("moments/{0}_moment0.fits".format(fname))
    m1fn = dpath("moments/{0}_moment1.fits".format(fname))
    m2fn = dpath("moments/{0}_moment2.fits".format(fname))

    if os.path.exists(m0fn):
        m0 = load_projection(m0fn)
        m1 = load_projection(m1fn)
        m2 = load_projection(m2fn)
    else:
        cube = SpectralCube.read(dpath(fn)).minimal_subcube()
        vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
        vcube = vcube.spectral_slab(35*u.km/u.s, 75*u.km/u.s)
        try:
            m0 = vcube.with_mask(vcube>0.05*u.Jy).moment0()
            m1 = vcube.with_mask(vcube>0.05*u.Jy).moment1()
            m2 = vcube.with_mask(vcube>0.05*u.Jy).moment2()
        except ValueError as ex:
            print(ex)
            vcube.beam_threshold = 1
            m0 = vcube.with_mask(vcube>0.05*u.Jy).moment0()
            m1 = vcube.with_mask(vcube>0.05*u.Jy).moment1()
            m2 = vcube.with_mask(vcube>0.05*u.Jy).moment2()


    m0.quicklook()
    m0.hdu.writeto(m0fn, clobber=True)

    m1.quicklook()
    m1.hdu.writeto(m1fn, clobber=True)

    m2.quicklook()
    m2.hdu.writeto(m2fn, clobber=True)

    try:
        m0.FITSFigure.save(fpath("moments/{0}_moment0.png".format(fname)))
        m1.FITSFigure.show_colorscale(cmap='viridis', vmin=45, vmax=68)
        m1.FITSFigure.show_contour(m0.hdu, levels=[4], colors=['k'])
        m1.FITSFigure.save(fpath("moments/{0}_moment1.png".format(fname)))
        m2.FITSFigure.show_colorscale(cmap='viridis', vmin=0, vmax=30)
        m2.FITSFigure.show_contour(m0.hdu, levels=[4], colors=['k'])
        m2.FITSFigure.save(fpath("moments/{0}_moment2.png".format(fname)))

        for ra,dec,rad,name in ((290.93253, 14.508016, 0.00311407, 'e2e8'),
                                (290.9118, 14.512366, 0.00311407, 'southeast'),
                                (290.9166, 14.518094, 0.00311407, 'north')):

            m0.FITSFigure.recenter(ra, dec, rad)
            m0.FITSFigure.save(fpath("moments/{0}_{1}zoom_moment0.png".format(fname,name)))
            m1.FITSFigure.recenter(ra, dec, rad)
            m1.FITSFigure.save(fpath("moments/{0}_{1}zoom_moment1.png".format(fname,name)))
            m2.FITSFigure.recenter(ra, dec, rad)
            m2.FITSFigure.save(fpath("moments/{0}_{1}zoom_moment2.png".format(fname,name)))

        m0.FITSFigure.close()
        m1.FITSFigure.close()
        m2.FITSFigure.close()
    except AttributeError:
        continue
