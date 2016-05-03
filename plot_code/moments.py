import os
from astropy import units as u
from spectral_cube import SpectralCube

try:
    from paths import fpath,dpath
    files = [
        '12m/w51_H2CO_303_202_contsub.image.pbcor.fits',
        '12m/w51_C18O_21_contsub.image.pbcor.fits',
        '12m/w51_12CO_21_contsub_hires.image.pbcor.fits',
        '12m/w51_H2CO_303_202_merge7m12m_contsub.image.pbcor.fits',
        '12m/w51_SO_65-54_contsub.fits',
        '12m/w51_h41alpha_contsub.image.pbcor.fits',
        '12m/w51_H2CO_321_220_contsub.image.pbcor.fits',
    ]
except ImportError:
    import glob
    files = glob.glob("W51_b6*.image.pbcor.fits")
    dpath = lambda x: x
    fpath = lambda x: os.path.join('moments',x)

for fn in files:
    cube = SpectralCube.read(dpath(fn)).minimal_subcube()
    vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
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

    fname = os.path.splitext(os.path.basename(fn))[0]

    m0.quicklook()
    m0.hdu.writeto(dpath("moments/{0}_moment0.fits".format(fname)), clobber=True)

    m1.quicklook()
    m1.hdu.writeto(dpath("moments/{0}_moment1.fits".format(fname)), clobber=True)

    m2.quicklook()
    m2.hdu.writeto(dpath("moments/{0}_moment2.fits".format(fname)), clobber=True)

    try:
        m0.FITSFigure.save(fpath("moments/{0}_moment0.png".format(fname)))
        m0.FITSFigure.close()
        m1.FITSFigure.show_colorscale(vmin=45, vmax=68)
        m1.FITSFigure.save(fpath("moments/{0}_moment1.png".format(fname)))
        m1.FITSFigure.close()
        m2.FITSFigure.show_colorscale(vmin=0, vmax=6)
        m2.FITSFigure.save(fpath("moments/{0}_moment2.png".format(fname)))
        m2.FITSFigure.close()
    except AttributeError:
        continue
