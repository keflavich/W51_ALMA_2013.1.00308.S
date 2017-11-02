import numpy as np
import pylab as pl
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import coordinates

checkregions = ['19:23:43.961 +14:30:34.693',
                '19:23:43.975 +14:30:34.655',
                '19:23:43.961 +14:30:34.249',
                '19:23:43.969 +14:30:34.516',
               ]

e2_cutoutregion = [
    '19:23:43.998 +14:30:34.043',
    '19:23:43.880 +14:30:35.114',
]
e2e_cutoutregion = [
    '19:23:43.998 +14:30:34.043',
    '19:23:43.944 +14:30:34.872',
]

def getmoments(cubefilename, vmin=40*u.km/u.s, vmax=75*u.km/u.s,
               signal_range=(50,65)*u.km/u.s,
               cutoutregion=e2e_cutoutregion,
               threshold=-6*u.mJy,
               do_check_results=True,
              ):

    cube = SpectralCube.read(cubefilename).with_spectral_unit(u.km/u.s)
    print(cube)
    print(cube.spectral_extrema)

    bl,tr = coordinates.SkyCoord(cutoutregion, frame='fk5', unit=(u.hour,
                                                                  u.deg)).transform_to(cube.wcs.wcs.radesys.lower())
    print(bl,tr)

    xxbl, yybl = map(lambda x: int(np.round(x)), cube.wcs.celestial.wcs_world2pix(bl.ra.deg, bl.dec.deg, 0))
    xxtr, yytr = map(lambda x: int(np.round(x)), cube.wcs.celestial.wcs_world2pix(tr.ra.deg, tr.dec.deg, 0))
    
    ccube = cube[:, yybl:yytr, xxbl:xxtr].mask_out_bad_beams(0.1)
    if any([x==0 for x in ccube.shape]):
        raise ValueError("Cutout out of range")
    print(ccube.spectral_extrema)

    vmask = ((ccube.spectral_axis < signal_range[0]) |
             (ccube.spectral_axis > signal_range[1]))[:,None,None]

    cont = ccube.with_mask(vmask).spectral_slab(vmin,vmax).percentile(50, axis=0)

    cscube = ccube - cont

    smask = cscube.min(axis=0) < threshold

    #m1 = cscube.with_mask(smask).moment1()
    m1 = cscube.spectral_slab(*signal_range).with_mask(smask).moment1()

    #check_results(m1, cscube, vmin=signal_range[0], vmax=signal_range[1])
    if do_check_results:
        check_results(m1, cscube, vmin=vmin, vmax=vmax)

    return m1, cscube

def check_results(m1, cscube, vmin=40*u.km/u.s, vmax=65*u.km/u.s,
                  signal_range=(50,65)*u.km/u.s):

    fig1 = pl.figure(1)
    fig1.clf()

    m1.quicklook(aplpy_kwargs={'figure':fig1})

    m1.FITSFigure.show_colorscale(vmin=vmin.to(m1.unit).value,
                                  vmax=vmax.to(m1.unit).value,
                                  cmap='spectral')

    fig2 = pl.figure(2)
    fig2.clf()

    for ii,loc in enumerate(checkregions):
        coord = coordinates.SkyCoord(loc, frame='fk5', unit=(u.hour, u.deg)).transform_to(cscube.wcs.wcs.radesys.lower())

        xx, yy = map(lambda x: int(np.round(x)), cscube.wcs.celestial.wcs_world2pix(coord.ra.deg, coord.dec.deg, 0))

        try:
            ffax = m1.FITSFigure.ax
        except AttributeError:
            ffax = m1.FITSFigure._ax1
        ffax.text(xx, yy, str(ii), color='w')
        ffax.plot(xx, yy, color='w', marker='+')

        ax = fig2.add_subplot(2,2,ii+1)
        
        cscube[:, yy, xx].quicklook()
        ax.set_title(str(ii))
        ax.plot((vmin.to(m1.unit).value, vmax.to(m1.unit).value),
                [0, 0], 'k--')
        ylim = ax.get_ylim()
        ax.vlines(signal_range.to(m1.unit).value, *ylim,
                  linestyle='--', color='k')
        ax.vlines(m1[yy,xx].value, *ylim,
                  linestyle=':', color='r')
        ax.set_xlim(vmin.to(m1.unit).value, vmax.to(m1.unit).value)
        ax.set_ylim(*ylim)

def multicube_fit(m1, cscube, signal_range=(50,65)*u.km/u.s):
    import multicube

    sc = multicube.SubCube(cscube[:,::3,::3].hdu)

    npeaks = 2
    sc.update_model('gaussian')
    sc.specfit.fitter.npeaks = npeaks
    sc.specfit.parinfo = sc.specfit.fitter.make_parinfo(npeaks=npeaks)

    minpars = [-0.03,   signal_range[0].value, 1.0]*2
    maxpars = [-0.0005, signal_range[1].value, 10.0]*2
    finesse = [5, 10, 5] + [5, 10, 5]

    sc.make_guess_grid(minpars, maxpars, finesse)
    sc.generate_model()
    sc.get_snr_map()
    sc.best_guess()

    sc.fiteach(fittype=sc.fittype,
               guesses=sc.best_guesses,
               #multicore=4,
               signal_cut=0,
               maskmap=np.isfinite(m1[::3,::3]),
               errmap=sc._rms_map,
               verbose_level=0,
               **sc.fiteach_args)

    return sc
