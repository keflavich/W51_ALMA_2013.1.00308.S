import paths
import pylab as pl
from astropy.io import fits
from astropy.wcs import WCS
from astropy import wcs
from astropy.nddata import Cutout2D
from spectral_cube import SpectralCube
import numpy as np
from astropy import units as u
from outflow_meta import e2e, between_e2e_and_e8, e8fil, north
import aplpy
import pylab as pl
pl.matplotlib.rc_file('pubfiguresrc')

from get_m0 import get_mom0

#'CH3OH808-716':'CH$_3$OH $8_{0,8}-7_{1,6}$',

fn = paths.dpath('merge/W51_b6_7M_12M.CH3OH808-716.image.pbcor.fits')
cube = SpectralCube.read(fn)
m0 = get_mom0(fn)
mx = cube.max(axis=0)
avbm = cube._average_beams(1, mask=None)
mx_K = mx.to(u.K, u.brightness_temperature(avbm,
                                           cube.with_spectral_unit(u.GHz).spectral_axis.mean()))


blue_fits_fn = paths.dpath('moments/w51_12co2-1_blue0to45_masked.fits')
red_fits_fn = paths.dpath('moments/w51_12co2-1_red73to130_masked.fits')
red_fits = fits.open(red_fits_fn)
blue_fits = fits.open(blue_fits_fn)

redhead = red_fits[0].header
bluehead = blue_fits[0].header
red_wcs = wcs.WCS(redhead)
blue_wcs = wcs.WCS(bluehead)


#ax = fig1.add_axes([0.15, 0.1, 0.8, 0.8], projection=cube.wcs.celestial)
#ax.imshow(cutout_m0.data, cmap='gray')
#
#ax.contourf(cutout_blue.data,
#            transform=ax.get_transform(WCS(bluehead)),
#            levels=[0.5,1,2,3,4,5,6,7,8,9,10],
#            colors=[(0,0,1,ii) for ii in np.arange(0, 1, 0.1)])
#ax.contourf(cutout_red.data,
#            transform=ax.get_transform(WCS(redhead)),
#            levels=[0.5,1,2,3,4,5,6,7,8,9,10],
#            colors=[(1,0,0,ii) for ii in np.arange(0, 1, 0.1)])


for name, center, size in (('e2e', e2e, 7.5*u.arcsec),
                           ('between_e2e_and_e8', between_e2e_and_e8, 15*u.arcsec),
                           ('e8fil', e8fil, 12.5*u.arcsec),
                           ('north', north, 7.5*u.arcsec)):

    cutout_red = Cutout2D(red_fits[0].data, center, size, wcs=red_wcs)
    cutout_blue = Cutout2D(blue_fits[0].data, center, size, wcs=blue_wcs)

    cutout_m0 = Cutout2D(m0.value, center, size, wcs=cube.wcs.celestial)
    cutout_mx_K = Cutout2D(mx_K.value, center, size, wcs=cube.wcs.celestial)

    fig1 = pl.figure(1)
    fig1.clf()
    cutout_m0_hdu = fits.PrimaryHDU(data=cutout_m0.data, header=cutout_m0.wcs.to_header())
    cutout_mx_K_hdu = fits.PrimaryHDU(data=cutout_mx_K.data, header=cutout_mx_K.wcs.to_header())
    #FF = aplpy.FITSFigure(cutout_m0_hdu, figure=fig1)
    FF = aplpy.FITSFigure(cutout_mx_K_hdu, figure=fig1)
    FF.show_grayscale()
    FF.save(paths.fpath("methanol_{0}.png".format(name)))
    FF.set_tick_xspacing(1.8/3600.)
    
    levels = [0.5,1,2,3,4,5,6,7,8,9,10,15]
    alphas = np.linspace(0.1, 0.9, len(levels))

    FF.add_scalebar((0.1*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
    FF.scalebar.set_label('0.1 pc')
    FF.scalebar.set_color('w')

    FF.save(paths.fpath("outflows/methanol_scalebar_{0}.png".format(name)))

    FF.show_contour(fits.PrimaryHDU(data=cutout_blue.data, header=cutout_blue.wcs.to_header()),
                    levels=levels,
                    filled=True,
                    colors=[(0,0,1,ii) for ii in alphas])
    FF.show_contour(fits.PrimaryHDU(data=cutout_red.data, header=cutout_red.wcs.to_header()),
                    filled=True,
                    levels=levels,
                    colors=[(1,0,0,ii) for ii in alphas])

    FF.savefig(paths.fpath("outflows/outflows_over_methanol_{0}.png".format(name)))
