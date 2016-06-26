"""
REQUIRES aplpy branch my_master_mar2016
"""
import numpy as np
from astropy import units as u
import os
import pylab as pl
import aplpy
import paths
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs
from outflow_meta import e2e, e8fil, north, e8, e8south, d2

from get_m0 import get_mom0

fn = paths.dpath('merge/W51_b6_7M_12M.CH3OH808-716.image.pbcor.fits')
m0ch3oh = get_mom0(fn, iterate=False) # iteration takes longer but doesn't risk eating 100% memory

cont_fits = fits.open(paths.dpath('w51_te_continuum_best.fits'))


def zoomfigure(target=e2e, targetname='e2e', radius=7.5*u.arcsec, cutout='e2e8',
               zoom_radius=3*u.arcsec, tick_spacing=1.8*u.arcsec):

    fn = paths.dpath('merge/cutouts/W51_b6_7M_12M.HNCO10010-909.image.pbcor_{0}cutout.fits'.format(cutout))
    m0hnco = get_mom0(fn, iterate=False)

    cutout_cont = Cutout2D(cont_fits[0].data, target, radius, wcs=wcs.WCS(cont_fits[0].header))
    cutout_ch3oh = Cutout2D(m0ch3oh.value, target, radius, wcs=wcs.WCS(m0ch3oh.header))
    cutout_hnco = Cutout2D(m0hnco.value, target, radius, wcs=wcs.WCS(m0hnco.header))

    cont_fits_cutout = fits.PrimaryHDU(data=cutout_cont.data, header=cutout_cont.wcs.to_header())
    ch3oh_fits_cutout = fits.PrimaryHDU(data=cutout_ch3oh.data, header=cutout_ch3oh.wcs.to_header())
    hnco_fits_cutout = fits.PrimaryHDU(data=cutout_hnco.data, header=cutout_hnco.wcs.to_header())
    cont_fits_fn = "rgb/continuum_{0}_cutout.fits".format(targetname)
    hnco_fits_fn = "rgb/hnco_{0}_cutout.fits".format(targetname)
    ch3oh_fits_fn = "rgb/ch3oh_{0}_cutout.fits".format(targetname)
    cont_fits_cutout.writeto(cont_fits_fn, clobber=True)
    ch3oh_fits_cutout.writeto(ch3oh_fits_fn, clobber=True)
    hnco_fits_cutout.writeto(hnco_fits_fn, clobber=True)


    rgb_cube_fits = '{0}_ch3oh_hnco_cont.fits'.format(targetname)
    if not os.path.exists(rgb_cube_fits):
        # does not return anything
        aplpy.make_rgb_cube([ch3oh_fits_fn, hnco_fits_fn, cont_fits_fn,], rgb_cube_fits)

    rgb_cube_png = rgb_cube_fits[:-5]+"_auto.png"
    rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                                  embed_avm_tags=True)

    rgb_cube_png = rgb_cube_fits[:-5]+"_logcont.png"
    rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                                  #vmin_b=0.005,
                                  #vmax_b=0.15,
                                  stretch_b='log', embed_avm_tags=True)

    #rgb_cube_png = rgb_cube_fits[:-5]+"_asinhgreen.png"
    #rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
    #                              vmax_g=0.017,
    #                              vmax_b=6.5,
    #                              vmax_r=7.0,
    #                              vmin_g=0.0001,
    #                              stretch_g='arcsinh', embed_avm_tags=True)
    #
    #
    pl.rcParams['font.size'] = 18
    fig1 = pl.figure(1)
    fig1.clf()
    F = aplpy.FITSFigure(rgb_cube_png, figure=fig1)
    F.show_rgb(rgb_cube_png)
    #F.recenter(290.93315, 14.509584, radius=0.00075)
    F.recenter(target.ra.deg, target.dec.deg, radius=zoom_radius.to(u.deg).value)
    F.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
    F.scalebar.set_label('5000 au / 0.025 pc')
    F.scalebar.set_color('w')
    F.set_tick_xspacing(tick_spacing.to(u.deg).value)
    F.add_label(0.05, 0.95, "CH$_3$OH", relative=True, color='r', horizontalalignment='left')
    F.add_label(0.05, 0.91, "HNCO", relative=True, color='g', horizontalalignment='left')
    F.add_label(0.05, 0.87, "Continuum", relative=True, color='b', horizontalalignment='left')
    F.save(paths.fpath("W51{0}_ch3oh_hnco_continuum_aplpy.png".format(targetname)))
    F.save(paths.fpath("W51{0}_ch3oh_hnco_continuum_aplpy.pdf".format(targetname)))

    F.show_contour(paths.vpath('data/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'),
                   levels=np.array([0.0015,0.0045,0.0135,0.0270,0.054,0.108])*1.25,
                   colors=['w']*7, layer='evla_cont')
    F.save(paths.fpath("W51{0}_ch3oh_hnco_continuum_aplpy_kucontours.png".format(targetname)))
    F.save(paths.fpath("W51{0}_ch3oh_hnco_continuum_aplpy_kucontours.pdf".format(targetname)))

zoomfigure(target=e2e, targetname='e2e', radius=7.5*u.arcsec, cutout='e2e8',
           zoom_radius=3*u.arcsec, tick_spacing=1.8*u.arcsec)
zoomfigure(target=e8fil, targetname='e8fil', radius=15*u.arcsec, cutout='e2e8',
           zoom_radius=10*u.arcsec, tick_spacing=5*u.arcsec)
zoomfigure(target=e8, targetname='e8', radius=8.5*u.arcsec, cutout='e2e8',
           zoom_radius=4*u.arcsec, tick_spacing=2*u.arcsec)
zoomfigure(target=e8south, targetname='e8south', radius=15*u.arcsec,
           cutout='e2e8', zoom_radius=5*u.arcsec, tick_spacing=2*u.arcsec)
zoomfigure(target=north, targetname='north', radius=10*u.arcsec,
           cutout='north', zoom_radius=4*u.arcsec, tick_spacing=2*u.arcsec)
zoomfigure(target=d2, targetname='d2', radius=10*u.arcsec,
           cutout='north', zoom_radius=4*u.arcsec, tick_spacing=2*u.arcsec)
