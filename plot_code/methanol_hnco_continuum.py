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
from outflow_meta import e2e, e8fil, north, e8, e8south

from get_m0 import get_mom0

fn = paths.dpath('merge/W51_b6_7M_12M.CH3OH808-716.image.pbcor.fits')
m0ch3oh = get_mom0(fn, iterate=False) # iteration takes longer but doesn't risk eating 100% memory

cont_fits = fits.open(paths.dpath('w51_te_continuum_best.fits'))


def zoomfigure(target=e2e, targetname='e2e', radius=7.5*u.arcsec, cutout='e2e8'):

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
    F.recenter(target.ra.deg, target.dec.deg, radius=0.00075)
    F.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
    F.scalebar.set_label('5000 au / 0.025 pc')
    F.scalebar.set_color('w')
    F.set_tick_xspacing(0.0005)
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

zoomfigure(target=e2e, targetname='e2e', radius=7.5*u.arcsec, cutout='e2e8')
zoomfigure(target=e8fil, targetname='e8fil', radius=15*u.arcsec, cutout='e2e8')
zoomfigure(target=e8, targetname='e8', radius=8.5*u.arcsec, cutout='e2e8')
zoomfigure(target=e8south, targetname='e8south', radius=10*u.arcsec, cutout='e2e8')
zoomfigure(target=north, targetname='north', radius=10*u.arcsec, cutout='north')

#F.recenter(e8fil.ra.deg, e8fil.dec.deg, radius=0.001)
#F.save(paths.fpath("W51e8_ch3oh_hnco_continuum_aplpy_kucontours.png"))
#F.save(paths.fpath("W51e8_ch3oh_hnco_continuum_aplpy_kucontours.pdf"))
#
#F.recenter(north.ra.deg, north.dec.deg, radius=0.001)
#F.save(paths.fpath("W51north_ch3oh_hnco_continuum_aplpy_kucontours.png"))
#F.save(paths.fpath("W51north_ch3oh_hnco_continuum_aplpy_kucontours.pdf"))

#F.save(paths.fpath("W51e8_cycle3green_outflows_aplpy_cmcontours.png"))
#F.save(paths.fpath("W51e8_cycle3green_outflows_aplpy_cmcontours.pdf"))
#
## zoom-in figures
#rgb_cube_naco_fits = 'outflow_co_redblue_naco_green.fits'
#rgb_cube_naco_png_fullstretch = rgb_cube_naco_fits[:-5]+"_asinhgreen_fullstretch.png"
#rgb_im = aplpy.make_rgb_image(data=rgb_cube_naco_fits,
#                              output=rgb_cube_naco_png_fullstretch, stretch_g='arcsinh',
#                              vmin_g=-0.1,
#                              vmax_g=200,
#                              vmin_r=0.0,
#                              vmax_r=9,
#                              vmax_b=7,
#                              vmin_b=0.0,
#                              embed_avm_tags=True)
#fig1 = pl.figure(1)
#fig1.clf()
#F = aplpy.FITSFigure(rgb_cube_naco_png_fullstretch, figure=fig1)
#F.recenter(290.91689, 14.518196, radius=0.0005)
#F.show_rgb(rgb_cube_naco_png_fullstretch)
#F.add_scalebar((5000*u.au / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
#F.scalebar.set_label('5000 au / 0.025 pc')
#F.scalebar.set_color('w')
#F.show_contour('../../alma/cycle3goddi/W51n.cont.image.pbcor.fits', levels=[0.001,
#                                                                            0.0020,
#                                                                            0.004,
#                                                                            0.008,
#                                                                            0.012],
#               colors=['w']*6, layer='alma_cont_cycle3hires')
#F.save(paths.fpath("NACO_green_outflows_aplpy_zoomNorth_cycle3hires.png"))
#F.save(paths.fpath("NACO_green_outflows_aplpy_zoomNorth_cycle3hires.pdf"))
#F.show_contour(h77a_green, levels=[0.0075, 0.015], colors=['b']*6,
#               layer='h77a_outflow')
#F.save(paths.fpath("NACO_green_outflows_aplpy_zoomNorth_cycle3hires_h77acontour.png"))
#F.save(paths.fpath("NACO_green_outflows_aplpy_zoomNorth_cycle3hires_h77acontour.pdf"))
#F.remove_layer('h77a_outflow')
#
#rgb_cube_naco_png_fullstretch = rgb_cube_naco_fits[:-5]+"_asinhgreen_fullstretch_ALMAmm31.png"
#rgb_im = aplpy.make_rgb_image(data=rgb_cube_naco_fits,
#                              output=rgb_cube_naco_png_fullstretch, stretch_g='arcsinh',
#                              vmin_g=-0.1,
#                              vmax_g=50,
#                              vmin_r=-0.1,
#                              vmax_r=5.5,
#                              vmax_b=3,
#                              vmin_b=-0.1,
#                              embed_avm_tags=True)
#fig1.clf()
#F = aplpy.FITSFigure(rgb_cube_naco_png_fullstretch, figure=fig1)
#F.recenter(290.91689, 14.518196, radius=0.0005)
#F.show_rgb(rgb_cube_naco_png_fullstretch)
#F.add_scalebar((5000*u.au / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
#F.scalebar.set_label('5000 au / 0.025 pc')
#F.scalebar.set_color('w')
#F.show_contour('../../alma/cycle3goddi/W51n.cont.image.pbcor.fits', levels=[0.001,
#                                                                            0.0020,
#                                                                            0.004,
#                                                                            0.008,
#                                                                            0.012],
#               colors=['w']*6, layer='alma_cont_cycle3hires')
#F.recenter(290.91564, 14.518128, radius=0.0005)
#F.save(paths.fpath("NACO_green_outflows_aplpy_zoomALMAmm31_cycle3hires.png"))
#F.save(paths.fpath("NACO_green_outflows_aplpy_zoomALMAmm31_cycle3hires.pdf"))
#
