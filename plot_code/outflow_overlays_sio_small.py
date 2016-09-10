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
from outflow_meta import e2e

e2_green_fits = paths.dpath('longbaseline/W51e2cax.cont.image.pbcor.fits')
north_green_fits = paths.dpath('longbaseline/W51ncax.cont.image.pbcor.fits')
blue_fits_fn = paths.dpath('longbaseline/SiO_m32to55kms_e2.fits')
red_fits_fn = paths.dpath('longbaseline/SiO_74to118kms_e2.fits')
#h77a
h77a_green = paths.dpath('W51_H77a_LacyJetOutflow_Sum.fits')

red_fits = fits.open(red_fits_fn)
blue_fits = fits.open(blue_fits_fn)
redhead = red_fits[0].header
bluehead = blue_fits[0].header
red_wcs = wcs.WCS(redhead)
blue_wcs = wcs.WCS(bluehead)

cutout_red = Cutout2D(red_fits[0].data, e2e, 0.0072*u.deg, wcs=red_wcs)
cutout_blue = Cutout2D(blue_fits[0].data, e2e, 0.0072*u.deg, wcs=blue_wcs)

red_fits_cutoute2e_fn = '/Users/adam/work/w51/alma/FITS/moments/w51_LB_SiO_red74to118_masked_cutoute2e.fits'
redhead.update(cutout_red.wcs.to_header())
red_fits_co = fits.PrimaryHDU(data=cutout_red.data, header=redhead)
red_fits_co.writeto(red_fits_cutoute2e_fn, clobber=True)

blue_fits_cutoute2e_fn = '/Users/adam/work/w51/alma/FITS/moments/w51_LB_SiO_bluem32to55_masked_cutoute2e.fits'
bluehead.update(cutout_blue.wcs.to_header())
blue_fits_co = fits.PrimaryHDU(data=cutout_blue.data, header=bluehead)
blue_fits_co.writeto(blue_fits_cutoute2e_fn, clobber=True)


e2_rgb_cube_fits = 'e2e_outflow_SiO_redblue_cycle3green.fits'
if not os.path.exists(e2_rgb_cube_fits):
    # does not return anything
    aplpy.make_rgb_cube([red_fits_cutoute2e_fn, e2_green_fits, blue_fits_cutoute2e_fn], e2_rgb_cube_fits)

rgb_cube_png = e2_rgb_cube_fits[:-5]+"_auto.png"
rgb_im = aplpy.make_rgb_image(data=e2_rgb_cube_fits, output=rgb_cube_png,
                              vmax_g=0.017,
                              vmax_b=0.3,
                              vmax_r=0.6,
                              vmin_g=0.0001,
                              embed_avm_tags=True)

rgb_cube_png = e2_rgb_cube_fits[:-5]+"_loggreen.png"
rgb_im = aplpy.make_rgb_image(data=e2_rgb_cube_fits, output=rgb_cube_png,
                              vmax_g=0.017,
                              vmax_b=0.3,
                              vmax_r=0.6,
                              vmin_g=0.0001,
                              stretch_g='log', embed_avm_tags=True)

rgb_cube_png = e2_rgb_cube_fits[:-5]+"_asinhgreen.png"
rgb_im = aplpy.make_rgb_image(data=e2_rgb_cube_fits, output=rgb_cube_png,
                              vmax_g=0.017,
                              vmax_b=0.5,
                              vmax_r=0.2,
                              vmin_r=-0.0025,
                              vmin_b=-0.0025,
                              vmin_g=0.0001,
                              stretch_g='arcsinh', embed_avm_tags=True)

rgb_cube_png_faintercont = e2_rgb_cube_fits[:-5]+"_lineargreen.png"
rgb_im = aplpy.make_rgb_image(data=e2_rgb_cube_fits, output=rgb_cube_png_faintercont,
                              vmax_g=0.017,
                              vmax_b=0.6,
                              vmax_r=0.2,
                              vmin_r=-0.0025,
                              vmin_b=-0.0025,
                              vmin_g=0.0001,
                              stretch_g='linear', embed_avm_tags=True)


pl.rcParams['font.size'] = 18
fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_png, figure=fig1)
F.show_rgb(rgb_cube_png)
F.recenter(290.93315, 14.509584, radius=0.0004)
F.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('5000 au / 0.025 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_aplpy.png"))
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_aplpy.pdf"))


F.recenter(290.93318, 14.509591, radius=0.0001)
F.scalebar.set_length((0.005*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('1000 au / 0.005 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_aplpy_zoom.png"))
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_aplpy_zoom.pdf"))


fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_png_faintercont, figure=fig1)
F.show_rgb(rgb_cube_png_faintercont)
F.recenter(290.93315, 14.509584, radius=0.0004)
F.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('5000 au / 0.025 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_faintercontinuum_aplpy.png"))
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_faintercontinuum_aplpy.pdf"))


F.recenter(290.93318, 14.509591, radius=0.0001)
F.scalebar.set_length((0.005*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('1000 au / 0.005 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_faintercontinuum_aplpy_zoom.png"))
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_faintercontinuum_aplpy_zoom.pdf"))





blue_fits_e8_fn = paths.dpath('longbaseline/SiO_m32to55kms_e8.fits')
red_fits_e8_fn = paths.dpath('longbaseline/SiO_74to118kms_e8.fits')

e8_rgb_cube_fits = 'e8_outflow_SiO_redblue_cycle3green.fits'
if not os.path.exists(e8_rgb_cube_fits):
    # does not return anything
    aplpy.make_rgb_cube([red_fits_e8_fn, e2_green_fits, blue_fits_e8_fn], e8_rgb_cube_fits)


e8_rgb_cube_png = e8_rgb_cube_fits[:-5]+"_asinhgreen.png"
rgb_im = aplpy.make_rgb_image(data=e8_rgb_cube_fits, output=e8_rgb_cube_png,
                              vmax_g=0.019,
                              vmin_g=-0.0001,
                              vmax_b=0.3,
                              vmin_b=-0.025,
                              vmax_r=0.2,
                              vmin_r=-0.025,
                              stretch_g='arcsinh', embed_avm_tags=True)

fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_png, figure=fig1)
F.show_rgb(e8_rgb_cube_png)
F.recenter(290.93291, 14.507857, radius=0.00025)
F.add_scalebar((0.01*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('2000 au / 0.01 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("outflows/W51e8_cycle3green_SiO_outflows_aplpy.png"))
F.save(paths.fpath("outflows/W51e8_cycle3green_SiO_outflows_aplpy.pdf"))

F.show_contour('/Users/adam/work/w51/paper_w51_evla/data/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits',
               levels=np.array([0.0015,0.0045,0.0135,0.0270])*1.25,
               colors=['w']*7, layer='evla_cont_lores')

F.save(paths.fpath("outflows/W51e8_cycle3green_outflows_aplpy_cmcontours.png"))
F.save(paths.fpath("outflows/W51e8_cycle3green_outflows_aplpy_cmcontours.pdf"))


F.hide_layer('evla_cont_lores')
F.recenter(290.93291, 14.507857, radius=0.00010)
F.scalebar.set_length((0.005*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('1000 au / 0.005 pc')
F.save(paths.fpath("outflows/W51e8_cycle3green_SiO_outflows_aplpy_zoom.png"))
F.save(paths.fpath("outflows/W51e8_cycle3green_SiO_outflows_aplpy_zoom.pdf"))


north_blue_fits_fn = paths.dpath('longbaseline/SiO_m32to55kms_north.fits')
north_red_fits_fn = paths.dpath('longbaseline/SiO_74to118kms_north.fits')
north_rgb_cube_fits = 'north_outflow_SiO_redblue_cycle3green.fits'
if not os.path.exists(north_rgb_cube_fits):
    # does not return anything
    aplpy.make_rgb_cube([north_red_fits_fn, north_green_fits,
                         north_blue_fits_fn], north_rgb_cube_fits)

north_rgb_cube_png = north_rgb_cube_fits[:-5]+"_asinhgreen.png"
north_rgb_im = aplpy.make_rgb_image(data=north_rgb_cube_fits,
                                    output=north_rgb_cube_png, vmax_g=0.019,
                                    vmax_b=0.6, vmax_r=0.3, vmin_g=0.0001,
                                    stretch_g='arcsinh', embed_avm_tags=True)

fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(north_rgb_cube_png, figure=fig1)
F.show_rgb(north_rgb_cube_png)
F.add_scalebar((0.01*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('2000 au / 0.01 pc')
F.scalebar.set_color('w')
F.recenter(290.91688, 14.518189, radius=0.000284)
F.save(paths.fpath("outflows/W51north_cycle3green_SiO_outflows_aplpy.png"))
F.save(paths.fpath("outflows/W51north_cycle3green_SiO_outflows_aplpy.pdf"))

F.show_contour('/Users/adam/work/w51/paper_w51_evla/data/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits',
               levels=np.array([0.0015,0.0045,0.0135,0.0270])*1.25,
               colors=['w']*7, layer='evla_cont_lores')

F.save(paths.fpath("outflows/W51north_cycle3green_outflows_aplpy_cmcontours.png"))
F.save(paths.fpath("outflows/W51north_cycle3green_outflows_aplpy_cmcontours.pdf"))

F.hide_layer('evla_cont_lores')
F.scalebar.set_length((0.005*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('1000 au / 0.005 pc')
F.scalebar.set_color('w')
F.recenter(290.91688, 14.518189, radius=0.0001)
F.save(paths.fpath("outflows/W51north_cycle3green_SiO_outflows_aplpy_zoom.png"))
F.save(paths.fpath("outflows/W51north_cycle3green_SiO_outflows_aplpy_zoom.pdf"))





naco_green = '/Users/adam/work/w51/paper_w51_evla/data/naco_Kband_W51.fits'
rgb_cube_naco_fits = 'outflow_sio_redblue_naco_green.fits'
if not os.path.exists(rgb_cube_naco_fits):
    aplpy.make_rgb_cube([north_red_fits_fn, naco_green, north_blue_fits_fn],
                        rgb_cube_naco_fits)

rgb_cube_naco_png_fullstretch = rgb_cube_naco_fits[:-5]+"_asinhgreen_fullstretch.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_naco_fits,
                              output=rgb_cube_naco_png_fullstretch, stretch_g='arcsinh',
                              vmin_g=-0.1,
                              vmax_g=200,
                              vmin_r=0.0,
                              vmax_r=0.3,
                              vmax_b=0.3,
                              vmin_b=0.0,
                              embed_avm_tags=True)
fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_naco_png_fullstretch, figure=fig1)
F.recenter(290.91689, 14.518196, radius=0.0005)
F.show_rgb(rgb_cube_naco_png_fullstretch)
F.add_scalebar((5000*u.au / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('5000 au / 0.025 pc')
F.scalebar.set_color('w')
F.show_contour(paths.dpath('longbaseline/W51ncax.cont.image.pbcor.fits'),
               levels=[0.001, 0.0020, 0.004, 0.008, 0.012],
               colors=['w']*6, layer='alma_cont_cycle3hires')
F.save(paths.fpath("outflows/NACO_green_SiO_outflows_aplpy_zoomNorth_cycle3hires.png"))
F.save(paths.fpath("outflows/NACO_green_SiO_outflows_aplpy_zoomNorth_cycle3hires.pdf"))
F.show_contour(h77a_green, levels=[0.0075, 0.015], colors=['b']*6,
               layer='h77a_outflow')
F.save(paths.fpath("outflows/NACO_green_SiO_outflows_aplpy_zoomNorth_cycle3hires_h77acontour.png"))
F.save(paths.fpath("outflows/NACO_green_SiO_outflows_aplpy_zoomNorth_cycle3hires_h77acontour.pdf"))
F.remove_layer('h77a_outflow')

rgb_cube_naco_png_fullstretch = rgb_cube_naco_fits[:-5]+"_asinhgreen_fullstretch_ALMAmm31.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_naco_fits,
                              output=rgb_cube_naco_png_fullstretch, stretch_g='arcsinh',
                              vmin_g=-0.1,
                              vmax_g=50,
                              vmin_r=-0.1,
                              vmax_r=5.5,
                              vmax_b=3,
                              vmin_b=-0.1,
                              embed_avm_tags=True)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_naco_png_fullstretch, figure=fig1)
F.recenter(290.91689, 14.518196, radius=0.0005)
F.show_rgb(rgb_cube_naco_png_fullstretch)
F.add_scalebar((5000*u.au / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('5000 au / 0.025 pc')
F.scalebar.set_color('w')
F.show_contour(paths.dpath('longbaseline/W51ncax.cont.image.pbcor.fits'),
               levels=[0.001, 0.0020, 0.004, 0.008, 0.012],
               colors=['w']*6, layer='alma_cont_cycle3hires')
F.recenter(290.91564, 14.518128, radius=0.0005)
F.save(paths.fpath("outflows/NACO_green_SiO_outflows_aplpy_zoomALMAmm31_cycle3hires.png"))
F.save(paths.fpath("outflows/NACO_green_SiO_outflows_aplpy_zoomALMAmm31_cycle3hires.pdf"))



F.scalebar.set_length((0.005*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('1000 au / 0.005 pc')
F.scalebar.set_color('w')
F.recenter(290.91688, 14.518189, radius=0.0001)
F.save(paths.fpath("outflows/NACO_green_SiO_outflows_aplpy_zoomALMAmm31_cycle3hires_zoom.png"))
F.save(paths.fpath("outflows/NACO_green_SiO_outflows_aplpy_zoomALMAmm31_cycle3hires_zoom.pdf"))
