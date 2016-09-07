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
import paths
from outflow_meta import e2e

e2_green_fits = '/Users/adam/work/w51/alma/cycle3goddi/W51e2.cont.image.pbcor.fits'
north_green_fits = '/Users/adam/work/w51/alma/cycle3goddi/W51n.cont.image.pbcor.fits'
blue_fits_fn = '/Users/adam/work/w51/alma/FITS/moments/longbaseline/SiO_m32to55kms_e2.fits'
red_fits_fn = '/Users/adam/work/w51/alma/FITS/moments/longbaseline/SiO_74to118kms_e2.fits'
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


rgb_cube_fits = 'e2e_outflow_SiO_redblue_cycle3green.fits'
if not os.path.exists(rgb_cube_fits):
    # does not return anything
    aplpy.make_rgb_cube([red_fits_cutoute2e_fn, e2_green_fits, blue_fits_cutoute2e_fn], rgb_cube_fits)

rgb_cube_png = rgb_cube_fits[:-5]+"_auto.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              vmax_g=0.017,
                              vmax_b=0.3,
                              vmax_r=0.6,
                              vmin_g=0.0001,
                              embed_avm_tags=True)

rgb_cube_png = rgb_cube_fits[:-5]+"_loggreen.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              vmax_g=0.017,
                              vmax_b=0.3,
                              vmax_r=0.6,
                              vmin_g=0.0001,
                              stretch_g='log', embed_avm_tags=True)

rgb_cube_png = rgb_cube_fits[:-5]+"_asinhgreen.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              vmax_g=0.017,
                              vmax_b=0.3,
                              vmax_r=0.6,
                              vmin_g=0.0001,
                              stretch_g='arcsinh', embed_avm_tags=True)


pl.rcParams['font.size'] = 18
fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_png, figure=fig1)
F.show_rgb(rgb_cube_png)
F.recenter(290.93315, 14.509584, radius=0.0005)
F.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('5000 au / 0.025 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_aplpy.png"))
F.save(paths.fpath("outflows/W51e2_cycle3green_SiO_outflows_aplpy.pdf"))

e8_rgb_cube_png = rgb_cube_fits[:-5]+"_asinhgreen.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=e8_rgb_cube_png,
                              vmax_g=0.017,
                              vmin_g=-0.0001,
                              vmax_b=0.3,
                              vmin_b=-0.05,
                              vmax_r=0.3,
                              vmin_r=-0.05,
                              stretch_g='arcsinh', embed_avm_tags=True)

fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_png, figure=fig1)
F.show_rgb(e8_rgb_cube_png)
F.recenter(290.93288, 14.507858, radius=0.00025)
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

# zoom-in figures
rgb_cube_naco_fits = 'outflow_co_redblue_naco_green.fits'
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
F.show_contour(paths.lbpath('W51n.cont.image.pbcor.fits'), levels=[0.001,
                                                                   0.0020,
                                                                   0.004,
                                                                   0.008,
                                                                   0.012],
               colors=['w']*6, layer='alma_cont_cycle3hires')
F.save(paths.fpath("outflows/NACO_green_outflows_aplpy_zoomNorth_cycle3hires.png"))
F.save(paths.fpath("outflows/NACO_green_outflows_aplpy_zoomNorth_cycle3hires.pdf"))
F.show_contour(h77a_green, levels=[0.0075, 0.015], colors=['b']*6,
               layer='h77a_outflow')
F.save(paths.fpath("outflows/NACO_green_outflows_aplpy_zoomNorth_cycle3hires_h77acontour.png"))
F.save(paths.fpath("outflows/NACO_green_outflows_aplpy_zoomNorth_cycle3hires_h77acontour.pdf"))
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
F.show_contour(paths.lbpath('W51n.cont.image.pbcor.fits'), levels=[0.001,
                                                                   0.0020,
                                                                   0.004,
                                                                   0.008,
                                                                   0.012],
               colors=['w']*6, layer='alma_cont_cycle3hires')
F.recenter(290.91564, 14.518128, radius=0.0005)
F.save(paths.fpath("outflows/NACO_green_outflows_aplpy_zoomALMAmm31_cycle3hires.png"))
F.save(paths.fpath("outflows/NACO_green_outflows_aplpy_zoomALMAmm31_cycle3hires.pdf"))
