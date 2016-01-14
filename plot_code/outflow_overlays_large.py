from astropy import units as u
import aplpy
import os

green_fits = '/Users/adam/work/w51/paper_w51_evla/data/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'
blue_fits = '/Users/adam/work/w51/alma/FITS/moments/w51_12co2-1_blue0to45_masked.fits'
red_fits = '/Users/adam/work/w51/alma/FITS/moments/w51_12co2-1_red73to130_masked.fits'
rgb_cube_fits = 'outflow_co_redblue_kucont_green.fits'

if not os.path.exists(rgb_cube_fits):
    # does not return anything
    aplpy.make_rgb_cube([red_fits, green_fits, blue_fits], rgb_cube_fits)

rgb_cube_png = rgb_cube_fits[:-5]+"_auto.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              embed_avm_tags=True)

rgb_cube_png = rgb_cube_fits[:-5]+"_loggreen.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              stretch_g='log', embed_avm_tags=True)

rgb_cube_png = rgb_cube_fits[:-5]+"_asinhgreen.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              stretch_g='arcsinh', embed_avm_tags=True)

naco_green = '/Users/adam/work/w51/paper_w51_evla/data/naco_Kband_W51.fits'
rgb_cube_naco_fits = 'outflow_co_redblue_naco_green.fits'
if not os.path.exists(rgb_cube_naco_fits):
    aplpy.make_rgb_cube([red_fits, naco_green, blue_fits], rgb_cube_naco_fits)

rgb_cube_naco_png = rgb_cube_naco_fits[:-5]+"_asinhgreen.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_naco_fits,
                              output=rgb_cube_naco_png, stretch_g='arcsinh',
                              embed_avm_tags=True)

import pylab as pl
fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_naco_png, figure=fig1)
F.show_rgb(rgb_cube_naco_png)
F.recenter(290.91665, 14.518691, radius=0.003)
F.add_scalebar((0.1*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('w')
F.save("NACO_green_outflows_aplpy.png")

F.show_contour('../../alma/FITS/w51_spw3_continuum.image.fits', levels=[0.15,
                                                                        0.30,
                                                                        0.45,
                                                                        0.60,
                                                                        0.75],
               colors=['w']*4)
F.save("NACO_green_outflows_aplpy_CONTours.png")
