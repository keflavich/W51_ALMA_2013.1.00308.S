from astropy import units as u
import pylab as pl
from spectral_cube import SpectralCube
import aplpy
import os
import paths

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

rgb_cube_png = rgb_cube_fits[:-5]+"_loggreen_max.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              pmax_g=99.99,
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


fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_naco_png, figure=fig1)
F.show_rgb(rgb_cube_naco_png)
F.recenter(290.91665, 14.518691, radius=0.003)
F.add_scalebar((0.1*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("NACO_green_outflows_aplpy.png"))
F.save(paths.fpath("NACO_green_outflows_aplpy.pdf"))

F.show_contour(paths.dpath('12m/continuum/w51_spw3_continuum.image.fits'),
               levels=[0.15, 0.30, 0.45, 0.60, 0.75],
               colors=['w']*4, layer='alma_cont_lores')
F.save(paths.fpath("NACO_green_outflows_aplpy_CONTours.png"))
F.save(paths.fpath("NACO_green_outflows_aplpy_CONTours.pdf"))

F.hide_layer('alma_cont_lores')
F.show_contour(paths.dpath('W51_te_continuum_best.fits'),
               levels=[0.04,
                       #0.05,
                       0.06,
                       #0.08,
                       0.10,
                       #0.15,
                       0.20,
                       #0.30,
                       0.40,], colors=['w']*12,
               layer='alma_cont_hires')
F.save(paths.fpath("NACO_green_outflows_aplpy_CONTours_hires.png"))
F.save(paths.fpath("NACO_green_outflows_aplpy_CONTours_hires.pdf"))



h77a = SpectralCube.read(paths.vpath('data/W51north_H77_Outflow_cutout.fits'))
h77a_outflow = h77a.spectral_slab(-16*u.km/u.s, -60*u.km/u.s).sum(axis=0)
h77a_green = paths.dpath('W51_H77a_LacyJetOutflow_Sum.fits')
h77a_outflow.hdu.writeto(h77a_green, clobber=True)

F.show_contour(h77a_green, levels=[0.0075, 0.015, 0.030], colors=['b']*6,
               layer='h77a_outflow')
F.save(paths.fpath("NACO_green_outflows_aplpy_CONTours_hires_h77acontour.png"))
F.save(paths.fpath("NACO_green_outflows_aplpy_CONTours_hires_h77acontour.pdf"))


F.hide_layer('h77a_outflow')
F.hide_layer('alma_cont_hires')


F.show_contour(paths.dpath('longbaseline/W51ncax.cont.image.pbcor.fits'),
               levels=[0.0015, 0.0030, 0.0045, 0.0060, 0.0075],
               colors=['w']*4, layer='alma_cont_cycle3hires')
F.save(paths.fpath("NACO_green_outflows_aplpy_CONTours_cycle3hires.png"))
F.save(paths.fpath("NACO_green_outflows_aplpy_CONTours_cycle3hires.pdf"))


rgb_cube_h77a_fits = 'outflow_co_redblue_h77a_green.fits'
if not os.path.exists(rgb_cube_h77a_fits):
    aplpy.make_rgb_cube([red_fits, h77a_green, blue_fits], rgb_cube_h77a_fits)

rgb_cube_h77a_png = rgb_cube_h77a_fits[:-5]+"_asinhgreen.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_h77a_fits,
                              output=rgb_cube_h77a_png,
                              vmin_g=0.003,
                              vmax_g=0.03,
                              vmax_r=4.0,
                              vmin_r=0.00,
                              vmax_b=2.0,
                              vmin_b=0.00,
                              embed_avm_tags=True)

fig2 = pl.figure(2)
fig2.clf()
F = aplpy.FITSFigure(rgb_cube_h77a_png, figure=fig2)
F.show_rgb(rgb_cube_h77a_png)
F.recenter(290.91665, 14.518691, radius=0.005)
F.add_scalebar((0.1*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("H77a_green_outflows_aplpy.png"))
F.save(paths.fpath("H77a_green_outflows_aplpy.pdf"))



#green_fits = '/Users/adam/work/w51/alma/FITS/w51_spw3_continuum.image.fits'
green_fits = paths.dpath('W51_te_continuum_best.fits')
blue_fits = paths.dpath('moments/w51_12co2-1_blue0to45_masked.fits')
red_fits = paths.dpath('moments/w51_12co2-1_red73to130_masked.fits')
rgb_cube_fits = 'outflow_co_redblue_1.4mmcont_green.fits'

if not os.path.exists(rgb_cube_fits):
    # does not return anything
    aplpy.make_rgb_cube([red_fits, green_fits, blue_fits], rgb_cube_fits)

rgb_cube_png = rgb_cube_fits[:-5]+"_asinhgreen.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              vmax_g=0.6,
                              vmin_g=-0.005,
                              vmax_r=3,
                              vmin_r=-0.05,
                              vmax_b=3,
                              vmin_b=-0.05,
                              stretch_g='arcsinh', embed_avm_tags=True)


fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_png, figure=fig1)
F.show_rgb(rgb_cube_png)
F.recenter(290.93291, 14.508188, radius=0.003)
F.add_scalebar((0.1*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('w')
F.save(paths.fpath("ALMAcont_green_outflows_e2_aplpy.png"))
F.save(paths.fpath("ALMAcont_green_outflows_e2_aplpy.pdf"))

F.recenter(290.917, 14.518211, radius=0.005)
F.save(paths.fpath("ALMAcont_green_outflows_north_aplpy.png"))
F.save(paths.fpath("ALMAcont_green_outflows_north_aplpy.pdf"))

F.show_contour('/Users/adam/work/w51/paper_w51_evla/data/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits',
               levels=[0.005/9., 0.005/3., 0.005, 0.015, 0.045, 0.135],
               colors=['w']*7, layer='evla_cont_lores')
F.save(paths.fpath("Alma1.4mmcont_green_outflows_north_aplpy_cmCONTours.png"))
F.save(paths.fpath("Alma1.4mmcont_green_outflows_north_aplpy_cmCONTours.pdf"))

F.recenter(290.93291, 14.508188, radius=0.003)
F.save(paths.fpath("Alma1.4mmcont_green_outflows_e2_aplpy_cmCONTours.png"))
F.save(paths.fpath("Alma1.4mmcont_green_outflows_e2_aplpy_cmCONTours.pdf"))

