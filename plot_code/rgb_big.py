"""
REQUIRES aplpy branch my_master_mar2016
"""
#import numpy as np
#from astropy import coordinates
#from astropy import units as u
import os
import aplpy
import paths
#from outflow_meta import e2e
from vla_cont_cutout import fnku, fitsKu_fn

fn303 = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO303_202.image.pbcor_medsub_max.fits')
fn321 = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO321_220.image.pbcor_medsub_max.fits')
fn322 = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO322_221.coarser.image.pbcor_medsub_max.fits')
fnc18o = paths.dpath('merge/moments/W51_b6_7M_12M.C18O2-1.image.pbcor_medsub_max.fits')
fnSO = paths.dpath('merge/moments/W51_b6_7M_12M.SO65-54.image.pbcor_medsub_max.fits')
fnhc3n = paths.dpath('merge/moments/W51_b6_7M_12M.HC3N24-23.image.pbcor_medsub_max.fits')
fnch3oh422 = paths.dpath('merge/moments/W51_b6_7M_12M.CH3OH422-312.image.pbcor_medsub_max.fits')
#fits303 = fits.open(fn303)
#fits321 = fits.open(fn321)
#fits322 = fits.open(fn322)

def make_rgb(outname, redline='H2CO303_202', greenline='H2CO321_220', blueline='H2CO322_221',
             fntemplate=paths.dpath('merge/moments/W51_b6_7M_12M.{0}.image.pbcor_medsub_max.fits'),
             suffix="_auto",
             **kwargs):

    print(outname, suffix)
    rgb_cube_fits = outname
    if not os.path.exists(rgb_cube_fits):
        # does not return anything
        aplpy.make_rgb_cube([fntemplate.format(redline) if 'fits' not in redline else redline,
                             fntemplate.format(greenline) if 'fits' not in greenline else greenline,
                             fntemplate.format(blueline) if 'fits' not in blueline else blueline,],
                            rgb_cube_fits)

    rgb_cube_png = rgb_cube_fits[:-5]+suffix+".png"
    rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                                  embed_avm_tags=True, **kwargs)
    return rgb_im

for suffix, pct in (('_auto', 99.75),
                    ('_99.99', 99.99),
                   ):
    make_rgb('full_h2co_12monly_rgb.fits',
             fntemplate=paths.dpath('12m/moments/W51_b6_12M.{0}.image.pbcor_medsub_max.fits'),
             suffix=suffix,
             pmax_g=pct, pmax_b=pct, pmax_r=pct,
            )
    make_rgb('h2co_hc3n_ch3oh_rgb.fits',
             blueline='HC3N24-23',
             greenline='CH3OH422-312',
             suffix=suffix,
             pmax_g=pct, pmax_b=pct, pmax_r=pct,
            )
    make_rgb('h2co_so_hc3n_rgb.fits',
             blueline='HC3N24-23',
             greenline='SO65-54',
             suffix=suffix,
             pmax_g=pct, pmax_b=pct, pmax_r=pct,
            )
    make_rgb('c18o_hc3n_ch3oh_rgb.fits',
             redline='C18O2-1',
             blueline='HC3N24-23',
             greenline='CH3OH422-312',
             suffix=suffix,
             pmax_g=pct, pmax_b=pct, pmax_r=pct,
            )
    make_rgb('c18o_hc3n_so_rgb.fits',
             redline='C18O2-1',
             blueline='HC3N24-23',
             greenline='SO65-54',
             suffix=suffix,
             pmax_g=pct, pmax_b=pct, pmax_r=pct,
            )
    make_rgb('ku_hc3n_ch3oh_rgb.fits',
             redline=fnku,
             blueline='HC3N24-23',
             greenline='CH3OH422-312',
             suffix=suffix,
             pmax_g=pct, pmax_b=pct, pmax_r=pct,
            )
    make_rgb('hc3n_ch3oh_ocs_rgb.fits',
             greenline='OCS18-17',
             redline='HC3N24-23',
             blueline='CH3OH422-312',
             suffix=suffix,
             pmax_g=pct, pmax_b=pct, pmax_r=pct,
            )

make_rgb('ku_hc3n_ch3oh_rgb.fits',
         suffix="_max",
         pmax_r=99.75,
         pmax_b=99.9,
         pmax_g=99.9999,
         redline=fnku,
         blueline='HC3N24-23',
         greenline='CH3OH422-312')
make_rgb('ku_so_c18o_rgb.fits',
         greenline='SO65-54',
         redline=fnku,
         blueline='C18O2-1',
         pmax_g=99.99,
         pmax_r=99.95,
         pmax_b=99.99,
        )


rgb_cube_fits = 'full_h2co_rgb.fits'
if not os.path.exists(rgb_cube_fits):
    # does not return anything
    aplpy.make_rgb_cube([fn303, fn321, fn322,], rgb_cube_fits)

rgb_cube_png = rgb_cube_fits[:-5]+"_auto.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              embed_avm_tags=True)

rgb_cube_png = rgb_cube_fits[:-5]+"_setlevels.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              vmin_b=-0.005,
                              vmax_b=0.4,
                              vmin_g=-0.005,
                              vmax_g=0.4,
                              vmin_r=-0.005,
                              vmax_r=0.4,
                              embed_avm_tags=True)

rgb_cube_fits = 'c18o_h2co_ku_rgb.fits'
if not os.path.exists(rgb_cube_fits):
    # does not return anything
    aplpy.make_rgb_cube([fitsKu_fn, fn303, fnc18o,], rgb_cube_fits)

rgb_cube_png = rgb_cube_fits[:-5]+"_auto.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              embed_avm_tags=True)

rgb_cube_png = rgb_cube_fits[:-5]+"_setlevels.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              vmin_b=-0.05,
                              vmax_b=0.65,
                              vmin_g=-0.005,
                              vmax_g=0.4,
                              vmin_r=-0.0005,
                              vmax_r=0.1,
                              embed_avm_tags=True)

rgb_cube_png = rgb_cube_fits[:-5]+"_logred.png"
rgb_im = aplpy.make_rgb_image(data=rgb_cube_fits, output=rgb_cube_png,
                              vmin_b=-0.05,
                              vmax_b=0.65,
                              vmin_g=-0.005,
                              vmax_g=0.4,
                              vmin_r=-0.0005,
                              vmax_r=0.01,
                              stretch_r='log',
                              embed_avm_tags=True)
