"""
REQUIRES aplpy branch my_master_mar2016
"""
import os
import aplpy
import paths
from vla_cont_cutout import fnku, fitsKu_fn

fn303 = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO303_202.image.pbcor_medsub_max.fits')
fn321 = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO321_220.image.pbcor_medsub_max.fits')
fn322 = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO322_221.coarser.image.pbcor_medsub_max.fits')
fnc18o = paths.dpath('merge/moments/W51_b6_7M_12M.C18O2-1.image.pbcor_medsub_max.fits')
fnSO = paths.dpath('merge/moments/W51_b6_7M_12M.SO65-54.image.pbcor_medsub_max.fits')
fnhc3n = paths.dpath('merge/moments/W51_b6_7M_12M.HC3N24-23.image.pbcor_medsub_max.fits')
fnch3oh422 = paths.dpath('merge/moments/W51_b6_7M_12M.CH3OH422-312.image.pbcor_medsub_max.fits')
fnSOFIA37 = paths.root+"../sofia/W51A_37um_single.fits"
fnSOFIA20 = paths.root+"../sofia/W51A_20um_single.fits"

def make_rgb(outname, redline='H2CO303_202', greenline='H2CO321_220',
             blueline='H2CO322_221',
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


make_rgb('W51_RGB_VLAKu_SOFIA37_SOFIA20.fits',
         redline=fnku,
         greenline=fnSOFIA37,
         blueline=fnSOFIA20,
         pmax_b=99.9,
         pmax_g=99.9,
         stretch_r='log',
         stretch_g='log',
         stretch_b='log',
         suffix='_log99',
         pmin_r=1,
        )
