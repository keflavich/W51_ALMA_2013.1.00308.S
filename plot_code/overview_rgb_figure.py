import numpy as np
from astropy import wcs
import aplpy
import paths
from astropy.io import fits
import reproject
from vla_cont_cutout import fnku
from matplotlib.colors import Normalize,LogNorm
from matplotlib.colors import rgb_to_hsv,hsv_to_rgb
import PIL
from PIL import ImageEnhance
import pyavm

if 'c18o' not in locals():
    fncont = paths.dpath('W51_te_continuum_best.fits')
    fnc18o = paths.dpath('merge/moments/W51_b6_7M_12M.C18O2-1.image.pbcor_medsub_max.fits')
    fnSO = paths.dpath('merge/moments/W51_b6_7M_12M.SO65-54.image.pbcor_medsub_max.fits')
    fnhc3n = paths.dpath('merge/moments/W51_b6_7M_12M.HC3N24-23.image.pbcor_medsub_max.fits')
    fnch3oh422 = paths.dpath('merge/moments/W51_b6_7M_12M.CH3OH422-312.image.pbcor_medsub_max.fits')
    fnch3oh808 = paths.dpath('merge/moments/W51_b6_7M_12M.CH3OH808-716.image.pbcor_medsub_moment0.fits')

    outhdr = fits.getheader(fnc18o)
    #c18o = reproject.reproject_interp(fnc18o, outhdr, order=1)[0]
    c18o = fits.getdata(fnc18o)
    mmcont = reproject.reproject_interp(fncont, outhdr, order=1)[0]
    SO = reproject.reproject_interp(fnSO, outhdr, order=1)[0]
    hc3n = reproject.reproject_interp(fnhc3n, outhdr, order=1)[0]
    ch3oh422 = reproject.reproject_interp(fnch3oh422, outhdr, order=1)[0]
    ch3oh808 = reproject.reproject_interp(fnch3oh808, outhdr, order=1)[0]

    ku_hdu = fits.open(fnku)
    ku_hdu[0].data = ku_hdu[0].data.squeeze()
    ku_hdu[0].header = wcs.WCS(ku_hdu[0].header).celestial.to_header()
    ku = reproject.reproject_interp(ku_hdu, outhdr, order=1)[0]




# BEGIN IMAGE MAKING HERE
red,green,blue,alpha = 0,1,2,3

# FIRST LAYER: C18O = Blue
rgb_im = np.zeros(c18o.shape + (4,))
rgb_im[:,:,alpha]=1.0
colorized_c18o = Normalize(vmin=np.nanpercentile(c18o,10),
                           vmax=np.nanpercentile(c18o,99.995), clip=True)(c18o)
rgb_im[:,:,blue] += np.nan_to_num(colorized_c18o)

# SECOND LAYER: SO = orangish?  ...
monochrome_SO = Normalize(vmin=np.nanpercentile(SO,10),
                          vmax=np.nanpercentile(SO,99.995), clip=True)(SO)
rgb_so = np.zeros(c18o.shape + (3,))
rgb_so[:,:,blue] = np.nan_to_num(monochrome_SO)
hsv_SO = rgb_to_hsv(rgb_so)
hue, saturation, value = 0,1,2
hsv_SO[:,:,hue] = 210/360.
rgb_so = hsv_to_rgb(hsv_SO)
#rgb_im[:,:,:alpha] += rgb_so

monochrome_hc3n = Normalize(vmin=np.nanpercentile(hc3n,10),
                            vmax=np.nanpercentile(hc3n,99.9995), clip=True)(hc3n)
rgb_hc3n = np.zeros(c18o.shape + (3,))
rgb_hc3n[:,:,blue] = np.nan_to_num(monochrome_hc3n)
hsv_hc3n = rgb_to_hsv(rgb_hc3n)
hsv_hc3n[:,:,hue] = 305/360.
rgb_hc3n = hsv_to_rgb(hsv_hc3n)
rgb_im[:,:,:alpha] += np.nan_to_num(rgb_hc3n)

monochrome_ch3oh422 = Normalize(vmin=np.nanpercentile(ch3oh422,10),
                                vmax=np.nanpercentile(ch3oh422,99.9995), clip=True)(ch3oh422)
rgb_ch3oh422 = np.zeros(c18o.shape + (3,))
#rgb_im[:,:,red] += np.nan_to_num(monochrome_ch3oh422)
rgb_ch3oh422[:,:,blue] = np.nan_to_num(monochrome_ch3oh422)
hsv_ch3oh422 = rgb_to_hsv(rgb_ch3oh422)
hsv_ch3oh422[:,:,hue] = 25/360.
rgb_ch3oh422 = hsv_to_rgb(hsv_ch3oh422)
rgb_im[:,:,:alpha] += rgb_ch3oh422

monochrome_ch3oh808 = Normalize(vmin=np.nanpercentile(ch3oh808,10),
                                vmax=np.nanpercentile(ch3oh808,99.9995), clip=True)(ch3oh808)
rgb_ch3oh808 = np.zeros(c18o.shape + (3,))
rgb_ch3oh808[:,:,blue] = np.nan_to_num(monochrome_ch3oh808)
hsv_ch3oh808 = rgb_to_hsv(rgb_ch3oh808)
hsv_ch3oh808[:,:,hue] = 45/360.
rgb_ch3oh808 = hsv_to_rgb(hsv_ch3oh808)
rgb_im[:,:,:alpha] += rgb_ch3oh808


monochrome_ku = LogNorm(vmin=0.0005-np.nanmin(ku),#np.nanpercentile(ku-np.nanmin(ku),30),
                        vmax=np.nanpercentile(ku-np.nanmin(ku),99.9995), clip=True)(ku-np.nanmin(ku))
rgb_ku = monochrome_ku[:,:,None]
rgb_im[:,:,:alpha] += np.nan_to_num(rgb_ku) * 0.75


vmin = 0.0025
monochrome_mmcont = LogNorm(vmin=vmin,
                            vmax=np.nanpercentile(mmcont,99.9995),
                            clip=True)(mmcont*(mmcont > vmin))
vmin = 0.001
monochrome_mmcont = Normalize(vmin=vmin,
                              vmax=np.nanpercentile(mmcont,99.9),
                              clip=True)(mmcont*(mmcont > vmin))
#monochrome_mmcont[monochrome_mmcont.mask] = 0.0
#monochrome_mmcont.mask[:] = False
rgb_mmcont = monochrome_mmcont[:,:,None]
rgb_mmcont = np.zeros(c18o.shape + (3,))
rgb_mmcont[:,:,blue] = np.nan_to_num(monochrome_mmcont)
hsv_mmcont = rgb_to_hsv(rgb_mmcont)
hue, saturation, value = 0,1,2
hsv_mmcont[:,:,hue] = 150/360.
rgb_mmcont = hsv_to_rgb(hsv_mmcont)
rgb_im[:,:,:alpha] += rgb_mmcont

#rgb_im[rgb_im>1] = 1
#rgb_im /= rgb_im.max()

avm = pyavm.AVM.from_header(outhdr)

im = PIL.Image.fromarray((rgb_im[:,:,:alpha]/rgb_im[:,:,:alpha].max()*255).astype('uint8')[::-1,:])
outname = paths.fpath("rgb_overview_default.png")
im.save(outname)
avm.embed(outname, outname)
im = ImageEnhance.Contrast(im).enhance(1.5)
outname = paths.fpath("rgb_overview_contrast.png")
im.save(outname)
avm.embed(outname, outname)
im = ImageEnhance.Brightness(im).enhance(1.5)
outname = paths.fpath("rgb_overview_brightness.png")
im.save(outname)
avm.embed(outname, outname)
#im = ImageEnhance.Brightness(im).enhance(1.5)
#outname = paths.fpath("rgb_overview_brightness_1.5.png")
#im.save(outname)
#avm.embed(outname, outname)


FF = aplpy.FITSFigure(outname)
FF.show_rgb(outname)
FF.show_regions(paths.rpath('overview_labels.reg'), layer='labels')
FF.save(paths.fpath("rgb_overview_aplpy_withlabels.png"), dpi=150)
