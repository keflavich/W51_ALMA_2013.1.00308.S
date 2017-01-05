from paths import dpath, fpath
from astropy.io import fits
import numpy as np
from astropy import wcs
import aplpy
from astropy.io import fits
import reproject
from vla_cont_cutout import fnku
from matplotlib.colors import Normalize,LogNorm
from matplotlib.colors import rgb_to_hsv,hsv_to_rgb
import PIL
from PIL import ImageEnhance
import pyavm



if 'nir_im' not in locals():

    nir_fn = (dpath('not_ALMA/naco_Kband_W51.fits'))
    nir_hdu = fits.open(nir_fn)
    cm_fn = dpath('not_ALMA/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits')
    cm_hdu = fits.open(cm_fn)
    mm_fn = (dpath('W51_te_continuum_best.fits'))
    mm_hdu = fits.open(mm_fn)

    outhdr = nir_hdu[0].header
    mm_im = reproject.reproject_interp(mm_fn, outhdr, order=1)[0]

    cm_hdu[0].data = cm_hdu[0].data.squeeze()
    cm_hdu[0].header = wcs.WCS(cm_hdu[0].header).celestial.to_header()
    cm_im = reproject.reproject_interp(cm_hdu, outhdr, order=1)[0]

    nir_im = nir_hdu[0].data


# BEGIN IMAGE MAKING HERE
red,green,blue,alpha = 0,1,2,3
hue, saturation, value = 0,1,2

# FIRST LAYER: NIR = white
rgb_im = np.zeros(nir_im.shape + (4,))
rgb_im[:,:,alpha]=1.0
nir_im_rezero = nir_im - np.nanmin(nir_im)
colorized_nir = LogNorm(vmin=np.nanpercentile(nir_im_rezero,30),
                        vmax=np.nanpercentile(nir_im_rezero,99.9995), clip=True)(nir_im_rezero)
colorized_nir.fill_value = 0.0
rgb_nir_im = np.nan_to_num(colorized_nir.filled())[:,:,None] * np.ones(3)
rgb_im[:,:,:alpha] += rgb_nir_im


cm_im_rezero = cm_im - np.nanmin(cm_im)
monochrome_cm_im = LogNorm(vmin=0.0003,#np.nanpercentile(cm_im-np.nanmin(cm_im),30),
                           vmax=np.nanpercentile(cm_im_rezero,99.9995), clip=True)(cm_im_rezero)
monochrome_cm_im.fill_value = 0.0
monochrome_cm_im = monochrome_cm_im.filled()
#rgb_cm_im = np.nan_to_num(monochrome_cm_im)
#rgb_im[:,:,red] += rgb_cm_im

rgb_cm_im = np.zeros(nir_im.shape + (3,))
rgb_cm_im[:,:,blue] = np.nan_to_num(monochrome_cm_im)
#rgb_cm_im[:,:,blue] = monochrome_cm_im.filled()
hsv_cm_im = rgb_to_hsv(rgb_cm_im)
hsv_cm_im[:,:,hue] = 220/360.
rgb_cm_im = hsv_to_rgb(hsv_cm_im)
rgb_im[:,:,:alpha] += rgb_cm_im * 0.8



vmin = 0.0035
monochrome_mm_im = LogNorm(vmin=vmin,
                           vmax=np.nanpercentile(mm_im,99.999995),
                           clip=True)(mm_im*(mm_im > vmin))
monochrome_mm_im.fill_value = 0.0
monochrome_mm_im = monochrome_mm_im.filled()

#vmin = 0.001
#monochrome_mm_im = Normalize(vmin=vmin,
#                             vmax=np.nanpercentile(mm_im,99.99),
#                             clip=True)(mm_im*(mm_im > vmin))

#monochrome_mm_im[monochrome_mm_im.mask] = 0.0
#monochrome_mm_im.mask[:] = False
rgb_mm_im = np.zeros(nir_im.shape + (3,))
rgb_mm_im[:,:,blue] = np.nan_to_num(monochrome_mm_im)
#rgb_mm_im[:,:,blue] = monochrome_mm_im.filled()
hsv_mm_im = rgb_to_hsv(rgb_mm_im)
hsv_mm_im[:,:,hue] = 20/360.
rgb_mm_im = hsv_to_rgb(hsv_mm_im)
rgb_im[:,:,:alpha] += rgb_mm_im

#rgb_im[rgb_im>1] = 1
#rgb_im /= rgb_im.max()

avm = pyavm.AVM.from_header(outhdr)

im = PIL.Image.fromarray((rgb_im[:,:,:alpha]/rgb_im[:,:,:alpha].max()*255).astype('uint8')[::-1,:])
outname = fpath("rgb_irs2_default.png")
im.save(outname)
avm.embed(outname, outname)
im = ImageEnhance.Contrast(im).enhance(1.5)
outname = fpath("rgb_irs2_contrast.png")
im.save(outname)
avm.embed(outname, outname)
im = ImageEnhance.Brightness(im).enhance(1.5)
outname = fpath("rgb_irs2_brightness.png")
im.save(outname)
avm.embed(outname, outname)
#im = ImageEnhance.Brightness(im).enhance(1.5)
#outname = paths.fpath("rgb_irs2_brightness_1.5.png")
#im.save(outname)
#avm.embed(outname, outname)


FF = aplpy.FITSFigure(outname)
FF.show_rgb(outname)
#FF.show_regions(rpath('irs2_labels.reg'), layer='labels')
FF.save(fpath("rgb_irs2_aplpy_withlabels.png"), dpi=150)
