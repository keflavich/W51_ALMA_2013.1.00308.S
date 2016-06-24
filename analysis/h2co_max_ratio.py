from astropy.io import fits
from astropy import units as u
from h2co_modeling import lte_model
from constants import distance
import numpy as np
import paths
import reproject
import pylab as pl
import aplpy
import matplotlib
cm = matplotlib.cm.RdYlBu_r
cm.set_bad('#888888')


fn303merge = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO303_202.image.pbcor_medsub_max.fits')
fn321merge = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO321_220.image.pbcor_medsub_max.fits')
fn322merge = paths.dpath('merge/moments/W51_b6_7M_12M.H2CO322_221.image.pbcor_medsub_max.fits')
fn303_12m = paths.dpath('12m/moments/W51_b6_12M.H2CO303_202.image.pbcor_medsub_max.fits')
fn321_12m = paths.dpath('12m/moments/W51_b6_12M.H2CO321_220.image.pbcor_medsub_max.fits')
fn322_12m = paths.dpath('12m/moments/W51_b6_12M.H2CO322_221.image.pbcor_medsub_max.fits')

fh321merge = fits.open(fn321merge)
fh322merge = fits.open(fn322merge)
fh303merge = fits.open(fn303merge)

reproj321merge = reproject.reproject_interp(fh321merge[0], fh303merge[0].header)
ratio321merge_fh = fits.PrimaryHDU(data=reproj321merge[0]/fh303merge[0].data, header=fh303merge[0].header)
ratio321merge_fh.data[fh303merge[0].data < 0.02] = np.nan
ratio321merge_fh.data[ratio321merge_fh.data < 0.01] = np.nan
ratio321merge_fh.writeto(paths.dpath('merge/moments/H2CO_321_to_303_max_ratio.fits'), clobber=True)

reproj322merge = reproject.reproject_interp(fh322merge[0], fh303merge[0].header)
ratio322merge_fh = fits.PrimaryHDU(data=reproj322merge[0]/fh303merge[0].data, header=fh303merge[0].header)
ratio322merge_fh.data[fh303merge[0].data < 0.02] = np.nan
ratio322merge_fh.data[ratio322merge_fh.data < 0.015] = np.nan
ratio322merge_fh.writeto(paths.dpath('merge/moments/H2CO_322_to_303_max_ratio.fits'), clobber=True)

fh321_12m = fits.open(fn321_12m)
fh322_12m = fits.open(fn322_12m)
fh303_12m = fits.open(fn303_12m)

reproj321_12m = reproject.reproject_interp(fh321_12m[0], fh303_12m[0].header)
ratio321_12m_fh = fits.PrimaryHDU(data=reproj321_12m[0]/fh303_12m[0].data, header=fh303_12m[0].header)
ratio321_12m_fh.data[fh303_12m[0].data < 0.01] = np.nan
ratio321_12m_fh.data[ratio321_12m_fh.data < 0.005] = np.nan
ratio321_12m_fh.writeto(paths.dpath('12m/moments/H2CO_321_to_303_max_ratio.fits'), clobber=True)

reproj322_12m = reproject.reproject_interp(fh322_12m[0], fh303_12m[0].header)
ratio322_12m_fh = fits.PrimaryHDU(data=reproj322_12m[0]/fh303_12m[0].data, header=fh303_12m[0].header)
ratio322_12m_fh.data[fh303_12m[0].data < 0.01] = np.nan
ratio322_12m_fh.data[ratio322_12m_fh.data < 0.015] = np.nan
ratio322_12m_fh.writeto(paths.dpath('12m/moments/H2CO_322_to_303_max_ratio.fits'), clobber=True)

modelratio321303 = lte_model.T_321/lte_model.T_303
modelratio322303 = lte_model.T_322/lte_model.T_303

def ratio_to_temperature(data, modelratio=modelratio321303):
    vals = data[np.isfinite(data)]
    tems = np.interp(vals, modelratio[np.isfinite(modelratio)],
                     np.array(lte_model.tem)[np.isfinite(modelratio)])
    newr = data.copy()
    newr[np.isfinite(data)] = tems
    return newr

t321303_12m = ratio_to_temperature(ratio321_12m_fh.data, modelratio=modelratio321303)
t322303_12m = ratio_to_temperature(ratio322_12m_fh.data, modelratio=modelratio322303)
t321303merge = ratio_to_temperature(ratio321merge_fh.data, modelratio=modelratio321303)
t322303merge = ratio_to_temperature(ratio322merge_fh.data, modelratio=modelratio322303)

hdutem321_12m = fits.PrimaryHDU(data=t321303_12m, header=ratio321_12m_fh.header)
hdutem322_12m = fits.PrimaryHDU(data=t322303_12m, header=ratio322_12m_fh.header)
hdutem321merge = fits.PrimaryHDU(data=t321303merge, header=ratio321merge_fh.header)
hdutem322merge = fits.PrimaryHDU(data=t322303merge, header=ratio322merge_fh.header)

hdutem321_12m.writeto(paths.dpath('12m/moments/temperature_LTE_321to303_mom0.fits'), clobber=True)
hdutem322_12m.writeto(paths.dpath('12m/moments/temperature_LTE_322to303_mom0.fits'), clobber=True)
hdutem321merge.writeto(paths.dpath('merge/moments/temperature_LTE_321to303_mom0.fits'), clobber=True)
hdutem322merge.writeto(paths.dpath('merge/moments/temperature_LTE_322to303_mom0.fits'), clobber=True)

for hdu, rhdu, label in ((hdutem321_12m, ratio321_12m_fh, '321_to_303_max_12m'),
                         (hdutem322_12m, ratio322_12m_fh, '322_to_303_max_12m'),
                         (hdutem321merge, ratio321merge_fh, '321_to_303_max_merge'),
                         (hdutem322merge, ratio322merge_fh, '322_to_303_max_merge'),
                        ):
        
    fig1 = pl.figure(1)
    fig1.clf()
    FF1 = aplpy.FITSFigure(rhdu, figure=fig1)
    FF1.show_colorscale(cmap=cm, vmin=0, vmax=1, stretch='log', vmid=-50)
    FF1.show_colorbar()
    FF1.colorbar.set_axis_label_text("Ratio $E_U=68/E_U=23$")
    FF1.add_scalebar((1*u.pc/distance).to(u.deg, u.dimensionless_angles()).value,)
    FF1.scalebar.set_label('1 pc')
    FF1.save(paths.fpath('H2CO_{0}_RatioMap.png'.format(label)))


    fig2 = pl.figure(2)
    pl.clf()
    FF2 = aplpy.FITSFigure(hdu, figure=fig2)
    FF2.show_colorscale(cmap=cm, vmin=10, vmax=120, stretch='log', vmid=-50)
    FF2.show_colorbar()
    FF2.colorbar.set_axis_label_text("Temperature [K]")
    FF2.add_scalebar((1*u.pc/distance).to(u.deg, u.dimensionless_angles()).value,)
    FF2.scalebar.set_label('1 pc')
    FF2.save(paths.fpath('H2CO_{0}_LTEtemperaturemap.png'.format(label)))
    FF2.show_contour(paths.dpath('evla/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'),
                     colors=['k'], levels=[0.001], layer='black_contours')
    FF2.save(paths.fpath('H2CO_{0}_LTEtemperaturemap_withCMcontours.png'.format(label)))
    FF2.hide_layer('black_contours')
    FF2.show_contour(paths.dpath('evla/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'),
                     colors=['w'], levels=[0.001], layer='white_contours')
    FF2.save(paths.fpath('H2CO_{0}_LTEtemperaturemap_withwhiteCMcontours.png'.format(label)))
    FF2.hide_layer('white_contours')
    FF2.show_regions(paths.rpath('cores.reg'), layer='cores')
    FF2.save(paths.fpath('H2CO_{0}_LTEtemperaturemap_withcores.png'.format(label)))
