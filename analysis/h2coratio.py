import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
import paths
from constants import distance

p303 = paths.dpath('w51_H2CO_303_202_contsub.image.pbcor.fits')
p321 = paths.dpath('w51_H2CO_321_220_contsub.image.pbcor.fits')
cube303 = SpectralCube.read(p303).with_spectral_unit(u.km/u.s,
                                                     velocity_convention='radio')
min_slices = cube303.subcube_slices_from_mask(cube303.mask)
cube321 = SpectralCube.read(p321).with_spectral_unit(u.km/u.s,
                                                     velocity_convention='radio')
# tight cropping
cube303 = cube303[min_slices]
cube321 = cube321[min_slices]

std = cube303[:10].std(axis=0)
mask = cube303 > 3*std

# sad hacks: these are the same to very high but not infinite precision
cube321._wcs = cube303._wcs
cube321.mask._wcs = cube321.wcs

int303 = cube303.with_mask(mask).moment0()
int321 = cube321.with_mask(mask).moment0()
int303.hdu.writeto(paths.dpath('moments/w51_H2CO_303_202_contsub.mom0.fits'),
                   clobber=True)
int321.hdu.writeto(paths.dpath('moments/w51_H2CO_321_220_contsub.mom0.fits'),
                   clobber=True)
int303.quicklook()
int321.quicklook()
r = int321/int303

hdu = int303.hdu
hdu.data = r.value
hdu.writeto(paths.dpath('moments/ratio_321to303_mom0.fits'), clobber=True)

import pylab as pl
import matplotlib
cm = matplotlib.cm.RdYlBu_r
cm.set_bad('#888888')

import aplpy
fig1 = pl.figure(1)
pl.clf()
FF = aplpy.FITSFigure(hdu, figure=fig1)
FF.show_colorscale(cmap=cm, vmin=0, vmax=1)
FF.show_colorbar()
FF.save(paths.fpath('H2CO_321_to_303_ratiomap.png'))
FF.show_contour(paths.dpath('evla/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'),
                colors=['k'], #levels=[0.001,0.002,0.004,0.008,0.016,0.032,0.064],
                levels=np.logspace(-3,-1),
                linewidth=0.5,
                alpha=0.2, layer='black_contours')
FF.save(paths.fpath('H2CO_321_to_303_ratiomap_withcontours.png'))
FF.recenter(290.91644,14.518939,radius=0.15/60.)
FF.save(paths.fpath('H2CO_321_to_303_ratiomap_withcontours_IRS2.png'))
FF.recenter(290.93268,14.508363,radius=0.15/60.)
FF.save(paths.fpath('H2CO_321_to_303_ratiomap_withcontours_e1e2.png'))
FF.show_contour(paths.dpath("w51_te_continuum_best.fits"),
                levels=[0.02, 0.04, 0.08, 0.16],
                colors=['g'],
                layer='almate_cont_ours')
FF.save(paths.fpath('H2CO_321_to_303_ratiomap_withcontours_e1e2_almacont.png'))
FF.recenter(290.91644,14.518939,radius=0.15/60.)
FF.save(paths.fpath('H2CO_321_to_303_ratiomap_withcontours_IRS2_almacont.png'))

from h2co_modeling import lte_model
ratio = lte_model.T_321/lte_model.T_303
vals = r.value[np.isfinite(r.value)]
tems = np.interp(vals, ratio[np.isfinite(ratio)], np.array(lte_model.tem)[np.isfinite(ratio)])
newr = r.value.copy()
newr[np.isfinite(r.value)] = tems

hdu2 = int303.hdu
hdu2.data = newr
hdu2.writeto(paths.dpath('moments/temperature_LTE_321to303_mom0.fits'), clobber=True)
fig2 = pl.figure(2)
pl.clf()
FF2 = aplpy.FITSFigure(hdu2, figure=fig2)
FF2.show_colorscale(cmap=cm, vmin=10, vmax=200, stretch='log', vmid=-50)
FF2.show_colorbar()
FF2.colorbar.set_axis_label_text("Temperature [K]")
FF2.add_scalebar((1*u.pc/distance).to(u.deg, u.dimensionless_angles()).value,)
FF2.scalebar.set_label('1 pc')
FF2.save(paths.fpath('H2CO_321_to_303_LTEtemperaturemap.png'))
FF2.show_contour(paths.dpath('evla/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'),
                 colors=['k'], levels=[0.001], layer='black_contours')
FF2.save(paths.fpath('H2CO_321_to_303_LTEtemperaturemap_withCMcontours.png'))
FF2.hide_layer('black_contours')
FF2.show_contour(paths.dpath('evla/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'),
                 colors=['w'], levels=[0.001], layer='white_contours')
FF2.save(paths.fpath('H2CO_321_to_303_LTEtemperaturemap_withwhiteCMcontours.png'))
FF2.hide_layer('white_contours')
FF2.show_regions(paths.rpath('cores.reg'), layer='cores')
FF2.save(paths.fpath('H2CO_321_to_303_LTEtemperaturemap_withcores.png'))
