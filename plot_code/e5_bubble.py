import aplpy
from astropy import units as u
from astropy import coordinates
import paths
import pylab as pl
pl.matplotlib.rc_file('pubfiguresrc')

center = coordinates.SkyCoord('19:23:41.849 +14:30:56.70', frame='fk5', unit=(u.hour, u.deg))

ku_fn = '/Users/adam/work/w51/paper_w51_evla/data/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'

fig = pl.figure(1)
fig.clf()
FF = aplpy.FITSFigure(paths.dpath('12m/continuum/selfcal_allspw_selfcal_3_mfs_deeper_r0.0.image.pbcor.fits'),
                      figure=fig)
FF.show_colorscale(cmap='bone_r', vmin=-0.005, vmax=0.025, stretch='log', vmid=-0.006)
FF.recenter(center.ra.deg, center.dec.deg, radius=0.0025)
FF.show_colorbar()
FF.add_scalebar((0.1*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
FF.scalebar.set_color('w')
FF.scalebar.set_linewidth(4)
FF.scalebar.set_label("0.1 pc")
FF.scalebar.set_font_size(20)
FF.add_beam()
FF.beam.set_facecolor('none')
FF.beam.set_linewidth(2)
FF.beam.set_edgecolor('w')
FF.save(paths.fpath("e5_bubble.png"))

fig = pl.figure(1)
fig.clf()
FF = aplpy.FITSFigure(paths.dpath('12m/continuum/selfcal_allspw_selfcal_3_mfs_deeper_r2.0.image.pbcor.fits'),
                      figure=fig)
FF.show_colorscale(cmap='bone_r', vmin=-0.015, vmax=0.095, stretch='log', vmid=-0.016)
FF.recenter(center.ra.deg, center.dec.deg, radius=0.0025)
FF.show_colorbar()
FF.add_scalebar((0.1*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
FF.scalebar.set_color('w')
FF.scalebar.set_linewidth(4)
FF.scalebar.set_label("0.1 pc")
FF.scalebar.set_font_size(20)
FF.show_contour(ku_fn, levels=[0.0015, 0.003, 0.006], colors=['w']*5)
FF.add_beam()
FF.beam.set_facecolor('none')
FF.beam.set_linewidth(2)
FF.beam.set_edgecolor('w')
FF.save(paths.fpath("e5_bubble_robust2.png"))

pl.draw()
pl.show()
