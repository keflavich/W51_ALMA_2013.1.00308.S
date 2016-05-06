import aplpy
from astropy import units as u
from astropy import coordinates
import paths
import pylab as pl

center = coordinates.SkyCoord('19:23:41.849 +14:30:56.70', frame='fk5', unit=(u.hour, u.deg))
fig = pl.figure(1)
fig.clf()
FF = aplpy.FITSFigure(paths.dpath('12m/continuum/selfcal_allspw_selfcal_3_mfs_deeper_r0.0.image.pbcor.fits'),
                      figure=fig)
FF.show_colorscale(cmap='bone_r', vmin=-0.005, vmax=0.025, stretch='log', vmid=-0.006)
FF.recenter(center.ra.deg, center.dec.deg, radius=0.0025)
FF.show_colorbar()
FF.save(paths.fpath("e5_bubble.png"))

pl.draw()
pl.show()
