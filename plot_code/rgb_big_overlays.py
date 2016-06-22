import pylab as pl
import paths
import aplpy
from astropy import coordinates
from astropy import units as u
from astropy.table import Table

rgb_cube_fits = 'c18o_h2co_ku_rgb.fits'
rgb_cube_png = rgb_cube_fits[:-5]+"_logred.png"

pl.rcParams['font.size'] = 18
fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(rgb_cube_png, figure=fig1)
F.show_rgb(rgb_cube_png)
#F.recenter(290.93315, 14.509584, radius=0.0075)
F.add_scalebar((0.5*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
F.scalebar.set_label('0.5 pc')
F.scalebar.set_color('w')
#F.set_tick_xspacing(0.0005)
F.add_label(0.05, 0.95, "14.5 GHz Continuum", relative=True, color='r', horizontalalignment='left')
F.add_label(0.05, 0.91, "H$_2$CO", relative=True, color='g', horizontalalignment='left')
F.add_label(0.05, 0.87, "C$^{18}$O", relative=True, color='b', horizontalalignment='left')
F.save(paths.fpath("W51e2_ku_h2co_c18o_aplpy.png"))
F.save(paths.fpath("W51e2_ku_h2co_c18o_aplpy.pdf"))

cmcontsrc = Table.read(paths.vpath('tables/EVLA_VLA_PointSourcePhotometry.ipac'),
                       format='ascii.ipac')
cmok = (cmcontsrc['Frequency']==5.9) & (cmcontsrc['Epoch']=='3')
cmcoords = coordinates.SkyCoord(cmcontsrc['gracen'][cmok],
                                cmcontsrc['gdeccen'][cmok], frame='fk5')


core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
cores = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                             frame='fk5')

F.show_markers(cmcoords.ra, cmcoords.dec, edgecolor='w', marker='*', alpha=0.75,
               zorder=1, facecolor='w', layer='hiiregions')
#F.hide_layer('hiiregions_text')
F.save(paths.fpath("W51e2_ku_h2co_c18o_aplpy_hiiregions.png"))
F.save(paths.fpath("W51e2_ku_h2co_c18o_aplpy_hiiregions.pdf"))

F.show_markers(cores.ra, cores.dec, edgecolor='y', marker='.', alpha=0.9,
               zorder=1, facecolor='y', layer='cores')
F.save(paths.fpath("W51e2_ku_h2co_c18o_aplpy_hiiregions_cores.png"))
F.save(paths.fpath("W51e2_ku_h2co_c18o_aplpy_hiiregions_cores.pdf"))

F.hide_layer('hiiregions')
F.save(paths.fpath("W51e2_ku_h2co_c18o_aplpy_cores.png"))
F.save(paths.fpath("W51e2_ku_h2co_c18o_aplpy_cores.pdf"))
