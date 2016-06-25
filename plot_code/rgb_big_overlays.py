import pylab as pl
import paths
import aplpy
from astropy import coordinates
from astropy import units as u
from astropy.table import Table

for (rgb_cube_fits, rgb_cube_png, star_color, core_color, rlabel, glabel,
     blabel) in (('full_h2co_rgb.fits', 'full_h2co_rgb_auto.png', 'y', 'b',
                  'H$_2$CO $3_{0,3}-2_{0,2}$', 'H$_2$CO $3_{2,1}-2_{2,0}$',
                  'H$_2$CO $3_{2,2}-2_{2,1}$',),
                 ('c18o_h2co_ku_rgb.fits', 'c18o_h2co_ku_rgb_logred.png', 'w',
                  'y', "14.5 GHz Continuum", "H$_2$CO", "C$^{18}$O",
                 )):
    name = rgb_cube_fits[:-9]

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
    F.add_label(0.05, 0.95, rlabel, relative=True, color='r', horizontalalignment='left')
    F.add_label(0.05, 0.91, glabel, relative=True, color='g', horizontalalignment='left')
    F.add_label(0.05, 0.87, blabel, relative=True, color='b', horizontalalignment='left')
    F.save(paths.fpath("W51e2_{0}_aplpy.png".format(name)))
    F.save(paths.fpath("W51e2_{0}_aplpy.pdf".format(name)))

    cmcontsrc = Table.read(paths.vpath('tables/EVLA_VLA_PointSourcePhotometry.ipac'),
                           format='ascii.ipac')
    cmok = (cmcontsrc['Frequency']==5.9) & (cmcontsrc['Epoch']=='3')
    cmcoords = coordinates.SkyCoord(cmcontsrc['gracen'][cmok],
                                    cmcontsrc['gdeccen'][cmok], frame='fk5')


    core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
    cores = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                                 frame='fk5')

    F.show_markers(cmcoords.ra, cmcoords.dec, edgecolor=star_color, marker='*', alpha=0.75,
                   zorder=1, facecolor=star_color, layer='hiiregions')
    #F.hide_layer('hiiregions_text')
    F.save(paths.fpath("W51e2_{0}_aplpy_hiiregions.png".format(name)))
    F.save(paths.fpath("W51e2_{0}_aplpy_hiiregions.pdf".format(name)))

    F.show_markers(cores.ra, cores.dec, edgecolor=core_color, marker='.', alpha=0.9,
                   zorder=1, facecolor=core_color, layer='cores')
    F.save(paths.fpath("W51e2_{0}_aplpy_hiiregions_cores.png".format(name)))
    F.save(paths.fpath("W51e2_{0}_aplpy_hiiregions_cores.pdf".format(name)))

    F.hide_layer('hiiregions')
    F.save(paths.fpath("W51e2_{0}_aplpy_cores.png".format(name)))
    F.save(paths.fpath("W51e2_{0}_aplpy_cores.pdf".format(name)))
