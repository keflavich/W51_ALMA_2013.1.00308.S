import pylab as pl
import paths
import aplpy
from astropy import coordinates
from astropy import units as u
from astropy.table import Table
import os
from astropy import log
from vla_cont_cutout import fitsKu_cutout
from outflow_meta import e1

for suffix in ('auto','99.99'):
    # note to self: rgb_cube_fits doesn't have to exist, it's just a prefix
    # for "name"
    for (rgb_cube_fits, rgb_cube_png, star_color, core_color, rlabel, glabel,
         blabel) in (('full_h2co_rgb.fits', 'full_h2co_rgb_auto.png', 'y', 'b',
                      'H$_2$CO $3_{0,3}-2_{0,2}$', 'H$_2$CO $3_{2,1}-2_{2,0}$',
                      'H$_2$CO $3_{2,2}-2_{2,1}$',),
                     ('c18o_h2co_ku_rgb.fits', 'c18o_h2co_ku_rgb_logred.png', 'w',
                      'r', "14.5 GHz Continuum", "H$_2$CO", "C$^{18}$O",
                     ),
                     ('hc3n_ch3oh_ocs_rgb.fits', 'hc3n_ch3oh_ocs_rgb_auto.png', 'y', 'b',
                      'HC$_3$N', 'CH$_3$OH', 'OCS'),
                     ('h2co_hc3n_ch3oh_rgb.fits', 'h2co_hc3n_ch3oh_rgb_auto.png', 'w', 'b',
                      'H$_2$CO', 'HC$_3$N', 'CH$_3$OH'),
                     ('h2co_so_hc3n_rgb.fits', 'h2co_so_hc3n_rgb_auto.png', 'w', 'b',
                      'H$_2$CO', 'SO', 'HC$_3$N'),
                     ('ku_hc3n_ch3oh_rgb.fits', 'ku_hc3n_ch3oh_rgb_auto.png', 'w', 'b',
                      'Ku', 'HC$_3$N', 'CH$_3$OH'),
                     ('ku_so_c18o_rgb.fits', 'ku_so_c18o_rgb_auto.png', 'w', 'b',
                      'Ku', 'SO', 'C$^{18}$O',),
                    ):
        name = rgb_cube_fits[:-9]
        # stupid hack, should just rewrite from scratch but it's easier to edit inplace...
        if suffix != 'auto':
            name += suffix
            rgb_cube_png = rgb_cube_png.replace('auto', suffix)
            if not os.path.exists(rgb_cube_png):
                log.warn("Skipping {0}: does not exist".format(rgb_cube_png))
                continue

        pl.rcParams['font.size'] = 18
        fig1 = pl.figure(1)
        fig1.clf()
        F = aplpy.FITSFigure(rgb_cube_png, figure=fig1)
        F.show_rgb(rgb_cube_png)
        #F.recenter(290.93315, 14.509584, radius=0.0075)
        F.add_scalebar((0.05*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
        #F.scalebar.set_length((0.05*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
        F.scalebar.set_label('0.05 pc')
        F.scalebar.set_color('w')

        # zoom on "e1 cluster"
        F.recenter(e1.ra.deg, e1.dec.deg, radius=11./3600.)
        F.set_tick_xspacing(6./3600)

        #F.set_tick_xspacing(0.0005)
        #F.add_label(0.05, 0.95, rlabel, relative=True, color='r', horizontalalignment='left')
        #F.add_label(0.05, 0.91, glabel, relative=True, color='g', horizontalalignment='left')
        #F.add_label(0.05, 0.87, blabel, relative=True, color='b', horizontalalignment='left')
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


        F.scalebar.set_label('0.05 pc')
        F.scalebar.set_color('w')

        F.show_markers(cmcoords.ra, cmcoords.dec, edgecolor=star_color, marker='*', alpha=0.75,
                       zorder=1, facecolor=star_color, layer='hiiregions')
        F.save(paths.fpath("rgb_zooms/W51e1_{0}_aplpy_hiiregions.png".format(name)), dpi=300)
        F.save(paths.fpath("rgb_zooms/W51e1_{0}_aplpy_hiiregions.pdf".format(name)), dpi=300)

        F.show_markers(cores.ra, cores.dec, edgecolor=core_color, marker='.', alpha=0.9,
                       zorder=1, facecolor=core_color, layer='cores')
        F.save(paths.fpath("rgb_zooms/W51e1_{0}_aplpy_hiiregions_cores.png".format(name)), dpi=300)
        F.save(paths.fpath("rgb_zooms/W51e1_{0}_aplpy_hiiregions_cores.pdf".format(name)), dpi=300)

        F.hide_layer('hiiregions')
        F.save(paths.fpath("rgb_zooms/W51e1_{0}_aplpy_cores.png".format(name)), dpi=300)
        F.save(paths.fpath("rgb_zooms/W51e1_{0}_aplpy_cores.pdf".format(name)), dpi=300)

        F.hide_layer('cores')
        F.show_contour(fitsKu_cutout, colors=['w'], levels=[0.001, 0.005, 0.010, 0.05, 0.1])
        F.save(paths.fpath("rgb_zooms/W51e1_{0}_aplpy_cmcontours.png".format(name)), dpi=300)
        F.save(paths.fpath("rgb_zooms/W51e1_{0}_aplpy_cmcontours.pdf".format(name)), dpi=300)
