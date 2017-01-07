import paths
import pylab as pl
import os
from spectral_cube import SpectralCube
import numpy as np
from astropy import units as u
import matplotlib.gridspec as gridspec
from astropy.utils.console import ProgressBar
import radio_beam
from spectral_cube.lower_dimensional_structures import Projection
from astropy.io import fits
from astropy import wcs
from astropy.stats import mad_std
from line_to_image_list import labeldict

import re
import glob

# Sometimes debugging is necessary to prevent abort traps
from astropy import log
log.setLevel('DEBUG')

# PN 5-4 looks terrible and doesn't have 12m data; it is screwy.
if 'PN5-4' in labeldict:
    labeldict.pop('PN5-4')

def chem_plot(linere, yslice=slice(367,467), xslice=slice(114,214), vrange=[51,60]*u.km/u.s,
              sourcename='e2', filelist=glob.glob(paths.dpath('12m/cutouts/*e2e8*fits')),
              # 5,8 -> 12.8,8
              suffix="", plotgrid=(6,8), figsize=(12.0,8.6),
              vmax_m0=5.0,
              vmax_max=150,
              maxbeam=0.5*u.arcsec,
              contourlevels=None,
             ):
    nplots = np.product(plotgrid)

    for ii in (1,2):
        if not all(pl.figure(ii, figsize=figsize).get_size_inches() == figsize):
            pl.close(ii)

    fig1 = pl.figure(1, figsize=figsize)
    fig1.clf()
    gs1 = gridspec.GridSpec(*plotgrid)
    gs1.update(wspace=0.0, hspace=0.0)

    fig2 = pl.figure(2, figsize=figsize)
    fig2.clf()
    gs2 = gridspec.GridSpec(*plotgrid)
    gs2.update(wspace=0.0, hspace=0.0)

    fig3 = pl.figure(3, figsize=figsize)
    fig3.clf()
    gs3 = gridspec.GridSpec(*plotgrid)
    gs3.update(wspace=0.0, hspace=0.0)

    fig4 = pl.figure(4, figsize=figsize)
    fig4.clf()
    gs4 = gridspec.GridSpec(*plotgrid)
    gs4.update(wspace=0.0, hspace=0.0)

    fig5 = pl.figure(5, figsize=figsize)
    fig5.clf()
    gs5 = gridspec.GridSpec(*plotgrid)
    gs5.update(wspace=0.0, hspace=0.0)

    fig6 = pl.figure(6, figsize=figsize)
    fig6.clf()
    gs6 = gridspec.GridSpec(*plotgrid)
    gs6.update(wspace=0.0, hspace=0.0)


    figcounter = 0

    for ii,fn in enumerate(ProgressBar(filelist)):

        linename = linere.search(fn).groups()[0]
        if linename not in labeldict:
            print("Skipping {0} because it's not in the label dict".format(linename))
            continue
        label = labeldict[linename]

        # cache the results for use in other work, later use, ...
        m0fitsfn = paths.dpath("chemslices/chemical_m0_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
        m1fitsfn = paths.dpath("chemslices/chemical_m1_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
        m2fitsfn = paths.dpath("chemslices/chemical_m2_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
        maxfitsfn = paths.dpath("chemslices/chemical_max_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
        maxsubfitsfn = paths.dpath("chemslices/chemical_max_sub_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
        madstdfitsfn = paths.dpath("chemslices/chemical_madstd_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
        if not (os.path.exists(m0fitsfn) and
                os.path.exists(m2fitsfn) and
                os.path.exists(maxfitsfn) and
                os.path.exists(madstdfitsfn)):
            print("Extracting max/m0/m1/m2 for {0}".format(fn))
            cube = SpectralCube.read(fn)[:,yslice,xslice]
            goodbeams = np.array([bm.major < maxbeam for bm in cube.beams], dtype='bool')
            cube = cube.with_mask(goodbeams[:,None,None])
            cube = cube.minimal_subcube()

            if cube.shape[0] == 0:
                print("Skipping {0} because it was masked out".format(fn))
                continue

            bm = cube.beams[0]
            restfreq = cube.wcs.wcs.restfrq
            cube = cube.to(u.K, bm.jtok_equiv(restfreq*u.Hz))

            slab = cube.spectral_slab(*vrange)
            cube.beam_threshold = 1
            #contguess = cube.spectral_slab(0*u.km/u.s, 40*u.km/u.s).percentile(50, axis=0)
            #contguess = cube.spectral_slab(70*u.km/u.s, 100*u.km/u.s).percentile(50, axis=0)
            mask = (cube.spectral_axis<40*u.km/u.s) | (cube.spectral_axis > 75*u.km/u.s)
            try:
                contguess = cube.with_mask(mask[:,None,None]).percentile(30, axis=0)
            except ValueError as ex:
                print("skipping {0}".format(fn))
                print(ex)
                continue
            slabsub = (slab-contguess)
            slab.beam_threshold = 0.25
            slabsub.beam_threshold = 0.25
            m0 = slabsub.moment0()
            m1 = slabsub.moment1()
            m2 = slabsub.moment2()
            max_sub = slabsub.max(axis=0)
            max = slab.max(axis=0)
            madstd = cube.with_mask(mask[:,None,None]).apply_function(mad_std,
                                                                      axis=0,
                                                                      projection=True,
                                                                      progressbar=True,
                                                                      unit=cube.unit,)

            m0.write(m0fitsfn, overwrite=True)
            m1.write(m1fitsfn, overwrite=True)
            max.write(maxfitsfn, overwrite=True)
            max_sub.write(maxsubfitsfn, overwrite=True)
            madstd.write(madstdfitsfn, overwrite=True)
        else:
            m0fh = fits.open(m0fitsfn)
            m1fh = fits.open(m1fitsfn)
            m2fh = fits.open(m2fitsfn)
            maxfh = fits.open(maxfitsfn)
            maxsubfh = fits.open(maxsubfitsfn)
            madstdfh = fits.open(madstdfitsfn)

            m0 = Projection(value=m0fh[0].data, header=m0fh[0].header,
                            wcs=wcs.WCS(m0fh[0].header),
                            unit=u.Unit(m0fh[0].header['BUNIT']),)
            m1 = Projection(value=m1fh[0].data, header=m1fh[0].header,
                            wcs=wcs.WCS(m1fh[0].header),
                            unit=u.Unit(m1fh[0].header['BUNIT']),)
            m2 = Projection(value=m2fh[0].data, header=m2fh[0].header,
                            wcs=wcs.WCS(m2fh[0].header),
                            unit=u.Unit(m2fh[0].header['BUNIT']),)
            max = Projection(value=maxfh[0].data, header=maxfh[0].header,
                             wcs=wcs.WCS(maxfh[0].header),
                             unit=u.Unit(maxfh[0].header['BUNIT']),)
            max_sub = Projection(value=maxsubfh[0].data,
                                 header=maxsubfh[0].header,
                                 wcs=wcs.WCS(maxsubfh[0].header),
                                 unit=u.Unit(maxsubfh[0].header['BUNIT']),)
            madstd = Projection(value=madstdfh[0].data,
                                header=madstdfh[0].header,
                                wcs=wcs.WCS(madstdfh[0].header),
                                unit=u.Unit(madstdfh[0].header['BUNIT']),)

            bm = radio_beam.Beam.from_fits_header(m0fh[0].header)
            restfreq = m0fh[0].header['RESTFRQ']

        jtok = bm.jtok(restfreq*u.Hz)

        if figcounter>=nplots:
            print("Skipping {0}".format(fn))
            break

        ax1 = fig1.add_subplot(gs1[figcounter])

        im1 = ax1.imshow(m0.value, vmin=-1.25*jtok.value, vmax=vmax_m0*jtok.value,
                         cmap=pl.cm.bone_r, interpolation='nearest')
        ax1.text(3, 0.87*m0.shape[0], label, fontsize=9)
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_aspect('equal')

        ax2 = fig2.add_subplot(gs2[figcounter])

        im2 = ax2.imshow(m1.value, vmin=vrange[0].value, vmax=vrange[1].value,
                         cmap=pl.cm.jet, interpolation='nearest')
        ax2.text(3, 0.87*m0.shape[0], label, fontsize=9, color='w')
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_aspect('equal')

        ax3 = fig3.add_subplot(gs3[figcounter])

        im3 = ax3.imshow(max_sub.value, vmin=-10, vmax=vmax_max,
                         cmap=pl.cm.bone_r, interpolation='nearest')
        # add a contour to show the regions that are "saturated" above T_max
        if contourlevels is None:
            contourlevels = [vmax_max, 300, 400, 500]
        qcs = ax3.contour(max_sub.value, levels=contourlevels, colors=['r','g','b','y'])
        #print("levels: {0} = {1}".format(qcs.levels, contourlevels))
        ax3.text(3, 0.87*m0.shape[0], label, fontsize=9, color='r')
        ax3.set_xticklabels([])
        ax3.set_yticklabels([])
        ax3.set_aspect('equal')

        ax5 = fig5.add_subplot(gs5[figcounter])

        im5 = ax5.imshow(max.value, vmin=-10, vmax=vmax_max,
                         cmap=pl.cm.bone_r, interpolation='nearest')
        # add a contour to show the regions that are "saturated" above T_max
        qcs = ax5.contour(max.value, levels=contourlevels, colors=['r','g','b','y'])
        #print("levels: {0} = {1}".format(qcs.levels, contourlevels))
        ax5.text(3, 0.87*m0.shape[0], label, fontsize=9, color='r')
        ax5.set_xticklabels([])
        ax5.set_yticklabels([])
        ax5.set_aspect('equal')

        ax4 = fig4.add_subplot(gs4[figcounter])

        im4 = ax4.imshow(madstd.value,
                         cmap=pl.cm.bone_r, interpolation='nearest')
        ax4.text(3, 0.87*m0.shape[0], label, fontsize=9, color='r')
        ax4.set_xticklabels([])
        ax4.set_yticklabels([])
        ax4.set_aspect('equal')

        ax6 = fig6.add_subplot(gs6[figcounter])

        im6 = ax6.imshow(m2.value, vmin=0, vmax=10,
                         cmap=pl.cm.jet, interpolation='nearest')
        ax6.text(3, 0.87*m0.shape[0], label, fontsize=9, color='w')
        ax6.set_xticklabels([])
        ax6.set_yticklabels([])
        ax6.set_aspect('equal')

        figcounter += 1

    cbs = {}
    for ii,fig, im, gs in ((1,fig1,im1,gs1), (2,fig2,im2,gs2),
                           (3,fig3,im3,gs3), (4,fig4,im4,gs4),
                           (5,fig5,im5,gs5), (6,fig6,im6,gs6),
                          ):
        bottom,top,left,right = gs.get_grid_positions(fig)
        cbar_ax = fig.add_axes([np.max(right)+0.01, np.min(bottom),
                                0.05, np.max(top)-np.min(bottom)])
        cbs[ii] = pl.colorbar(mappable=im, cax=cbar_ax)
        cbs[ii].ax.tick_params(labelsize=12)

    cbs[1].set_label("Flux Density (K km s$^{-1}$)", fontsize=12)
    cbs[2].set_label("Velocity (km s$^{-1}$)", fontsize=12)
    cbs[3].set_label("Peak Brightness (K)", fontsize=12)
    cbs[5].set_label("Peak Brightness (K)", fontsize=12)
    cbs[4].set_label("MAD StdDev (K)", fontsize=12)
    cbs[6].set_label("Velocity Dispersion (km s$^{-1}$)", fontsize=12)

    pl.draw()
    pl.show()

    fig1.savefig(paths.fpath("chemical_m0_slabs_{0}{1}.png".format(sourcename,
                                                                   suffix)),
                 bbox_inches='tight', dpi=150)

    fig2.savefig(paths.fpath("chemical_m1_slabs_{0}{1}.png".format(sourcename,
                                                                   suffix)),
                 bbox_inches='tight', dpi=150)

    fig3.savefig(paths.fpath("chemical_max_contsub_slabs_{0}{1}.png"
                             .format(sourcename, suffix)),
                 bbox_inches='tight', dpi=150)

    fig4.savefig(paths.fpath("chemical_madstd_slabs_{0}{1}.png".format(sourcename,
                                                                       suffix)),
                 bbox_inches='tight', dpi=150)

    fig5.savefig(paths.fpath("chemical_max_slabs_{0}{1}.png"
                             .format(sourcename, suffix)),
                 bbox_inches='tight', dpi=150)

    fig6.savefig(paths.fpath("chemical_m2_slabs_{0}{1}.png".format(sourcename,
                                                                   suffix)),
                 bbox_inches='tight', dpi=150)


def test_layout(plotgrid=(6,8), figsize=(12.0,8.6), fignum=1):

    # test layout
    pl.close(fignum)
    fig = pl.figure(fignum, figsize=figsize)
    fig.clf()
    gs = gridspec.GridSpec(*plotgrid)
    gs.update(wspace=0.0, hspace=0.0)

    for ii in range(np.product(plotgrid)):
        ax = fig.add_subplot(gs[ii])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_aspect('equal')
        im = pl.imshow([[1,1,],[1,2]])

    bottom,top,left,right = gs.get_grid_positions(fig)
    cbar_ax = fig.add_axes([np.max(right)+0.01, np.min(bottom),
                            0.05, np.max(top)-np.min(bottom)])
    cb = pl.colorbar(mappable=im, cax=cbar_ax)
    cb.ax.tick_params(labelsize=12)



if __name__ == "__main__":



    linere = re.compile("W51_b6_12M.(.*).image.pbcor")

    chem_plot(linere, yslice=slice(357,477), xslice=slice(104,224),
              vrange=[51,60]*u.km/u.s, sourcename='e2',
              vmax_max=200,
              contourlevels=[150,200,250,300],
              filelist=glob.glob(paths.dpath('12m/cutouts/W51_b6_12M.*e2e8*fits')))

    chem_plot(linere, yslice=slice(227,347), xslice=slice(119,239),
              vrange=[52,63]*u.km/u.s, sourcename='e8',
              vmax_max=200,
              contourlevels=[150,200,250,300],
              filelist=glob.glob(paths.dpath('12m/cutouts/W51_b6_12M.*e2e8*fits')))

    chem_plot(linere, yslice=slice(31,231), xslice=slice(152,350),
              vrange=[54,64]*u.km/u.s, sourcename='north',
              vmax_max=200,
              contourlevels=[150,200,250,300],
              filelist=glob.glob(paths.dpath('12m/cutouts/W51_b6_12M.*north*fits')))

    chem_plot(linere, yslice=slice(50,150), xslice=slice(80,180),
              vrange=[58,67]*u.km/u.s, sourcename='ALMAmm14',
              vmax_max=100,
              contourlevels=[125,150,175,200],
              filelist=glob.glob(paths.dpath('12m/cutouts/W51_b6_12M.*ALMAmm14*fits')))

    ### Below here are not used in the paper

    chem_plot(linere, yslice=slice(83,133), xslice=slice(243,293),
              vrange=[58,67]*u.km/u.s, sourcename='d2',
              vmax_max=200,
              contourlevels=[150,200,250,300],
              filelist=glob.glob(paths.dpath('12m/cutouts/W51_b6_12M.*north*fits')))


    linere = re.compile("W51_b6_7M_12M.(.*).image.pbcor")

    chem_plot(linere, yslice=slice(357,477), xslice=slice(104,224),
              vrange=[51,60]*u.km/u.s, sourcename='e2',
              filelist=glob.glob(paths.dpath('merge/cutouts/W51_b6_7M_12M.*e2e8*fits')),
              suffix="_merge", vmax_max=200)
    # not implemented chem_plot(linere, yslice=slice(357,477), xslice=slice(104,224),
    # not implemented           vrange=[51,60]*u.km/u.s, sourcename='e2',
    # not implemented           filelist=glob.glob(paths.dpath('merge/cutouts/W51_b6_7M_12M.*e2e8*fits')),
    # not implemented           continuum=paths.dpath('W51_te_continuum_best.fits'),
    # not implemented           continuum_levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
    # not implemented                             0.282361, 0.4, ],
    # not implemented           suffix="_merge_cont")
    chem_plot(re.compile("W51_b6_7M_12M_natural.(.*).image.pbcor"), yslice=slice(168,249), xslice=slice(42,118),
              vrange=[51,60]*u.km/u.s, sourcename='e2',
              filelist=glob.glob(paths.dpath('merge/cutouts/W51_b6_7M_12M_natural.*e2e8*fits')),
              suffix="_merge_natural",
              maxbeam=1.5*u.arcsec,
              vmax_m0=7.5)

    chem_plot(linere, yslice=slice(227,347), xslice=slice(119,239),
              vrange=[52,63]*u.km/u.s, sourcename='e8',
              filelist=glob.glob(paths.dpath('merge/cutouts/W51_b6_7M_12M.*e2e8*fits')), suffix="_merge")

    chem_plot(linere, yslice=slice(31,231), xslice=slice(152,350),
              vrange=[54,64]*u.km/u.s, sourcename='north',
              filelist=glob.glob(paths.dpath('merge/cutouts/W51_b6_7M_12M.*north*fits')), suffix="_merge")

    chem_plot(linere, yslice=slice(50,150), xslice=slice(80,180),
              vrange=[58,67]*u.km/u.s, sourcename='ALMAmm14',
              filelist=glob.glob(paths.dpath('merge/cutouts/W51_b6_7M_12M.*ALMAmm14*fits')), suffix="_merge")

    chem_plot(linere, yslice=slice(83,133), xslice=slice(243,293),
              vrange=[58,67]*u.km/u.s, sourcename='d2',
              vmax_max=200,
              contourlevels=[150,200,250,300],
              filelist=glob.glob(paths.dpath('merge/cutouts/W51_b6_7M_12M.*north*fits')))


