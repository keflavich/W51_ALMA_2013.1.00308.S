import paths
import pylab as pl
import os
from spectral_cube import SpectralCube
import numpy as np
from astropy import units as u
import matplotlib.gridspec as gridspec
from astropy.utils.console import ProgressBar

import re
import glob

linere = re.compile("W51_b6_12M.(.*).image.pbcor")

labeldict = {
               'H2CO303_202':'H$_2$CO $3_{0,3}-2_{0,2}$',
               'H2CO321_220':'H$_2$CO $3_{2,1}-2_{2,0}$',
               'H2CO322_221':'H$_2$CO $3_{2,2}-2_{2,1}$',
               'HC3N24-23':'HC$_3$N 24-23',
               'OCS18-17':'OCS 18-17',
               'OCS19-18':'OCS 19-18',
               'SO65-54':'SO $6_5-5_4$',
               'HNCO10110-919':'HNCO $10_{1,10}-9_{1,9}$',
               'HNCO1028-927': 'HNCO $10_{2,8}-9_{2,7}$',
               'CH3OH423-514':  'CH$_3$OH $4_{2,3}-5_{1,4}$',
               'CH3OH5m42-6m43':'CH$_3$OH $5_{-4,2}-6_{-4,3}$',
               'CH3OH808-716':'CH$_3$OH $8_{0,8}-7_{1,6}$',
               'CH3OH1029-936':'CH$_3$OH $10_{2,9}-9_{3,6}$',
               'CH3OH1028-937':'CH$_3$OH $10_{2,8}-9_{3,7}$',
               'CH3OH422-312':'CH$_3$OH $4_{2,2}-3_{1,2}$',
               'CH3OHvt=01028-937++': 'CH$_3$OH $10_{2,8}-9_{3,7}$',
               '13CS5-4':'$^{13}$CS 5-4',
               'NH2CHO11210-1029':'NH$_2$CHO $11_{2,10}-10_{2,9}$',
               'NH2CHO1156-1055': 'NH$_2$CHO $11_{5,6}-10_{5,5}$',
               'HC3Nv7=124-23':'HC$_3$Nv$_7$=1 24-23a',
               'H30alpha':'H30$\\alpha$',
               'C18O2-1':'C$^{18}$O 2-1',
               'H2CCO11-10':'H$_2$CCO 11-10',
               'HCOOH431-524':'HCOOH $4_{3,1}-5_{2,4}$',
               'CH3CH2CN24321-23320':'CH$_3$CH$_2$CN $24_{3,21}-23_{3,20}$',
               'HC3Nv7=1_24-23':'HC$_3$Nv$_7$=1 24-23b',
               'Acetone21120-20219AE':'Acetone $21_{1,20}-20_{2,19}$AE',
               'Acetone21120-20119EE':'Acetone $21_{1,20}-20_{1,19}$EE',
               'CH3CH2CN24222-23221':'CH$_3$CH$_2$CN $24_{2,22}-23_{2,21}$',
               'H213CO312-211':'H$_2$$^{13}$CO $3_{1,2}-2_{1,1}$',
               'H2CN303-202_F3_2':'H$_2$CN $3_{0,3}-2{0,2}$ F$_{3/2-3/2}$',
               'H2CN322-221_F5_2':'H$_2$CN $3_{2,2}-2{2,1}$ F$_{5/2-3/2}$',
               'CH3OCHO17413-16412A':'CH$_3$OCHO $17_{4,13}-16_{4,12}$A',
               'CH3OCHO17314-16313E':'CH$_3$OCHO $17_{3,14}-16_{3,13}$E',
               'CH3OCHO17314-16313A':'CH$_3$OCHO $17_{3,14}-16_{3,13}$A',
               'CH3CH2OH550-541':'CH$_3$CH$_2$OH $5_{5,0}-5_{4,1}$',
               'CH3OCH323321-23222AA':'CH$_3$OCH$_3$ $23_{3,21}-23_{2,22}$AA',
               'CH3OCH313013-12112AA':'CH$_3$OCH$_3$ $13_{0,13}-12_{1,12}$AA',
               'CH3OCH323321-23222EE':'CH$_3$OCH$_3$ $23_{3,21}-23_{2,22}$EE',
               'HNCO10010-909':'HNCO $10_{0,10}-9_{0,9}$',
               'O13CS18-17':'O$^{13}$CS 18-17',
               'N2D+_3-2': 'N$_2$D+ 3-2',
}

def chem_plot(yslice=slice(367,467), xslice=slice(114,214), vrange=[51,60]*u.km/u.s,
              sourcename='e2', filelist=glob.glob(paths.dpath('12m/cutouts/*e2e8*fits'))
             ):
    for ii in (1,2):
        if not all(pl.figure(ii, figsize=(12.8,8)).get_size_inches() == (12.8,8)):
            pl.close(ii)

    fig1 = pl.figure(1, figsize=(12.8,8))
    fig1.clf()
    gs1 = gridspec.GridSpec(5,8)
    gs1.update(wspace=0.0, hspace=0.0)

    fig2 = pl.figure(2, figsize=(12.8,8))
    fig2.clf()
    gs2 = gridspec.GridSpec(5,8)
    gs2.update(wspace=0.0, hspace=0.0)

    for ii,fn in enumerate(ProgressBar(filelist)):

        cube = SpectralCube.read(fn)[:,yslice,xslice]
        bm = cube.beams[0]
        jtok = bm.jtok(cube.wcs.wcs.restfrq*u.Hz)
        cube = cube.to(u.K, bm.jtok_equiv(cube.wcs.wcs.restfrq*u.Hz))

        slab = cube.spectral_slab(*vrange)
        cube.beam_threshold = 1
        #contguess = cube.spectral_slab(0*u.km/u.s, 40*u.km/u.s).percentile(50, axis=0)
        #contguess = cube.spectral_slab(70*u.km/u.s, 100*u.km/u.s).percentile(50, axis=0)
        mask = (cube.spectral_axis<40*u.km/u.s) | (cube.spectral_axis > 75*u.km/u.s)
        contguess = cube.with_mask(mask[:,None,None]).percentile(30, axis=0)
        slabsub = (slab-contguess)
        slabsub.beam_threshold = 0.15
        m0 = slabsub.moment0()

        label = labeldict[linere.search(fn).groups()[0]]

        ax1 = fig1.add_subplot(gs1[ii])

        im1 = ax1.imshow(m0.value, vmin=-1.25*jtok.value, vmax=5.0*jtok.value, cmap=pl.cm.bone_r,
                         interpolation='nearest')
        ax1.text(3, 0.87*m0.shape[0], label, fontsize=10)
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_aspect('equal')

        m1 = slabsub.moment1()
        ax2 = fig2.add_subplot(gs2[ii])

        im2 = ax2.imshow(m1.value, vmin=vrange[0].value, vmax=vrange[1].value,
                         cmap=pl.cm.viridis, interpolation='nearest')
        ax2.text(3, 0.87*m0.shape[0], label, fontsize=10, color='w')
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_aspect('equal')

    cbs = {}
    for ii,fig, im, gs in ((1,fig1,im1,gs1), (2,fig2,im2,gs2)):
        bottom,top,left,right = gs.get_grid_positions(fig)
        cbar_ax = fig.add_axes([np.max(right)+0.01, np.min(bottom), 0.05, np.max(top)-np.min(bottom)])
        cbs[ii] = pl.colorbar(mappable=im, cax=cbar_ax)
        cbs[ii].ax.tick_params(labelsize=12)

    cbs[1].set_label("Flux Density (K km s$^{-1}$)", fontsize=12)
    cbs[2].set_label("Velocity (km s$^{-1}$)", fontsize=12)

    fig1.savefig(paths.fpath("chemical_m0_slabs_{0}.png".format(sourcename)), bbox_inches='tight', dpi=150)

    fig2.savefig(paths.fpath("chemical_m1_slabs_{0}.png".format(sourcename)), bbox_inches='tight', dpi=150)

    pl.draw()
    pl.show()


chem_plot(yslice=slice(367,467), xslice=slice(114,214),
          vrange=[51,60]*u.km/u.s, sourcename='e2',
          filelist=glob.glob(paths.dpath('12m/cutouts/*e2e8*fits')))

chem_plot(yslice=slice(227,347), xslice=slice(119,239),
          vrange=[52,63]*u.km/u.s, sourcename='e8',
          filelist=glob.glob(paths.dpath('12m/cutouts/*e2e8*fits')))

chem_plot(yslice=slice(31,231), xslice=slice(152,350),
          vrange=[54,64]*u.km/u.s, sourcename='north',
          filelist=glob.glob(paths.dpath('12m/cutouts/*north*fits')))

chem_plot(yslice=slice(50,150), xslice=slice(80,180),
          vrange=[58,67]*u.km/u.s, sourcename='ALMAmm14',
          filelist=glob.glob(paths.dpath('12m/cutouts/*ALMAmm14*fits')))
