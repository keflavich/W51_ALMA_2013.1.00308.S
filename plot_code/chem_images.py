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

filelist = glob.glob(paths.dpath('12m/cutouts/*e2e8*fits'))

labeldict = {
               'H2CO303_202':'H$_2$CO $3_{0,3}-2_{0,2}$',
               'H2CO321_220':'H$_2$CO $3_{2,1}-2_{2,0}$',
               'H2CO322_221':'H$_2$CO $3_{2,2}-2_{2,1}$',
               'CH3OH422-312':'CH$_3$OH $4_{2,2}-3_{1,2}$',
               'HC3N24-23':'HC$_3$N 24-23',
               'OCS18-17':'OCS 18-17',
               'OCS19-18':'OCS 19-18',
               'SO65-54':'SO $6_5-5_4$',
               'HNCO10110-919':'HNCO $10_{1,10}-9_{1,9}$',
               'HNCO1028-927': 'HNCO $10_{2,8}-9_{2,7}$',
               'CH3OH423-514':  'CH$_3$OH $4_{2,3}-5_{1,4}$',
               'CH3OH5m42-6m43':'CH$_3$OH $5_{-4,2}-6_{-4,3}$',
               'CH3OH808-716':'CH$_3$OH $8_{0,8}-7_{1,6}$',
               '13CS5-4':'$^{13}$CS 5-4',
               'CH3OCH3_13013-12112':'CH$_3$OCH$_3$ $13_{0,13}-12_{1,12}$',
               'NH2CHO11210-1029':'NH$_2$CHO $11_{2,10}-10_{2,9}$',
               'NH2CHO1156-1055': 'NH$_2$CHO $11_{5,6}-10_{5,5}$',
               'HC3Nv7=124-23':'HC$_3$Nv$_7$=1 24-23a',
               'H30alpha':'H30$\\alpha$',
               'C18O2-1':'C$^{18}$O 2-1',
               'H2CCO11-10':'H$_2$CCO 11-10',
               'HCOOH431-524':'HCOOH $4_{3,1}-5_{2,4}$',
               'CH3OCHO17314-16313E':'CH$_3$OCHO $17_{3,14}-16_{3,13}$E',
               'CH3CH2CN24321-23320':'CH$_3$CH$_2$CN $24_{3,21}-23_{3,20}$',
               'HC3Nv7=1_24-23':'HC$_3$Nv$_7$=1 24-23b',
               'Acetone21120-20219AE':'Acetone $21_{1,20}-20_{2,19}$AE',
               'Acetone21120-20119EE':'Acetone $21_{1,20}-20_{1,19}$EE',
               'CH3CH2CN24222-23221':'CH$_3$CH$_2$CN $24_{2,22}-23_{2,21}$',
               'H213CO312-211':'H$_2$$^{13}$CO $3_{1,2}-2_{1,1}$',
               'H2CN303-202_F3_2':'H$_2$CN $3_{0,3}-2{0,2} F$_{3/2-3/2}$',
               'H2CN322-221_F5_2':'H$_2$CN $3_{2,2}-2{2,1} F$_{5/2-3/2}$',
               'CH3OCHO17413-16412A':'CH$_3$OCHO $17_{4,13}-16_{4,12}$A',
               'CH3CH2OH550-541':'CH$_3$CH$_2$OH $5_{5,0}-5_{4,1}$',
               'CH3OCH313013-12112AA':'CH$_3$OCH$_3$ $13_{0,13}-12_{1,12}$AA',
               'NH2CHO1156-1055':'NH$_2$CHO $11_{5,6}-10_{5,5}$',
}

if not all(pl.figure(1, figsize=(8,8)).get_size_inches() == 8):
    pl.close(1)

fig = pl.figure(1, figsize=(8,8))
fig.clf()
gs = gridspec.GridSpec(5,5)
gs.update(wspace=0.0, hspace=0.0)

for ii,fn in enumerate(ProgressBar(filelist)):

    cube = SpectralCube.read(fn)[:,367:467,114:214]
    bm = cube.beams[0]
    jtok = bm.jtok(cube.wcs.wcs.restfrq*u.Hz)
    cube = cube.to(u.K, bm.jtok_equiv(cube.wcs.wcs.restfrq*u.Hz))

    slab = cube.spectral_slab(51*u.km/u.s, 60*u.km/u.s)
    cube.beam_threshold = 1
    #contguess = cube.spectral_slab(0*u.km/u.s, 40*u.km/u.s).percentile(50, axis=0)
    #contguess = cube.spectral_slab(70*u.km/u.s, 100*u.km/u.s).percentile(50, axis=0)
    mask = (cube.spectral_axis<40*u.km/u.s) | (cube.spectral_axis > 75*u.km/u.s)
    contguess = cube.with_mask(mask[:,None,None]).percentile(30, axis=0)
    slabsub = (slab-contguess)
    slabsub.beam_threshold = 0.02
    m0 = slabsub.moment0()

    label = labeldict[linere.search(fn).groups()[0]]

    ax = fig.add_subplot(gs[ii])

    im = ax.imshow(m0.value, vmin=-1.25*jtok.value, vmax=5.0*jtok.value, cmap=pl.cm.bone_r,
                   interpolation='nearest')
    ax.text(3, 87, label, fontsize=10)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_aspect('equal')

bottom,top,left,right = gs.get_grid_positions(fig)
cbar_ax = fig.add_axes([np.max(right)+0.01, np.min(bottom), 0.05, np.max(top)-np.min(bottom)])
cb = pl.colorbar(mappable=im, cax=cbar_ax)
cb.set_label("Flux Density (K km s$^{-1}$)", fontsize=12)
cb.ax.tick_params(labelsize=12)

pl.savefig(paths.fpath("chemical_m0_slabs_e2e8.png"), bbox_inches='tight', dpi=150)

pl.draw()
pl.show()
