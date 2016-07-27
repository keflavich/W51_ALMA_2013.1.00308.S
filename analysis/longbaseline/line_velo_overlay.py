import glob
import numpy as np
import paths
import pyspeckit
from astropy import constants
from astropy import units as u
from astropy.table import Table
from astropy import table
from astropy.utils.console import ProgressBar
import pylab as pl

ch3cn = Table.read('ch3cn_v=0_lines.csv')
ch3cn_v = Table.read('ch3cn_v8=1_lines.csv')
all_ch3cn = table.vstack([ch3cn, ch3cn_v])

all_ch3cn.add_column(table.Column(name='FittedAmplitude', data=np.zeros(len(all_ch3cn))))
all_ch3cn.add_column(table.Column(name='FittedCenter', data=np.zeros(len(all_ch3cn))))
all_ch3cn.add_column(table.Column(name='FittedWidth', data=np.zeros(len(all_ch3cn))))

spectra_center = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_center_W51e2_spw*fits')))
spectra_left = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_left_W51e2_spw*fits')))
spectra_right = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_right_W51e2_spw*fits')))

spectra=spectra_right

velo = 56*u.km/u.s

pl.figure(1).clf()
ax = pl.gca()

pb = ProgressBar(len(spectra) * len(all_ch3cn))

ii = 0
for sp in spectra:
    sp.xarr.convert_to_unit(u.GHz)
    mid = np.median(sp.data)
    for line in all_ch3cn:
        frq = line['Freq-GHz']*u.GHz
        if sp.xarr.in_range(frq*(1-velo/constants.c)):
            offset = ii*0.000 + mid
            ii += 1
            sp.xarr.convert_to_unit(u.km/u.s, refX=frq)
            sp.plotter(axis=ax, clear=False, offset=offset)
            sp.specfit(fittype='vheightgaussian',
                       guesses=[mid, -0.01, 56, 2],)
            line['FittedAmplitude'] = sp.specfit.parinfo['AMPLITUDE0'].value
            line['FittedCenter'] = sp.specfit.parinfo['SHIFT0'].value
            line['FittedWidth'] = sp.specfit.parinfo['WIDTH0'].value
            sp.xarr.convert_to_unit(u.GHz)
        pb.update()

pl.xlim(42,70)
pl.ylim(0, offset)
pl.draw()
pl.show()

pl.figure(2).clf()
pl.plot(all_ch3cn['E_U (K)'], all_ch3cn['FittedWidth'], 'o')
pl.xlabel("E$_U$ (K)")
pl.ylabel("$\sigma$ (km/s)")
pl.ylim(1,3.5)


pl.figure(3).clf()
pl.plot(all_ch3cn['E_U (K)'], all_ch3cn['FittedCenter'], 'o')
pl.xlabel("E$_U$ (K)")
pl.ylabel("$v_{lsr}$ (km/s)")
pl.ylim(53, 58)
