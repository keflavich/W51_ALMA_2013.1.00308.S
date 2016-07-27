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
ch3cn_v = Table.read('ch3cn_v8=1_lines.csv') # a splatalogue query with the overlapping hyperfines removed

def fit_ch3cn_lines(spectra, save_prefix, velo=56*u.km/u.s, ampguess=-0.01):
    all_ch3cn = table.vstack([ch3cn, ch3cn_v])

    all_ch3cn.add_column(table.Column(name='FittedAmplitude', data=np.zeros(len(all_ch3cn))))
    all_ch3cn.add_column(table.Column(name='FittedCenter', data=np.zeros(len(all_ch3cn))))
    all_ch3cn.add_column(table.Column(name='FittedWidth', data=np.zeros(len(all_ch3cn))))
    all_ch3cn.add_column(table.Column(name='FittedAmplitudeError', data=np.zeros(len(all_ch3cn))))
    all_ch3cn.add_column(table.Column(name='FittedCenterError', data=np.zeros(len(all_ch3cn))))
    all_ch3cn.add_column(table.Column(name='FittedWidthError', data=np.zeros(len(all_ch3cn))))

    vkms = velo.to(u.km/u.s).value
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
                           guesses=[mid, ampguess, vkms, 2],)
                line['FittedAmplitude'] = sp.specfit.parinfo['AMPLITUDE0'].value
                line['FittedCenter'] = sp.specfit.parinfo['SHIFT0'].value
                line['FittedWidth'] = sp.specfit.parinfo['WIDTH0'].value
                line['FittedAmplitudeError'] = sp.specfit.parinfo['AMPLITUDE0'].error
                line['FittedCenterError'] = sp.specfit.parinfo['SHIFT0'].error
                line['FittedWidthError'] = sp.specfit.parinfo['WIDTH0'].error
                sp.xarr.convert_to_unit(u.GHz)
            pb.update()

    pl.xlim(vkms-14, vkms+14)
    #pl.ylim(0, offset)
    pl.draw()
    pl.show()
    pl.savefig(save_prefix+"_spectra_overlay.png")

    pl.figure(2).clf()
    pl.plot(all_ch3cn['E_U (K)'], all_ch3cn['FittedWidth'], 'o')
    pl.xlabel("E$_U$ (K)")
    pl.ylabel("$\sigma$ (km/s)")
    pl.ylim(0,3.5)
    pl.savefig(save_prefix+"_sigma_vs_eupper.png")


    pl.figure(3).clf()
    pl.plot(all_ch3cn['E_U (K)'], all_ch3cn['FittedCenter'], 'o')
    pl.xlabel("E$_U$ (K)")
    pl.ylabel("$v_{lsr}$ (km/s)")
    pl.ylim(vkms-3, vkms+3)
    pl.savefig(save_prefix+"_vcen_vs_eupper.png")


spectra_center = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_center_W51e2_spw*fits')))
spectra_left = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_left_W51e2_spw*fits')))
spectra_right = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_right_W51e2_spw*fits')))
spectra_se_emission = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_se_emission_W51e2_spw*fits')))
spectra_ALMAmm24 = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/ALMAmm24_W51n_spw*.fits')))
spectra_d2 = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/d2_W51n_spw*.fits')))
spectra_e2e = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_W51e2_spw*.fits')))
spectra_e2nw = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2nw_W51e2_spw*.fits')))
spectra_e2w = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2w_W51e2_spw*.fits')))
spectra_e8 = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e8_W51e2_spw*.fits')))
spectra_north = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/north_W51n_spw*.fits')))

fit_ch3cn_lines(spectra_center, paths.fpath('longbaseline/ch3cn_fits/e2e_center'), velo=56*u.km/u.s)
fit_ch3cn_lines(spectra_right, paths.fpath('longbaseline/ch3cn_fits/e2e_right'), velo=56*u.km/u.s)
fit_ch3cn_lines(spectra_left, paths.fpath('longbaseline/ch3cn_fits/e2e_left'), velo=56*u.km/u.s)
fit_ch3cn_lines(spectra_se_emission, paths.fpath('longbaseline/ch3cn_fits/e2e_se_emission'), velo=59*u.km/u.s, ampguess=0.01)
fit_ch3cn_lines(spectra_d2, paths.fpath('longbaseline/ch3cn_fits/d2'), velo=56*u.km/u.s)
fit_ch3cn_lines(spectra_e8, paths.fpath('longbaseline/ch3cn_fits/e8'), velo=60*u.km/u.s)
fit_ch3cn_lines(spectra_north, paths.fpath('longbaseline/ch3cn_fits/north'), velo=61*u.km/u.s)
fit_ch3cn_lines(spectra_ALMAmm24, paths.fpath('longbaseline/ch3cn_fits/ALMAmm24'), velo=61*u.km/u.s)
