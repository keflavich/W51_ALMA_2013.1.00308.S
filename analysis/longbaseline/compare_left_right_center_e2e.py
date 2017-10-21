import re
import pyspeckit
import glob
import paths
import aplpy
import pyregion
import pylab as pl
from astropy import units as u

pl.figure(1).clf()
pl.figure(2).clf()

# sanity check: are the regions in the right place?
F = aplpy.FITSFigure(paths.dpath('longbaseline/W51e2cax.cont.image.pbcor.fits'),
                     figure=pl.figure(1), slices=[0,1])
F.show_grayscale(vmax=0.015)
region_list = pyregion.open(paths.rpath("cores_longbaseline_spectralextractionregions.reg"))
F.show_regions(region_list)
F.recenter(290.9332, 14.509589, 0.5/3600.)

spectra_se_emission = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_se_emission_W51e2_spw*fits')))
spectra_center = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_center_W51e2_spw*fits')))
spectra_left = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_left_W51e2_spw*fits')))
spectra_right = pyspeckit.Spectra(glob.glob(paths.dpath('longbaseline/spectra/e2e_right_W51e2_spw*fits')))

spectra_center.xarr.convert_to_unit(u.GHz)
spectra_right.xarr.convert_to_unit(u.GHz)
spectra_left.xarr.convert_to_unit(u.GHz)

spectra_center.plotter(figure=pl.figure(2))
spectra_left.plotter(axis=spectra_center.plotter.axis, clear=False, color='r')
spectra_right.plotter(axis=spectra_center.plotter.axis, clear=False, color='g')

for sp in spectra_center:
    pl.xlim(sp.xarr.min().to(u.GHz).value, sp.xarr.max().to(u.GHz).value)
    pl.ylim(0.006, 0.027)
    spwnum = re.compile("spw([0-9])").search(sp.fileprefix).groups()[0]
    pl.savefig(paths.fpath('longbaseline/kinematics/e2_center_left_right_spw{0}.png'.format(spwnum)))
