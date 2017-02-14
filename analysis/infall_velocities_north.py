import numpy as np
from astropy import units as u
import pyspeckit
import radio_beam
import paths

F = False
T = True


sp = pyspeckit.Spectrum(paths.spath('north_radial_bin_0.00to0.38_spw1.fits'))
# CH3OH 8-7
sp.xarr.refX=220.07849*u.GHz
sp.xarr.convert_to_unit(u.km/u.s)
sp.data = sp.data_quantity.to(u.K,
                              u.brightness_temperature(radio_beam.Beam.from_fits_header(sp.header),
                                                       sp.xarr.to(u.GHz))).value
# hack =(
sp.data = np.ma.masked_where(np.isnan(sp.data), sp.data)
sp.unit = u.K

sp.plotter(xmin=20,xmax=110, reset_ylimits=True)
# Hack - xarr is backwards
sp.baseline.selectregion(exclude=[-100000,37,43,96,103,100000][::-1])
sp.baseline.highlight_fitregion()
sp.baseline(order=1, subtract=False, plot=True, reset_selection=False,
            selectregion=False, fit_plotted_area=False,
            highlight_fitregion=True)

sp.specfit.selectregion(xmin=48, xmax=63.5)
sp.specfit.register_fitter('hill5', pyspeckit.spectrum.models.hill5infall.hill5_fitter, 5)
sp.specfit(fittype='hill5', guesses=[0.3, 56.5, 3.63, 2.048, 100],
           fixed=[F,F,F,F,F], reset_selection=False)



sp = pyspeckit.Spectrum(paths.spath('north_radial_bin_0.00to0.38_spw0.fits'))
# CH3OH 4-3
sp.xarr.refX= 218.44005 * u.GHz
sp.xarr.convert_to_unit(u.km/u.s)
sp.data = sp.data_quantity.to(u.K,
                              u.brightness_temperature(radio_beam.Beam.from_fits_header(sp.header),
                                                       sp.xarr.to(u.GHz))).value
# hack =(
sp.data = np.ma.masked_where(np.isnan(sp.data), sp.data)
sp.unit = u.K

sp.plotter(xmin=20,xmax=110, reset_ylimits=True)
sp.baseline.basespec[:]=78
sp.baseline.baselinepars=[0,78]
#sp.baseline.plot_baseline()

sp.specfit.selectregion(xmin=48, xmax=63.5)
sp.specfit.register_fitter('hill5', pyspeckit.spectrum.models.hill5infall.hill5_fitter, 5)
sp.specfit(fittype='hill5', guesses=[0.3, 56.5, 3.63, 2.048, 100],
           fixed=[F,F,F,F,F], reset_selection=False)


sp = pyspeckit.Spectrum(paths.spath('north_radial_bin_0.00to0.38_spw2.fits'))
# CH3OH 10-9
sp.xarr.refX= 231.28115 * u.GHz
sp.xarr.convert_to_unit(u.km/u.s)
sp.data = sp.data_quantity.to(u.K,
                              u.brightness_temperature(radio_beam.Beam.from_fits_header(sp.header),
                                                       sp.xarr.to(u.GHz))).value
# hack =(
sp.data = np.ma.masked_where(np.isnan(sp.data), sp.data)
sp.unit = u.K

sp.plotter(xmin=20,xmax=110, reset_ylimits=True)
sp.baseline.basespec[:]=82
sp.baseline.baselinepars=[0,82]
#sp.baseline.plot_baseline()

sp.specfit.selectregion(xmin=48, xmax=63.5)
sp.specfit.register_fitter('hill5', pyspeckit.spectrum.models.hill5infall.hill5_fitter, 5)
sp.specfit(fittype='hill5', guesses=[0.3, 56.5, 3.63, 2.048, 100],
           fixed=[F,F,F,F,F], reset_selection=False)



sp = pyspeckit.Spectrum(paths.spath('north_radial_bin_0.00to0.38_spw3.fits'))
# CH3OH 18-17
sp.xarr.refX= 233.79580 * u.GHz
sp.xarr.convert_to_unit(u.km/u.s)
sp.data = sp.data_quantity.to(u.K,
                              u.brightness_temperature(radio_beam.Beam.from_fits_header(sp.header),
                                                       sp.xarr.to(u.GHz))).value
# hack =(
sp.data = np.ma.masked_where(np.isnan(sp.data), sp.data)
sp.unit = u.K

sp.plotter(xmin=20,xmax=110, reset_ylimits=True)
sp.baseline.basespec[:]=91
sp.baseline.baselinepars=[0,91]
#sp.baseline.plot_baseline()

sp.specfit.selectregion(xmin=48, xmax=63.5)
sp.specfit.register_fitter('hill5', pyspeckit.spectrum.models.hill5infall.hill5_fitter, 5)
sp.specfit(fittype='hill5', guesses=[0.3, 56.5, 3.63, 2.048, 100],
           fixed=[F,F,F,F,F], reset_selection=False)

