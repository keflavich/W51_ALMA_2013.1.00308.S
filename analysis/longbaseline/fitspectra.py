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
import pyregion
import velo_guesses

line_table = table.Table.read(paths.apath('full_line_table.csv'))
line_table.sort('Species')

regions = (pyregion.open(paths.rpath("cores_longbaseline_spectralextractionregions_pix.reg"))+
           pyregion.open(paths.rpath("cores_longbaseline_spectralextractionregions_pix_north.reg")))

for region in regions:
    name = region.attr[1]['text']

    spectral_files = glob.glob(paths.dpath('longbaseline/spectra/{0}_W51*_spw*.fits'.format(name)))
    spectra = pyspeckit.Spectra(spectral_files)
    stats = spectra.stats()
    err = stats['std'] # overly conservative guess
    med = stats['median']

    line_table.add_column(table.Column(name='{0}FittedAmplitude'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedCenter'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedWidth'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedAmplitudeError'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedCenterError'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedWidthError'.format(name), data=np.zeros(len(line_table))))

    vcen = velo_guesses.guesses[name]

    plotnum = 1

    pl.figure(1).clf()

    for ii,line_row in enumerate(line_table):

        frq = line_row['Freq-GHz'] * u.GHz

        linename = line_row['Species'] + line_row['Resolved QNs']

        # check if line in range
        for sp in spectra:
            if sp.xarr.in_range(frq):
                
                sp.xarr.convert_to_unit(u.km/u.s, refX=frq)

                # crop
                sp_sl = sp.slice((vcen-20)*u.km/u.s, (vcen+20)*u.km/u.s, unit=u.km/u.s)

                # convert original back
                sp.xarr.convert_to_unit(u.GHz)

                # make sure we're in km/s
                sp_sl.xarr.convert_to_unit(u.km/u.s, refX=frq)
                # guess at the errors
                sp_sl.error[:] = err
                
                # perform fit
                sp_sl.specfit(fittype='vheightgaussian',
                              guesses=[med, 1, vcen, 2],
                              fixed=[True,False,False,False], # use a fixed baseline
                              limited=[(True,True)]*4,
                              limits=[(0,1), (-5,5), (vcen-5, vcen+5), (0,10)],
                             )

                sp_sl.plotter(axis=pl.subplot(6,6,plotnum))
                plotnum += 1
                sp_sl.specfit.plot_fit(annotate=False)
                sp_sl.plotter.axis.set_ylabel('')
                sp_sl.plotter.axis.set_xlabel('')
                sp_sl.plotter.axis.annotate(linename, (0.5, 0.85), xycoords='axes fraction')

                # write fit to table
                line_row['{0}FittedAmplitude'.format(name)] = sp_sl.specfit.parinfo['AMPLITUDE0'].value
                line_row['{0}FittedCenter'.format(name)] = sp_sl.specfit.parinfo['SHIFT0'].value
                line_row['{0}FittedWidth'.format(name)] = sp_sl.specfit.parinfo['WIDTH0'].value
                line_row['{0}FittedAmplitudeError'.format(name)] = sp_sl.specfit.parinfo['AMPLITUDE0'].error
                line_row['{0}FittedCenterError'.format(name)] = sp_sl.specfit.parinfo['SHIFT0'].error
                line_row['{0}FittedWidthError'.format(name)] = sp_sl.specfit.parinfo['WIDTH0'].error

    pl.savefig(paths.fpath('longbaseline/spectral_fits/{0}_linefits.png'.format(name)))

    line_table.write(paths.tpath('longbaseline/spectral_lines_and_fits.csv'),
                     overwrite=True)
