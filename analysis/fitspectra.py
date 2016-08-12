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

import warnings
from numpy.ma.core import MaskedArrayFutureWarning
warnings.filterwarnings('ignore', category=MaskedArrayFutureWarning)

line_table = table.Table.read(paths.apath('full_line_table.csv'))
line_table.sort('Species')

regions = (pyregion.open(paths.rpath("cores.reg")))

for region in regions:
    name = region.attr[1]['text']

    spectral_files = glob.glob(paths.spath('{0}_spw*.fits'.format(name)))
    spectra = pyspeckit.Spectra(spectral_files)
    stats = spectra.stats()
    err = stats['std'] # overly conservative guess

    for sp in spectra:
        sp.data -= np.nanpercentile(sp.data, 25)
    #med = stats['median']
    #if med < 0:
    #    med = 0

    line_table.add_column(table.Column(name='{0}FittedAmplitude'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedCenter'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedWidth'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedAmplitudeError'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedCenterError'.format(name), data=np.zeros(len(line_table))))
    line_table.add_column(table.Column(name='{0}FittedWidthError'.format(name), data=np.zeros(len(line_table))))

    vcen = velo_guesses.guesses[name]

    plotnum = 1

    pl.figure(1).clf()

    for ii,line_row in enumerate(ProgressBar(line_table)):

        frq = line_row['Freq-GHz'] * u.GHz

        linename = line_row['Species'] + line_row['Resolved QNs']

        # check if line in range
        for sp in spectra:
            sp.xarr.convert_to_unit(u.GHz)
            if sp.xarr.in_range(frq):

                
                # not needed!
                #sp.xarr.convert_to_unit(u.km/u.s, refX=frq)
                assert sp.xarr.unit == u.GHz
                sp.xarr.refX = frq

                # crop
                sp_sl = sp.slice((vcen-20)*u.km/u.s, (vcen+20)*u.km/u.s, unit=u.km/u.s)

                if not np.any(np.isfinite(sp_sl.data)):
                    continue
                print(name, linename, plotnum)

                # convert original back
                # not needed! sp.xarr.convert_to_unit(u.GHz)

                # make sure we're in km/s
                sp_sl.xarr.convert_to_unit(u.km/u.s, refX=frq)
                # guess at the errors
                sp_sl.error[:] = err
                
                # perform fit
                sp_sl.specfit(fittype='gaussian',
                              guesses=[0.1, vcen, 2],
                              fixed=[False,False,False], # use a fixed baseline
                              limited=[(True,True)]*3,
                              limits=[(-5,5), (vcen-2.5, vcen+2.5), (0,5)],
                             )

                if np.abs(sp_sl.specfit.parinfo[0].value) > sp_sl.specfit.parinfo[0].error*2:
                    if plotnum <= 49:
                        sp_sl.plotter(axis=pl.subplot(7,7,plotnum))
                        plotnum += 1
                        sp_sl.specfit.plot_fit(annotate=False)
                        sp_sl.plotter.axis.set_ylabel('')
                        sp_sl.plotter.axis.set_xlabel('')
                        sp_sl.plotter.axis.annotate(linename, (0.05, 0.85), xycoords='axes fraction')
                    else:
                        print("Skipping line because too many plots")

                    # write fit to table
                    line_row['{0}FittedAmplitude'.format(name)] = sp_sl.specfit.parinfo['AMPLITUDE0'].value
                    line_row['{0}FittedCenter'.format(name)] = sp_sl.specfit.parinfo['SHIFT0'].value
                    line_row['{0}FittedWidth'.format(name)] = sp_sl.specfit.parinfo['WIDTH0'].value
                    line_row['{0}FittedAmplitudeError'.format(name)] = sp_sl.specfit.parinfo['AMPLITUDE0'].error
                    line_row['{0}FittedCenterError'.format(name)] = sp_sl.specfit.parinfo['SHIFT0'].error
                    line_row['{0}FittedWidthError'.format(name)] = sp_sl.specfit.parinfo['WIDTH0'].error

    ymin = 0
    ymax = 0
    for ii in range(1,plotnum):
        ax = pl.subplot(7,7,ii)
        ymin = min(ax.get_ylim()[0], ymin)
        ymax = max(ax.get_ylim()[1], ymax)

    for ii in range(1,plotnum):
        ax = pl.subplot(7,7,ii)
        #ax.set_ylim(ymin, ymax)
        #if ii % 7 != 1:
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_label("")
        if ii < plotnum-7:
            ax.xaxis.set_ticklabels([])
            ax.xaxis.set_label("")

    pl.savefig(paths.fpath('spectral_fits/{0}_linefits.png'.format(name)),
               bbox_inches='tight', bbox_extra_artists=[])

    line_table.write(paths.tpath('spectral_lines_and_fits.csv'),
                     overwrite=True)
