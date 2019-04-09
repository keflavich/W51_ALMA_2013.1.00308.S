import numpy as np
import os
import glob
import pyspeckit
from line_to_image_list import line_to_image_list
from astropy import units as u
import paths
import pylab as pl
from astroquery.splatalogue import Splatalogue, utils

plot_kwargs = {'color':'r', 'linestyle':'--'}
annotate_kwargs = {'color': 'r'}
pl.close(1)

velo = {'ALMAmm24_W51n': 60*u.km/u.s,
        'd2_W51n': 60*u.km/u.s,
        'e2e_W51e2': 60*u.km/u.s,
        'e2nw_W51e2': 60*u.km/u.s,
        'e2w_W51e2': 60*u.km/u.s,
        'e8_W51e2': 60*u.km/u.s,
        'north_W51n': 60*u.km/u.s,
       }

for target in velo:
    files = glob.glob(paths.dpath("longbaseline/spectra/{0}*.fits".format(target)))
    for fn in files:
        sp = pyspeckit.Spectrum(fn)
        sp.xarr.convert_to_unit(u.GHz)
        print(fn, sp.xarr.min(), sp.xarr.max())

        fpre = os.path.splitext(os.path.split(fn)[-1])[0]

        species_names = [x[0] for x in line_to_image_list]
        frequencies = u.Quantity([float(x[1].strip("GHz")) for x in line_to_image_list],
                                 unit=u.GHz)

        splat = Splatalogue.query_lines(sp.xarr.min(), sp.xarr.max(),
                                        chemical_name='NaCl|KCl',
                                        line_lists=['CDMS'])
        if len(splat) > 0:
            splat = utils.minimize_table(splat)
            species_names = [row['Species'] + row['QNs'] for row in splat]
            frequencies = u.Quantity(splat['Freq'], u.GHz)
        else:
            species_names = []
            frequencies = []

        sp.plotter(figure=pl.figure(1))
        sp.plotter(figure=pl.figure(1))

        outname = paths.dpath('longbaseline/spectra/pngs/{target}_{0}.png'
                              .format(fpre, target=target))
        sp.plotter.savefig(outname, bbox_extra_artists=[])

        okfreqs = (frequencies > sp.xarr.min()) & (frequencies < sp.xarr.max())
        #print(list(zip(species_names, frequencies)))
        #print(frequencies[okfreqs])

        sp.plotter.line_ids(np.array(species_names)[okfreqs],
                            u.Quantity(frequencies)[okfreqs],
                            velocity_offset=velo[target],
                            plot_kwargs=plot_kwargs,
                            annotate_kwargs=annotate_kwargs)
        outname = paths.dpath('longbaseline/spectra/pngs/labeled_{target}_{0}.png'
                              .format(fpre, target=target))
        sp.plotter.savefig(outname, bbox_extra_artists=[])
