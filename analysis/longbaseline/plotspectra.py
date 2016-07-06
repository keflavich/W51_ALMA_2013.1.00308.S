import os
import glob
import pyspeckit
from line_to_image_list import line_to_image_list
from astropy import units as u
import paths
import pylab as pl

plot_kwargs = {'color':'r', 'linestyle':'--'}
annotate_kwargs = {'color': 'r'}

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

        fpre = os.path.splitext(os.path.split(fn)[-1])[0]

        species_names = [x[0] for x in line_to_image_list]
        frequencies = u.Quantity([float(x[1].strip("GHz")) for x in line_to_image_list],
                                 unit=u.GHz)

        sp.plotter(figure=pl.figure(1))

        sp.plotter.line_ids(species_names, u.Quantity(frequencies),
                            velocity_offset=velo[target],
                            plot_kwargs=plot_kwargs,
                            annotate_kwargs=annotate_kwargs)
        outname = paths.dpath('longbaseline/spectra/pngs/{target}_{0}.png'.format(fpre,
                                                                                  target=target)
                                                                                                                     )
        sp.plotter.savefig(outname, bbox_extra_artists=[])
