import glob
import pyspeckit
from astroquery.splatalogue import Splatalogue
from astropy import units as u
import paths
import pylab as pl
from line_to_image_list import line_to_image_list


plot_kwargs = {'color':'r', 'linestyle':'--'}
annotate_kwargs = {'color': 'r'}


snu_min = {'e8mm': 0.1,
           'e2e': 0.3,
           'ALMAmm14': 0.0,
           'ALMAmm41': -0.005,
           'e2nw': -0.005,
           'e2se': -0.005,
           'north': 0.2,
          }
velo = {'e8mm': 61*u.km/u.s,
        'e2e': 56*u.km/u.s,
        'ALMAmm14': 62*u.km/u.s,
        'ALMAmm41': 55*u.km/u.s,
        'e2nw': 55.626*u.km/u.s,
        'e2se': 54*u.km/u.s,
        'north': 58*u.km/u.s,
       }

pl.figure(1).clf()

for target in snu_min:
    files = glob.glob(paths.spath("*{0}*fits".format(target)))
    if len(files) == 0:
        print("No matches for {0}".format(target))
        continue

    spectra = pyspeckit.Spectra(files)

    for ii in range(4):
        spectra[ii].plotter(figure=pl.figure(1))

        species_names = [x[0] for x in line_to_image_list]
        frequencies = u.Quantity([float(x[1].strip("GHz")) for x in line_to_image_list],
                                 unit=u.GHz)

        spectra[ii].plotter.axis.set_ylim(snu_min[target],
                                          spectra[ii].plotter.axis.get_ylim()[1])
        spectra[ii].plotter.line_ids(species_names,
                                     u.Quantity(frequencies),
                                     velocity_offset=velo[target],
                                     plot_kwargs=plot_kwargs,
                                     annotate_kwargs=annotate_kwargs)
        outname = paths.fpath('line_id_spectra/{target}_spw{0}.png'.format(ii, target=target)
                                                                                                                    )
        spectra[ii].plotter.savefig(outname,
                                    bbox_extra_artists=[])
