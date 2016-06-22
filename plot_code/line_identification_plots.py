import glob
import pyspeckit
from astroquery.splatalogue import Splatalogue
from astropy import units as u
import paths
import pylab as pl


plot_kwargs = {'color':'r', 'linestyle':'--'}
annotate_kwargs = {'color': 'r'}


spectra_to_species = {'e8mm': [('vinyl_cyanide', 'CH2CHCN', 500),
                               ('ethyl_cyanide', 'CH3CH2CN', 500),
                               ('cyanoacetylene', 'HC3N', 1500),
                               ('formamide', 'NH2CHO', 500),
                              ],
                      'e2e': [('formamide', 'NH2CHO', 500),
                              ('cyanoacetylene', 'HC3N', 1500),
                             ],
                      'e2nw': [('formamide', 'NH2CHO', 500),
                               ('methylformate', 'CH3OCHO', 1500),
                               ('dimethylether', 'CH3OCH3', 1500),
                             ],
                      'ALMAmm14':[('methanol', 'CH3OH', 500),],

                     }
snu_min = {'e8mm': 0.3,
           'e2e': 0.3,
           'ALMAmm14': 0.0,
           'north': 0.3,
           'e2nw': 0.1,
          }
velo = {'e8mm': 61*u.km/u.s,
        'e2e': 56*u.km/u.s,
        'e2nw': 55.626*u.km/u.s,
        'ALMAmm14': 62*u.km/u.s,
        'north': 55*u.km/u.s,
       }

pl.figure(1).clf()

for target,species_list in spectra_to_species.items():
    spectra = pyspeckit.Spectra(glob.glob(paths.spath("*{0}*fits".format(target))))

    for species_tuple in species_list:
        species_name, chemid, tmax = species_tuple

        for ii in range(4):
            cat = Splatalogue.query_lines(spectra[ii].xarr.min(),
                                          spectra[ii].xarr.max(),
                                          chemical_name=chemid,
                                          energy_max=tmax,
                                          energy_type='eu_k', noHFS=True,
                                          line_lists=['SLAIM'])
            spectra[ii].plotter(figure=pl.figure(1))
            spectra[ii].plotter.axis.set_ylim(snu_min[target], spectra[ii].plotter.axis.get_ylim()[1])
            spectra[ii].plotter.line_ids(cat['Resolved QNs'],
                                         cat['Freq-GHz']*u.GHz,
                                         velocity_offset=velo[target],
                                         plot_kwargs=plot_kwargs,
                                         annotate_kwargs=annotate_kwargs)
            spectra[ii].plotter.savefig(paths.fpath('line_id_spectra/{target}_{species_name}_{chemid}_spw{0}.png'.format(ii,
                                                                                                                         target=target,
                                                                                                                         chemid=chemid,
                                                                                                                         species_name=species_name,
                                                                                                                        )),
                                        bbox_extra_artists=[])
