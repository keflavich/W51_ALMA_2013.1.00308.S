import paths
import pylab as pl
import os
import itertools
from spectral_cube import SpectralCube
import numpy as np
from astropy import units as u

import re
import glob

linere = re.compile("W51_b6_7M_12M.(.*).image.pbcor")

for fignum,(region,suffix,xslc,yslc,velo) in enumerate((('north','',209,slice(65,175),60),
                                                       ('e2e8','_1',slice(95,220),411,57.4),
                                                       ('e2e8','_2',slice(80,260),282,57.4),
                                                       #('ALMAmm14',slice(50,150),20,58),
                                                      )):
    fig = pl.figure(fignum)
    fig.clf()
    ax = pl.gca()
    ax.set_title(region)
    path_template = paths.dpath("merge/cutouts/W51_b6_7M_12M.*{0}*cutout.fits".format(region))

    slices = {}

    for ii,fn in enumerate(glob.glob(path_template)):
        label = linere.search(fn).groups()[0]
    
        if 'PN' in label or '13CH3OH515' in label:
            continue

        cube = SpectralCube.read(fn)
        slc = cube.filled_data[cube.closest_spectral_channel(velo*u.km/u.s), yslc, xslc]

        slices[label] = slc

    #sort_order = sorted(slices.keys(), key=(slc.max().value for slc in slices.values()))
    sort_order = sorted(slices.items(), key=lambda x: x[1].max().value)

    colors = itertools.cycle(('c','m','k','y'))
    linestyles = itertools.cycle(('--','-',':','-.'))

    for ii,(label,slc) in enumerate(sort_order):

        offset = ii*0.3*slc.unit - slc.min()

        if 'CH3OH' in label:
            ax.plot(slc + offset, label=label, linewidth=2, color='b', alpha=0.5,
                    linestyle=next(linestyles))
        elif 'H2CO' in label:
            ax.plot(slc + offset, label=label, linewidth=2, color='r', alpha=0.5,
                    linestyle=next(linestyles))
        elif 'CH3O' in label:
            ax.plot(slc + offset, label=label, linewidth=2, color='g', alpha=0.5,
                    linestyle=next(linestyles))
        elif 'PN' in label:
            continue
        else:
            ax.plot(slc + offset, label=label, alpha=0.5, color=next(colors),
                    linestyle=next(linestyles))

    pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    pl.savefig(paths.fpath('chemslices/{0}{1}_chemslice.png'.format(region,suffix)),
               bbox_inches='tight')
