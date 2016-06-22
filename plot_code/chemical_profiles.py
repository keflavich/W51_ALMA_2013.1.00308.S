import paths
import pylab as pl
import os
from spectral_cube import SpectralCube
import numpy as np
from astropy import units as u

import re
import glob

linere = re.compile("W51_b6_7M_12M.(.*).image.pbcor")

for fignum,(region,xslc,yslc,velo) in enumerate((('north',209,slice(65,175),60),
                                                 ('e2e8',slice(95,220),411,57.4),
                                                 ('e2e8',slice(80,260),282,57.4),
                                                 #('ALMAmm14',slice(50,150),20,58),
                                                )):
    fig = pl.figure(fignum)
    fig.clf()
    ax = pl.gca()
    ax.set_title(region)
    for fn in glob.glob(paths.dpath("merge/cutouts/W51_b6_7M_12M.*{0}*cutout.fits".format(region))):
        cube = SpectralCube.read(fn)
        slc = cube.filled_data[cube.closest_spectral_channel(velo*u.km/u.s), yslc, xslc]

        label = linere.search(fn).groups()[0]
        if 'CH3OH' in label:
            ax.plot(slc, label=label, linewidth=2, color='b', alpha=0.5)
        elif 'H2CO' in label:
            ax.plot(slc, label=label, linewidth=2, color='r', alpha=0.5)
        else:
            ax.plot(slc, label=label, alpha=0.5)

    pl.legend(loc='best')
