"""
Given the chemical maps made from chem_images, perform slice plots through them
"""
import paths
import pylab as pl
import os
import itertools
from astropy.io import fits
import numpy as np
from astropy import units as u

import re
import glob

linere = re.compile("chemical_m0_slabs_[^_]*_(.*?)(_merge.fits|.fits)")

for slicetype in ('m0','max'):
    for fignum,(region,suffix,xslc,yslc,velo) in enumerate((('north','',51,slice(None),60),
                                                           ('e2','_1',slice(None),45,57.4),
                                                           ('e2','_2',53,slice(None),57.4),
                                                           ('e8','_1',slice(None),58,57.4),
                                                           ('e8','_2',67,slice(None),57.4),
                                                           #('ALMAmm14',slice(50,150),20,58),
                                                          )):
        fig = pl.figure(fignum)
        fig.clf()
        ax = pl.gca()
        ax.set_title(region)
        path_template = paths.dpath("chemslices/chemical_{1}_slabs*_{0}_*.fits".format(region, slicetype))

        slices = {}

        for ii,fn in enumerate(glob.glob(path_template)):
            if 'natural' in fn or 'merge' in fn:
                continue

            label = linere.search(fn).groups()[0]
            print(fn, label)
            assert 'merge' not in label
        
            if 'PN' in label or '13CH3OH515' in label:
                continue

            slc = fits.getdata(fn)

            slices[label] = slc[yslc,xslc]

        #sort_order = sorted(slices.keys(), key=(slc.max().value for slc in slices.values()))
        assert slices
        sort_order = sorted(slices.items(), key=lambda x: x[1].max())

        colors = itertools.cycle(('c','m','k','y'))
        linestyles = itertools.cycle(('--','-',':','-.'))

        for ii,(label,slc) in enumerate(sort_order):

            print(label)

            offset = ii*30 - slc.min()

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

        pl.savefig(paths.fpath('chemslices/{0}{2}{1}_chemslice.png'.format(region,suffix,slicetype)),
                   bbox_inches='tight')
