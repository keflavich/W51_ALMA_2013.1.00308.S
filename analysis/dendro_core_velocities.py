"""
Determine velocities *and* peak_TB's
"""
import numpy as np
import paths
import pyspeckit
import pylab as pl
from astropy import constants
from astropy.io import fits
from astropy.table import Table,Column
from astropy import units as u
from astropy import log
import spectral_overlays

from line_parameters import frequencies, freq_name_mapping, yoffset

pruned_ppcat = Table.read(paths.tpath("dendrogram_continuum_catalog.ipac"),
                          format='ascii.ipac')
dendromask = fits.getdata(paths.apath('dendrograms_min1mJy_diff1mJy_mask_pruned.fits'))

minvelo = 45*u.km/u.s
maxvelo = 90*u.km/u.s

# array merged
data = {}

for row in pruned_ppcat:
    name = row['_idx']
    data[name] = {}

    fn = paths.merge_spath("dendro{name:03d}_spw{ii}_mean_7m12m.fits")
    bgfn = paths.merge_spath("dendro{name:03d}_spw{ii}_background_mean_7m12m.fits")

    data[name] = spectral_overlays.spectral_overlays(fn, name=name,
                                                     freq_name_mapping=freq_name_mapping,
                                                     frequencies=frequencies,
                                                     yoffset=yoffset,
                                                     minvelo=minvelo,
                                                     maxvelo=maxvelo,
                                                     suffix="_7m12m",
                                                     background_fn=bgfn,
                                                    )

firstentry = list(data.keys())[0]
colnames = list(data[firstentry].keys())
coltypes = {k:type(data[firstentry][k]) for k in colnames}
names = Column([name for name in data], name='SourceID')
data_list = [Column(u.Quantity([data[name][key] for name in names]), name=key)
             if coltypes[key] not in (str,)
             else Column([data[name][key] for name in names], name=key)
             for key in colnames]
data_list.insert(0, names)


tbl = Table(data_list)
tbl.sort('SourceID')
tbl.write(paths.tpath("dendro_core_velocities_7m12mspectra.ipac"), format="ascii.ipac")


# 12m only
data = {}

for row in pruned_ppcat:
    name = row['_idx']

    fn = paths.spath("dendro{name:03d}_spw{ii}_mean.fits")

    result = spectral_overlays.spectral_overlays(fn, name=name,
                                                 freq_name_mapping=freq_name_mapping,
                                                 frequencies=frequencies,
                                                 yoffset=yoffset,
                                                 minvelo=minvelo,
                                                 maxvelo=maxvelo,
                                                )
    if result:
        data[name] = result

firstentry = list(data.keys())[0]
colnames = list(data[firstentry].keys())
coltypes = {k:type(data[firstentry][k]) for k in colnames}
names = Column([name for name in data], name='SourceID')
data_list = [Column(u.Quantity([data[name][key] for name in names]), name=key)
             if coltypes[key] not in (str,)
             else Column([data[name][key] for name in names], name=key)
             for key in colnames]
data_list.insert(0, names)


tbl = Table(data_list)
tbl.sort('SourceID')
tbl.write(paths.tpath("dendro_core_velocities.ipac"), format="ascii.ipac")
