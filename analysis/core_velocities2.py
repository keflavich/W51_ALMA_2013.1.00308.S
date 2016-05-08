"""
Determine velocities *and* peak_TB's for the ds9 region cores

This is the 2nd version of core_velocities updated to use the tools created for
dendro_core_velocities
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
import pyregion

cores = pyregion.open(paths.rpath('cores.reg'))

minvelo = 45*u.km/u.s
maxvelo = 90*u.km/u.s

# # array merged
# data = {}
# 
# for corereg in cores:
#     name = corereg.attr[1]['text']
#     data[name] = {}
# 
#     fn = paths.merge_spath("{name}_spw{ii}_mean_7m12m.fits")
#     bgfn = paths.merge_spath("{name}_spw{ii}_background_mean_7m12m.fits")
# 
#     data[name] = spectral_overlays.spectral_overlays(fn, name=name,
#                                                      freq_name_mapping=freq_name_mapping,
#                                                      frequencies=frequencies,
#                                                      yoffset=yoffset,
#                                                      minvelo=minvelo,
#                                                      maxvelo=maxvelo,
#                                                      suffix="_7m12m",
#                                                      background_fn=bgfn,
#                                                     )
# 
# firstentry = list(data.keys())[0]
# colnames = list(data[firstentry].keys())
# coltypes = {k:type(data[firstentry][k]) for k in colnames}
# names = Column([name for name in data], name='SourceID')
# data_list = [Column(u.Quantity([data[name][key] for name in names]), name=key)
#              if coltypes[key] not in (str,)
#              else Column([data[name][key] for name in names], name=key)
#              for key in colnames]
# data_list.insert(0, names)
# 
# 
# tbl = Table(data_list)
# tbl.sort('SourceID')
# tbl.write(paths.tpath("core_velocities_7m12mspectra.ipac"), format="ascii.ipac")


# 12m only
data = {}

for corereg in cores:
    name = corereg.attr[1]['text']

    fn = paths.spath("{name}_spw{ii}_mean.fits")

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
tbl.write(paths.tpath("core_velocities.ipac"), format="ascii.ipac")
