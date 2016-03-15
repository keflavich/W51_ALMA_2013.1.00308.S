import numpy as np
import paths
import pyspeckit
import pylab as pl
from astropy import constants
from astropy.io import fits
from astropy.table import Table,Column
from astropy import units as u
from astropy import log
import resource
import spectral_overlays

#lines_to_overlay = ['OCS','H2CO', 'HNCO', 'SO']
frequencies = {'H2CO303_202': 218.22219*u.GHz,
               'H2CO321_220': 218.76007*u.GHz,
               'H2CO322_221': 218.47564*u.GHz,
               'CH3OH422-312': 218.44005*u.GHz,
               'HC3N24-23': 218.32471*u.GHz,
               'OCS18-17': 218.90336*u.GHz,
               'OCS19-18': 231.06099*u.GHz,
               'SO65-54': 219.94944*u.GHz,
               'HNCO10110-919': 218.98102*u.GHz,
               'HNCO1028-927': 219.73719*u.GHz,
               'CH3OH423-514': 234.68345*u.GHz,
               'CH3OH5m42-6m43': 234.69847*u.GHz,
               'CH3OH808-716': 220.07849*u.GHz,
               '13CS5-4': 231.22069*u.GHz,
               'CH3OCH3_13013-12112': 231.98792*u.GHz,
               'NH2CHO11210-1029': 232.27363*u.GHz,
               'NH2CHO1156-1055': 233.59451*u.GHz,
               'HC3Nv7=124-23': 219.17358*u.GHz,
               #'H30alpha': 231.90093*u.GHz,
               #'C18O2-1': 219.56036*u.GHz,
              }
freq_name_mapping = {v:k for k,v in frequencies.items()}
yoffset = {'H2CO303_202': 0,
           'H2CO321_220': 1,
           'H2CO322_221': 2,
           'OCS18-17': 3,
           'SO65-54': 4,
           'CH3OH423-514': 5,
           'CH3OH5m42-6m43': 6,
           'OCS19-18': 7,
           '13CS5-4': 8,
           'CH3OCH3_13013-12112': 9,
           'HNCO1028-927': 10,
           'HNCO10110-919': 11,
           'HC3N24-23': 12,
           'HC3Nv7=124-23': 13,
           'NH2CHO11210-1029': 14,
           'NH2CHO1156-1055': 15,
           'CH3OH422-312': 16,
           'CH3OH808-716': 17,
           #'H30alpha': 4.5,
           #'C18O2-1': 3.5,
          }

pruned_ppcat = Table.read(paths.tpath("dendrogram_continuum_catalog.ipac"),
                          format='ascii.ipac')
dendromask = fits.getdata(paths.apath('dendrograms_min1mJy_diff1mJy_mask_pruned.fits'))

minvelo = 45*u.km/u.s
maxvelo = 90*u.km/u.s

data = {}

for row in pruned_ppcat:
    name = row['_idx']

    fn = paths.spath("dendro{name:03d}_spw{ii}_mean.fits")

    data[name] = spectral_overlays.spectral_overlays(fn, name=name,
                                                     freq_name_mapping=freq_name_mapping,
                                                     frequencies=frequencies,
                                                     yoffset=yoffset,
                                                     minvelo=minvelo,
                                                     maxvelo=maxvelo)

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



# now for the merged data
data = {}

for row in pruned_ppcat:
    name = row['_idx']
    data[name] = {}

    fn = paths.merge_spath("dendro{name:03d}_spw{ii}_mean_7m12m.fits")

    data[name] = spectral_overlays.spectral_overlays(fn, name=name,
                                                     freq_name_mapping=freq_name_mapping,
                                                     frequencies=frequencies,
                                                     yoffset=yoffset,
                                                     minvelo=minvelo,
                                                     maxvelo=maxvelo,
                                                     suffix="_7m12m",
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
