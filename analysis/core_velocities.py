import numpy as np
import paths
import pyspeckit
import pylab as pl
from astropy import coordinates, constants
from astropy.table import Table,Column
from astropy import units as u
import glob
from astropy import log
import pyregion

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

cores = pyregion.open(paths.rpath('cores.reg'))

minvelo = 45*u.km/u.s
maxvelo = 90*u.km/u.s

data = {}

for corereg in cores:
    name = corereg.attr[1]['text']
    data[name] = {}

    fn = "{name}_spw{ii}_mean.fits"
    spectra = pyspeckit.Spectra([paths.spath(fn.format(name=name, ii=ii))
                                 for ii in range(4)])
    spectra.data[(233.84*u.GHz<spectra.xarr) & (spectra.xarr>234.036*u.GHz)] = np.nan
    spectra.data[(230.00*u.GHz<spectra.xarr) & (spectra.xarr>230.523*u.GHz)] = np.nan
    scaling = np.nanmax(spectra.data) - np.nanpercentile(spectra.data, 20)
    print("Scaling for {fn} = {scaling}".format(fn=fn.format(name=name, ii=0), scaling=scaling))

    fig = pl.figure(1)
    fig.clf()

    for spwnum,sp in enumerate(spectra):
        if spwnum == 3:
            # flag out the middle section where apparently many antennae have been flagged
            sp.data[1600:1984] = np.nan
        elif spwnum == 2:
            # ignore 12CO if it's brightest: it's never representative
            # (and it probably screws up the plot scale if it's around)
            sp.data[:200] = np.nan
        peak = np.nanmax(sp.data)
        argmax = np.nanargmax(sp.data)
        cont = np.nanpercentile(sp.data, 20)

        sp.data -= cont
        sp.xarr.convert_to_unit(u.GHz)

        data[name]['continuum20pct'] = cont
        freqlist = list(freq_name_mapping.keys())
        peakfreq = sp.xarr[argmax]
        data[name]['peak{0}freq'.format(spwnum)] = peakfreq
        bestmatch = np.argmin(np.abs(peakfreq - u.Quantity(freqlist)))
        closest_freq = freqlist[bestmatch]

        peakvelo = ((closest_freq-peakfreq)/closest_freq *
                    constants.c).to(u.km/u.s)
        velo_OK = (minvelo < peakvelo) and (peakvelo < maxvelo)

        data[name]['peak{0}velo'.format(spwnum)] = peakvelo if velo_OK else np.nan*u.km/u.s
        peakspecies = (freq_name_mapping[closest_freq] if velo_OK else 'none')
        data[name]['peak{0}species'.format(spwnum)] = peakspecies
        data[name]['peak{0}'.format(spwnum)] = (peak if velo_OK else np.nan)*u.Jy/u.beam

        log.debug("spw{0} peak{0}={1} line={2}".format(spwnum,
                                                       peak,
                                                       peakspecies)
                 )

        for linename, freq in frequencies.items():
            if sp.xarr.in_range(freq):
                print("Plotting {0} for {1} from spw {2}".format(linename, name, spwnum))
                temp = np.array(sp.xarr)
                sp.xarr.convert_to_unit(u.km/u.s, refX=freq)
                sp.plotter(figure=fig, clear=False, offset=yoffset[linename]*scaling,
                           xmin=0, xmax=120)
                sp.plotter.axis.set_ylim(-0.1, max(yoffset.values())*scaling+scaling)
                sp.plotter.axis.text(122, yoffset[linename]*scaling, linename)
                sp.xarr.convert_to_unit(u.GHz, refX=freq)
                np.testing.assert_allclose(temp, np.array(sp.xarr))

    okvelos = [data[name]['peak{0}velo'.format(ii)] for ii in range(4) if
               (minvelo < data[name]['peak{0}velo'.format(ii)]) and
               (data[name]['peak{0}velo'.format(ii)] < maxvelo)]
    if okvelos:
        velo = np.mean(u.Quantity(okvelos))
        data[name]['mean_velo'] = velo
    else:
        data[name]['mean_velo'] = np.nan*u.km/u.s

    fig.savefig(paths.fpath("spectral_overlays/{name}_overlaid_spectra.png".format(name=name)))

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
