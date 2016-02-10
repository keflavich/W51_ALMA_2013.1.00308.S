import numpy as np
import paths
import pyspeckit
import pylab as pl
from astropy import coordinates, constants
from astropy.table import Table,Column
from astropy import units as u
import glob
import pyregion

#lines_to_overlay = ['OCS','H2CO', 'HNCO', 'SO']
frequencies = {'H2CO303_202': 218.22219*u.GHz,
               'H2CO321_220': 218.76007*u.GHz,
               'H2CO322_221': 218.47564*u.GHz,
               'OCS18-17': 218.90336*u.GHz,
               'OCS19-18': 231.06099*u.GHz,
               'SO65-54': 219.94944*u.GHz,
               'CH3OH423-514': 234.68345*u.GHz,
               'CH3OH5m42-6m43': 234.69847*u.GHz,
               '13CS5-4': 231.22069*u.GHz,
               'CH3OCH3_13013-12112': 231.98792*u.GHz,
              }
freq_name_mapping = {v:k for k,v in frequencies.items()}
yoffset = {'H2CO303_202': 0.4,
           'H2CO321_220': 0.3,
           'H2CO322_221': 0.2,
           'OCS18-17': 0.1,
           'SO65-54': 0.0,
           'CH3OH423-514': 0.5,
           'CH3OH5m42-6m43': 0.6,
           'OCS19-18': 0.7,
           '13CS5-4': 0.8,
           'CH3OCH3_13013-12112': 0.9,
          }

cores = pyregion.open(paths.rpath('cores.reg'))

data = {}

for corereg in cores:
    name = corereg.attr[1]['text']
    data[name] = {}

    fn = "{name}_spw{ii}_mean.fits"
    spectra = [pyspeckit.Spectrum(fn.format(name=name, ii=ii))
               for ii in range(4)]

    fig = pl.figure(1)
    fig.clf()

    for spwnum,sp in enumerate(spectra):
        if spwnum == 3:
            # flag out the middle section where apparently many antennae have been flagged
            sp.data[1618:1984] = np.nan
        peak = np.nanmax(sp.data)
        argmax = np.nanargmax(sp.data)
        cont = np.nanpercentile(sp.data, 20)

        sp.data -= cont
        sp.xarr.convert_to_unit(u.GHz)

        data[name]['continuum20pct'] = cont
        data[name]['peak{0}'.format(spwnum)] = peak
        data[name]['peak{0}freq'.format(spwnum)] = peakfreq = sp.xarr[argmax]
        freqlist = list(freq_name_mapping.keys())
        bestmatch = np.argmin(np.abs(peakfreq - u.Quantity(freqlist)))
        closest_freq = freqlist[bestmatch]
        data[name]['peak{0}velo'.format(spwnum)] = ((closest_freq-peakfreq)/closest_freq * constants.c).to(u.km/u.s)
        data[name]['peak{0}species'.format(spwnum)] = freq_name_mapping[closest_freq]


        for linename, freq in frequencies.items():
            if sp.xarr.in_range(freq):
                print("Plotting {0} for {1} from spw {2}".format(linename, name, spwnum))
                temp = np.array(sp.xarr)
                sp.xarr.convert_to_unit(u.km/u.s, refX=freq)
                sp.plotter(figure=fig, clear=False, offset=yoffset[linename],
                           xmin=0, xmax=120)
                sp.plotter.axis.set_ylim(-0.1, max(yoffset.values())+0.1)
                sp.plotter.axis.text(122, yoffset[linename], linename)
                sp.xarr.convert_to_unit(u.GHz, refX=freq)
                np.testing.assert_allclose(temp, np.array(sp.xarr))

    okvelos = [data[name]['peak{0}velo'.format(ii)] for ii in range(4) if
               (45*u.km/u.s < data[name]['peak{0}velo'.format(ii)]) and
               (data[name]['peak{0}velo'.format(ii)] < 75*u.km/u.s)]
    if okvelos:
        velo = np.mean(u.Quantity(okvelos))
        data[name]['mean_velo'] = velo
    else:
        data[name]['mean_velo'] = np.nan*u.km/u.s

    fig.savefig("{name}_overlaid_spectra.png".format(name=name))

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
