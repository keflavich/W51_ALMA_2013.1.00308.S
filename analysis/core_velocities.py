raise DeprecationWarning("Obsoleted by core_velocities2.py")
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

from line_parameters import frequencies, freq_name_mapping, yoffset

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
    assert not np.isnan(scaling)
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
