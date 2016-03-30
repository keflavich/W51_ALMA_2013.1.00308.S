import numpy as np
import paths
import pyspeckit
import pylab as pl
from astropy import constants
from astropy import units as u
from astropy import log

def quick_analyze(sp, freq_name_mapping, minvelo, maxvelo):
    """ get peak of spectrum, subtract continuum, etc. """
    argmax = np.nanargmax(sp.data)
    cont = np.nanpercentile(sp.data, 20)

    sp.data -= cont
    sp.xarr.convert_to_unit(u.GHz)
    peak = np.nanmax(sp.data)
    peakfreq = sp.xarr[argmax]
    freqlist = list(freq_name_mapping.keys())
    bestmatch = np.argmin(np.abs(peakfreq - u.Quantity(freqlist)))
    closest_freq = freqlist[bestmatch]
    peakvelo = ((closest_freq-peakfreq)/closest_freq *
                constants.c).to(u.km/u.s)
    velo_OK = (minvelo < peakvelo) and (peakvelo < maxvelo)
    peakspecies = (freq_name_mapping[closest_freq] if velo_OK else 'none')

    return cont, peak, peakfreq, bestmatch, peakvelo, velo_OK, peakspecies, argmax




def spectral_overlays(fn, name, freq_name_mapping, frequencies, yoffset,
                      minvelo, maxvelo, suffix="", background_fn=None):

    object_data_dict = {}

    spectra = pyspeckit.Spectra([fn.format(name=name, ii=ii)
                                 for ii in range(4)])
    bad_1 = (233.84*u.GHz<spectra.xarr.to(u.GHz)).value & (spectra.xarr.to(u.GHz)>234.036*u.GHz).value
    bad_2 = (230.00*u.GHz<spectra.xarr.to(u.GHz)).value & (spectra.xarr.to(u.GHz)>230.523*u.GHz).value
    spectra.data[bad_1] = np.nan
    spectra.data[bad_2] = np.nan

    # scaling: determine how much to separate spectra by vertically
    scaling = np.nanmax(spectra.data) - np.nanpercentile(spectra.data, 20)
    print("Scaling for {fn} = {scaling}".format(fn=fn.format(name=name, ii=0),
                                                scaling=scaling))
    if np.isnan(scaling):
        raise ValueError("All-nan slice encountered.  There is apparently no data in this file?")

    if background_fn is not None:
        bgspectra = pyspeckit.Spectra([background_fn.format(name=name, ii=ii)
                                     for ii in range(4)])
        bgspectra.data[(233.84*u.GHz<bgspectra.xarr).value &
                       (bgspectra.xarr>234.036*u.GHz).value] = np.nan
        bgspectra.data[(230.00*u.GHz<bgspectra.xarr).value &
                       (bgspectra.xarr>230.523*u.GHz).value] = np.nan
        bg = True
    else:
        bg = False

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
        
        # temporary hack for bad data
        if bg and all(bgspectra[spwnum].data == 0):
            bg = False

        (cont, peak, peakfreq, bestmatch, peakvelo, velo_OK,
         peakspecies, argmax) = quick_analyze(sp, freq_name_mapping, minvelo, maxvelo)

        object_data_dict['continuum20pct'] = cont
        object_data_dict['peak{0}freq'.format(spwnum)] = peakfreq

        object_data_dict['peak{0}velo'.format(spwnum)] = peakvelo if velo_OK else np.nan*u.km/u.s
        object_data_dict['peak{0}species'.format(spwnum)] = peakspecies
        object_data_dict['peak{0}'.format(spwnum)] = (peak if velo_OK else np.nan)*u.Jy/u.beam

        if bg:
            (bgcont, bgpeak, bgpeakfreq, bgbestmatch, bgpeakvelo, bgvelo_OK,
             bgpeakspecies, bgargmax) = quick_analyze(bgspectra[spwnum], freq_name_mapping, minvelo, maxvelo)
            object_data_dict['bgpeak{0}freq'.format(spwnum)] = bgpeakfreq

            object_data_dict['bgpeak{0}velo'.format(spwnum)] = bgpeakvelo if bgvelo_OK else np.nan*u.km/u.s
            object_data_dict['bgpeak{0}species'.format(spwnum)] = bgpeakspecies
            object_data_dict['bgpeak{0}'.format(spwnum)] = (bgpeak if bgvelo_OK else np.nan)*u.Jy/u.beam


        log.debug("spw{0} peak{0}={1} line={2}".format(spwnum,
                                                       peak,
                                                       peakspecies)
                 )

        for linename, freq in frequencies.items():
            if sp.xarr.in_range(freq):
                print("Plotting {0} for {1} from spw {2}.  Peakspecies={3}".format(linename, name, spwnum, peakspecies))
                temp = np.array(sp.xarr)
                sp.xarr.convert_to_unit(u.km/u.s, refX=freq)
                if np.isnan(yoffset[linename]) or np.isnan(scaling):
                    raise ValueError("NAN scaling is stupid.")
                if sp.slice(0*u.km/u.s, 120*u.km/u.s).data.max() > sp.slice(0*u.km/u.s, 120*u.km/u.s).data.min():
                    sp.plotter(figure=fig, clear=False, offset=yoffset[linename]*scaling,
                               xmin=0, xmax=120)
                    sp.plotter.axis.text(122, yoffset[linename]*scaling,
                                         "{0}: {1}".format(spwnum, linename))
                sp.xarr.convert_to_unit(u.GHz, refX=freq)
                np.testing.assert_allclose(temp, np.array(sp.xarr))

                if bg:
                    bgs = bgspectra[spwnum]
                    bgs.xarr.convert_to_unit(u.km/u.s, refX=freq)
                    if bgs.slice(0*u.km/u.s, 120*u.km/u.s).data.max() > bgs.slice(0*u.km/u.s, 120*u.km/u.s).data.min():
                        bgs.plotter(axis=sp.plotter.axis,
                                                  clear=False,
                                                  color='b',
                                                  offset=yoffset[linename]*scaling,
                                                  zorder=-100,
                                                  xmin=0,
                                                  xmax=120,
                                                 )
                    bgs.xarr.convert_to_unit(u.GHz, refX=freq)

                sp.plotter.axis.set_ylim(-0.1, max(yoffset.values())*scaling+scaling)

                if linename == peakspecies:
                    print('peakvelo, offset, maxdata: ', peakvelo,
                          yoffset[linename]*scaling, sp.data[argmax])
                    sp.plotter.axis.plot(peakvelo, yoffset[linename]*scaling +
                                         sp.data[argmax], 'rx')

    okvelos = [object_data_dict['peak{0}velo'.format(ii)] for ii in range(4) if
               (minvelo < object_data_dict['peak{0}velo'.format(ii)]) and
               (object_data_dict['peak{0}velo'.format(ii)] < maxvelo)]
    if okvelos:
        velo = np.mean(u.Quantity(okvelos))
        object_data_dict['mean_velo'] = velo
        sp.plotter.axis.vlines(velo.to(u.km/u.s).value,
                               sp.plotter.axis.get_ylim()[0],
                               sp.plotter.axis.get_ylim()[1],
                               linestyle='--', linewidth=2, color='r',
                               zorder=-50, alpha=0.5)
    else:
        object_data_dict['mean_velo'] = np.nan*u.km/u.s

    fig.savefig(paths.fpath("spectral_overlays/{name}_overlaid_spectra{suffix}.png".format(name=name, suffix=suffix)),
                bbox_inches='tight')

    return object_data_dict
