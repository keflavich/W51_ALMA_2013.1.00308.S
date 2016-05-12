import radio_beam
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

    shift = (minvelo+maxvelo)/2. / constants.c

    sp.data -= cont
    sp.xarr.convert_to_unit(u.GHz)
    peak = np.nanmax(sp.data)
    peakfreq = sp.xarr[argmax]
    assert sp.data[argmax] == peak
    peakfreq_shifted = peakfreq * (1+shift)
    freqlist = list(freq_name_mapping.keys())
    reverse_freq_name_mapping = {v:k for k,v in freq_name_mapping.items()}
    bestmatch = np.argmin(np.abs(peakfreq_shifted - u.Quantity(freqlist)))
    closest_freq = freqlist[bestmatch]
    peakvelo = ((closest_freq-peakfreq)/closest_freq *
                constants.c).to(u.km/u.s)
    velo_OK = (minvelo < peakvelo) and (peakvelo < maxvelo)
    peakspecies = (freq_name_mapping[closest_freq] if velo_OK else 'none')

    return (cont, peak, peakfreq, peakfreq_shifted, bestmatch, peakvelo,
            velo_OK, peakspecies, argmax)




def spectral_overlays(fn, name, freq_name_mapping, frequencies, yoffset,
                      minvelo, maxvelo, suffix="", background_fn=None,
                      return_spectra=False, plot_fullspec=True):

    object_data_dict = {}

    spectra = pyspeckit.Spectra([fn.format(name=name, ii=ii)
                                 for ii in range(4)])
    bad_1 = (233.74*u.GHz<spectra.xarr.to(u.GHz)) & (spectra.xarr.to(u.GHz)<234.036*u.GHz)
    bad_2 = (230.00*u.GHz<spectra.xarr.to(u.GHz)) & (spectra.xarr.to(u.GHz)<230.523*u.GHz)
    bad_3 = (spectra.xarr.to(u.GHz) < 218.11*u.GHz)
    spectra.data[bad_1 | bad_2 | bad_3] = np.nan
    beams = [radio_beam.Beam.from_fits_header(sp.header) for sp in spectra]

    # scaling: determine how much to separate spectra by vertically
    scaling = np.nanmax(spectra.data) - np.nanpercentile(spectra.data, 20)
    print("Scaling for {fn} = {scaling}".format(fn=fn.format(name=name, ii=0),
                                                scaling=scaling))
    if np.isnan(scaling):
        raise ValueError("All-nan slice encountered.  There is apparently no data in this file?")

    if background_fn is not None:
        bgspectra = pyspeckit.Spectra([background_fn.format(name=name, ii=ii)
                                       for ii in range(4)])
        bgspectra.data[bad_1 | bad_2 | bad_3] = np.nan
        bg = True
    else:
        bg = False

    fig = pl.figure(0)
    fig.clf()

    for spwnum,sp in enumerate(spectra):
        bad_1 = (233.74*u.GHz<sp.xarr.to(u.GHz)) & (sp.xarr.to(u.GHz)<234.036*u.GHz)
        bad_2 = (230.00*u.GHz<sp.xarr.to(u.GHz)) & (sp.xarr.to(u.GHz)<230.523*u.GHz)
        bad_3 = (sp.xarr.to(u.GHz) < 218.11*u.GHz)
        sp.data[bad_1 | bad_2 | bad_3] = np.nan
        
        # temporary hack for bad data
        if bg and all(bgspectra[spwnum].data == 0):
            bg = False

        (cont, peak, peakfreq, peakfreq_shifted, bestmatch, peakvelo, velo_OK,
         peakspecies, argmax) = quick_analyze(sp, freq_name_mapping, minvelo, maxvelo)

        object_data_dict['continuum20pct{0}'.format(spwnum)] = cont
        object_data_dict['peak{0}freq'.format(spwnum)] = peakfreq

        object_data_dict['peak{0}velo'.format(spwnum)] = peakvelo if velo_OK else np.nan*u.km/u.s
        object_data_dict['peak{0}species'.format(spwnum)] = peakspecies
        object_data_dict['peak{0}'.format(spwnum)] = (peak if velo_OK else np.nan)*u.Jy/u.beam
        object_data_dict['beam{0}area'.format(spwnum)] = beams[spwnum].sr.value

        if bg:
            (bgcont, bgpeak, bgpeakfreq, bgpeakfreq_shifted, bgbestmatch,
             bgpeakvelo, bgvelo_OK, bgpeakspecies,
             bgargmax) = quick_analyze(bgspectra[spwnum], freq_name_mapping,
                                       minvelo, maxvelo)
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
        velo = 60*u.km/u.s

    fig.savefig(paths.fpath("spectral_overlays/{name}_overlaid_spectra{suffix}.png".format(name=name, suffix=suffix)),
                bbox_inches='tight')

    # plot type #2: full spectrum, with lines ID'd
    if plot_fullspec:
        linenames = list(frequencies.keys())
        freqs_ghz = list(frequencies.values())
        plot_kwargs = {'color':'r', 'linestyle':'--'}
        annotate_kwargs = {'color': 'r'}

        for spwnum,sp in enumerate(spectra):
            fig = pl.figure(spwnum+1)
            fig.clf()
            sp.xarr.convert_to_unit(u.GHz)
            sp.plotter(figure=fig, axis=fig.gca())
            sp.plotter.line_ids(linenames,
                                u.Quantity(freqs_ghz),
                                velocity_offset=velo,
                                plot_kwargs=plot_kwargs,
                                annotate_kwargs=annotate_kwargs)

            if bg:
                bgs = bgspectra[spwnum]
                bgs.xarr.convert_to_unit(u.GHz)
                bgs.plotter(axis=sp.plotter.axis,
                                          clear=False,
                                          color='b',
                                          zorder=-100,
                                         )

            fig.savefig(paths.fpath("spectral_overlays/{name}_spw{spw}_fullspec{suffix}.png".format(name=name,
                                                                                                    spw=spwnum,
                                                                                                    suffix=suffix)),
                        bbox_inches='tight')

            

    if return_spectra:
        return object_data_dict, spectra
    else:
        return object_data_dict
