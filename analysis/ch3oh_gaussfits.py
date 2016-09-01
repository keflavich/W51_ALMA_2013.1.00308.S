import numpy as np
import line_to_image_list
from astropy import units as u
from astropy import constants
from astropy import log
import pylab as pl

from ch3oh_rotational_diagram_maps import nupper_of_kkms, fit_tex
from ch3oh_model import show_modelfit

from astroquery.splatalogue import Splatalogue, utils

slaim = Splatalogue.query_lines(218*u.GHz, 235*u.GHz, chemical_name=' CH3OH vt = 0 ',
                                energy_max=3500, energy_type='eu_k',
                                line_lists=['SLAIM'],
                                show_upper_degeneracy=True)
slaimfreqs = np.array(slaim['Freq-GHz'])*u.GHz

def fit_each_line_with_a_gaussian(spectra, ampguess, vguess_kms, widthguess,
                                  fixed_velo=False, fixed_width=False,
                                  wrange_GHz=(0,0.01), plot=False,
                                  separation_tolerance=0.005*u.GHz,
                                  fignum_start=1):
    assert spectra.unit == 'K'

    sp = spectra
    sp.plotter.plotkwargs = {}
    sp.xarr.convert_to_unit(u.GHz)
    # assume baselined already
    # sp.data -= np.nanpercentile(sp.data, 10)
    #sp.data *= beam.jtok(sp.xarr)
    #sp.unit = u.K

    linenamecen = [(x[0],float(x[1].strip('GHz')))
                   for x in line_to_image_list.line_to_image_list
                   if 'CH3OH' in x[0][:5]]

    freqs = u.Quantity([nu for ln,nu in linenamecen], u.GHz)
    # three checks:
    # 1) is the line in range?
    # 2) does the line correspond to finite data?
    # 3) is the closest spectral pixel actually near the line?
    #    (this is to check whether the line falls in SPW gaps)
    okfreqs = np.array([sp.xarr.in_range(nu) and
                        np.isfinite(sp.data[sp.xarr.x_to_pix(nu)]) and
                        np.min(np.abs(sp.xarr-nu)) < separation_tolerance
                        for nu in freqs], dtype='bool')

    redshift = (1-vguess_kms/constants.c.to(u.km/u.s).value)

    guesses = [x for nu in freqs[okfreqs]
               for x in (ampguess,
                         nu.value*redshift,
                         widthguess/constants.c.to(u.km/u.s).value*nu.value)]
    tied = ['','','']+[x for nu in freqs[okfreqs][1:] for x in
                       ('',
                        'p[1]+{0}'.format((nu.value-freqs[okfreqs][0].value)*redshift),
                        'p[2]')]
    #print(list(zip(guesses, tied)))
    fixed = [False,fixed_velo,fixed_width] * int((len(guesses)/3))
    limited=[(True,True)]*len(guesses)
    # velos in GHz
    limits=[(0,1000), (200,250), wrange_GHz]*int(len(guesses)/3)
    assert len(fixed) == len(guesses) == len(tied)
    sp.plotter(figure=pl.figure(fignum_start+1))
    sp.specfit(fittype='gaussian', guesses=guesses, tied=tied, fixed=fixed,
               limited=limited, limits=limits,
               annotate=False, verbose=True, renormalize=False)

    velocity_fit = -(spectra.specfit.parinfo[1].value - linenamecen[0][1])/linenamecen[0][1]*constants.c.to(u.km/u.s)

    print("{0}: v={1}".format(spectra.object_name, velocity_fit))

    sp.plotter.line_ids(line_names=[ln for ln,nu in linenamecen],
                        line_xvals=[nu*u.GHz for ln,nu in linenamecen],
                        velocity_offset=velocity_fit)

    qn_to_amp = {}

    # only one width error b/c of tied
    #width_error = sp.specfit.parinfo['WIDTH0'].error

    for (amp, freq, width) in zip(sp.specfit.parinfo[::3],
                                  sp.specfit.parinfo[1::3],
                                  sp.specfit.parinfo[2::3]):
        rfreq = freq * (1+velocity_fit/constants.c.to(u.km/u.s))
        closestind = np.argmin(np.abs(rfreq*u.GHz-slaimfreqs))
        closest = slaim[closestind]

        if amp.value < amp.error*3:
            amp.value = 0

        if abs(rfreq*u.GHz-slaimfreqs[closestind]) < separation_tolerance:
            qn_to_amp[closest['Resolved QNs']] = (amp, width, closest['Freq-GHz'],
                                                  closest['Log<sub>10</sub> (A<sub>ij</sub>)'],
                                                  closest['E_U (K)'],
                                                  closest['Upper State Degeneracy'],
                                                 )
        else:
            print("Skipped {0} for lack of fit.  Closest was {1}".format(rfreq, closest['Resolved QNs']))

    amps = u.Quantity([x[0].value for x in qn_to_amp.values()], u.K)
    amp_errs = u.Quantity([x[0].error for x in qn_to_amp.values()], u.K)
    widths = u.Quantity([x[1].value/x[2]*constants.c for x in qn_to_amp.values()])
    width_errs = u.Quantity([x[1].error/x[2]*constants.c for x in qn_to_amp.values()])
    thesefreqs = u.Quantity([x[2] for x in qn_to_amp.values()], u.GHz)
    these_aij = u.Quantity([10**x[3] for x in qn_to_amp.values()], u.s**-1)
    xaxis = u.Quantity([x[4] for x in qn_to_amp.values()], u.K)
    mydeg = [x[5] for x in qn_to_amp.values()]
    nupper = nupper_of_kkms(amps*widths*np.sqrt(2*np.pi), thesefreqs,
                            these_aij, mydeg).value
    integrated_error = (((amps*width_errs)**2+(widths*amp_errs)**2)*(2*np.pi))**0.5
    nupper_error = nupper_of_kkms(integrated_error, thesefreqs, these_aij,
                                  mydeg).value
    if plot:
        pl.figure(fignum_start).clf()

    log.info("nupper={0}, nupper_error={1}".format(nupper, nupper_error))
    result = fit_tex(xaxis, nupper, errors=nupper_error, plot=plot)
    pl.legend(loc='upper right')

    if plot:
        ntot, tex, _,_ = result
        peakamp = max(sp.specfit.parinfo.values[::3])
        width_fit = sp.specfit.parinfo['WIDTH0']*2.35 / sp.specfit.parinfo['SHIFT0'] * constants.c.to(u.km/u.s)

        log.info("vel_fit = {0},  width_fit = {1}".format(velocity_fit, width_fit))

        # critical that this is done *before* show_modelfit is called, since
        # that replaces the model
        show_gaussian_modelfit(sp,
                               vel=velocity_fit.value,
                               width=width_fit.value,
                               fignum=fignum_start+3,
                               ylim=(-0.1*peakamp,peakamp+2),
                              )

        show_modelfit(spectra, vel=velocity_fit.value, width=width_fit.value,
                      tem=tex, col=ntot,
                      fignum=fignum_start+2, ylim=(-0.1*peakamp,peakamp+2))

        # because there are 8 CH3OH lines, we can squeeze this into the bottom right
        ax = pl.subplot(3,3,9)
        result = fit_tex(xaxis, nupper, errors=nupper_error, plot=plot)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        pl.legend(loc='upper right', fontsize=10)

    return result, qn_to_amp


def show_gaussian_modelfit(spectra, vel, width, figsavename=None, fignum=1,
                           vscale=3, ylim=(-5,100), separation_tolerance=0.005*u.GHz,):
    """
    show the acutal model fits, one per panel for each of the ch3oh lines
    (this is a c&p of show_modelfit with the "sophisticated" ch3oh modeling
    excluded)

    width=fwhm

    vel, width are just to set the window limits; the actual plotted parameters
    come from the spectrum's parinfo
    """

    assert spectra.unit == 'K'

    pl.figure(fignum).clf()
    spectra.xarr.convert_to_unit(u.GHz)
    linenamecen = [(x[0],float(x[1].strip('GHz')))
                   for x in line_to_image_list.line_to_image_list
                   if 'CH3OH' in x[0][:5]]

    # three checks:
    # 1) is the line in range?
    # 2) does the line correspond to finite data?
    # 3) is the closest spectral pixel actually near the line?
    #    (this is to check whether the line falls in SPW gaps)
    okfreqs = np.array([spectra.xarr.in_range(nu) and
                        np.isfinite(spectra.data[spectra.xarr.x_to_pix(nu)]) and
                        np.min(np.abs(spectra.xarr-nu*u.GHz)) < separation_tolerance
                        for ln,nu in linenamecen], dtype='bool')


    nx,ny = 3,3
    plotnum = 1
    for ii,((linename,linecen),isOK) in enumerate(zip(linenamecen,okfreqs)):
        if not isOK:
            log.info("Skipped {0}:{1} because it is not in-band and finite"
                     .format(linename,linecen))
            continue
        ax = pl.subplot(nx,ny,plotnum)
        spectra.xarr.convert_to_unit(u.GHz)
        #spectra.xarr.convert_to_unit(u.km/u.s, refX=linecen*u.GHz)
        fcen = (1-vel/constants.c.to(u.km/u.s).value)*linecen
        dx = width/constants.c.to(u.km/u.s).value * linecen*3
        try:
            spectra.plotter(xmin=fcen-dx, xmax=fcen+dx,
                            axis=ax)
        except Exception as ex:
            if "Infinite recursion" in str(ex):
                ax.clear()
                ax.set_ylabel("")
                ax.get_yaxis().set_ticklabels([])
                ax.set_xlabel("")
                ax.get_xaxis().set_ticklabels([])
                print("Skipped {0}:{1} because it failed".format(linename,linecen))
                continue
            else:
                raise ex
        log.debug("%(module) parinfo for show_gaussian_modelfit: {0}"
                  .format(spectra.specfit.parinfo))
        spectra.specfit.plot_fit(annotate=False)
        ax.set_ylim(*ylim)
        ax.annotate(line_to_image_list.labeldict[linename],
                    (0.05, 0.85), horizontalalignment='left',
                    xycoords='axes fraction')

        if ((plotnum-1) % ny == 0) and (((plotnum-1) // nx) == 1):
            pl.ylabel("Brightness Temperature $T_B$ [K]")
            if (plotnum-1) != (ny*(nx-1)):
                ticks = pl.gca().get_yaxis().get_ticklocs()
                pl.gca().get_yaxis().set_ticks(ticks[1:])
        else:
            pl.ylabel("")
            pl.gca().get_yaxis().set_ticklabels([])
        if (plotnum-1) >= (ny*(nx-1)) and (((plotnum-1) % nx) == 1):
            pl.xlabel("Freq. [GHz]")
            #tl = pl.gca().get_yaxis().get_ticklabels()
            xax = pl.gca().get_xaxis()
            if (plotnum-1) == (nx*ny-1):
                pass
                #xax.set_ticks((0,200,400,600,800))
            else:
                pass
                #xax.set_ticks((0,200,400,600))
            xax.set_tick_params(labelsize=14)
            log.debug("Xlabel -> labeled: {0}".format(plotnum))
        else:
            log.debug("Xlabel -> blank: {0}".format(plotnum))
            pl.xlabel("")
            pl.gca().get_xaxis().set_ticklabels([])
        pl.subplots_adjust(hspace=0, wspace=0)
        plotnum += 1

    if figsavename is not None:
        pl.savefig(figsavename, bbox_inches='tight')



if __name__ == "__main__":

    import pyregion
    import paths
    import pyspeckit
    import glob
    from ch3oh_model import load_and_convert_spectra

    core_regions = pyregion.open(paths.rpath('cores.reg'))

    core_names = [r.attr[1]['text'] for r in core_regions]

    full_results = {}

    """
    NOTE: the J=25 lines are in the line wings of SO.  That can artificially
    inflate the fitted amplitude if a spw-wide continuum subtraction is done.
    """

    # now fit 'em all...
    #for cn in core_names:
    for cn in ['ALMAmm51',]: # debug
        spectra = load_and_convert_spectra('{corename}_spw[0-4]_mean.fits'
                                           .format(corename=cn))
        spectra.object_name = cn
        result, qn_to_amp = fit_each_line_with_a_gaussian(spectra, 20, 62.5, 6, plot=True)
        full_results[cn] = (result, qn_to_amp)
        pl.figure(4).savefig(paths.fpath('spectral_fits_ch3oh/ch3oh_spectra_gaussian_modelfits_{0}.png'
                                         .format(cn)))
        pl.figure(3).savefig(paths.fpath('spectral_fits_ch3oh/ch3oh_spectra_modelfits_{0}.png'
                                         .format(cn)))
        pl.figure(1).savefig(paths.fpath('spectral_fits_ch3oh/ch3oh_rotdiagram_{0}.png'
                                         .format(cn)))
