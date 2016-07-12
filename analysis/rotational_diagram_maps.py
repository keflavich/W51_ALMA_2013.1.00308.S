import paths
import itertools
from spectral_cube import SpectralCube
import numpy as np
import radio_beam
from astropy import coordinates
from astropy import units as u
from astropy import constants
from astropy.utils.console import ProgressBar
from astropy import wcs
from astropy import log
from astropy.io import fits
from astropy.nddata import Cutout2D

from pyspeckit.spectrum.models import lte_molecule

from astroquery.vamdc import Vamdc
from vamdclib import specmodel

from astropy import modeling

import re
import glob

from line_to_image_list import line_to_image_list
frequencies = u.Quantity([float(row[1].strip('GHz'))
                          for row in line_to_image_list], u.GHz)
name_to_freq = {row[0]:frq for frq, row in zip(frequencies, line_to_image_list)}
freq_to_name = {frq:row[0] for frq, row in zip(frequencies, line_to_image_list)}

linere = re.compile("W51_b6_12M.(.*).image.pbcor")

ch3oh = Vamdc.query_molecule('CH3OH')
rt = ch3oh.data['RadiativeTransitions']
frqs = u.Quantity([(float(rt[key].FrequencyValue)*u.MHz).to(u.GHz,
                                                            u.spectral())
                   for key in rt])

frqs_to_ids = {frq: key for frq,key in zip(frqs, rt)}

upperStateRefs = [rt[key].UpperStateRef for key in rt]
degeneracies = [int(ch3oh.data['States'][upperStateRef].TotalStatisticalWeight)
                for upperStateRef in upperStateRefs]
einsteinAij = u.Quantity([float(rt[key].TransitionProbabilityA) for key in rt], 1/u.s)

# http://www.astro.uni-koeln.de/cdms/catalog#equations
# the units are almost certainly wrong; I don't know how to compute line strength
# from aij =(
line_strengths_smu2 = u.Quantity([(frq**-3 * deg / 1.16395e-20 * Aij).value
                                  for frq,deg,Aij in zip(frqs, degeneracies, einsteinAij)],
                                 u.esu*u.cm)


def cutout_id_chem_map(yslice=slice(367,467), xslice=slice(114,214),
                       vrange=[51,60]*u.km/u.s, sourcename='e2',
                       filelist=glob.glob(paths.dpath('12m/cutouts/*e2e8*fits')),
                       source=None, radius=None,
                       chem_name='CH3OH',
                      ):
    assert filelist

    maps = {}
    energies = {}
    degeneracies = {}
    frequencies = {}
    indices = {}


    for ii,fn in enumerate(ProgressBar(filelist)):
        if chem_name not in fn:
            log.debug("Skipping {0} because it doesn't have {1}".format(fn, chem_name))
            continue
        if 'temperature' in fn or 'column' in fn:
            continue

        if 'moment0' in fn:
            m0 = fits.getdata(fn)
            header = fits.getheader(fn)
            cutout = Cutout2D(m0, source, 2*radius, wcs=wcs.WCS(header)).data
            m0 = cutout
            try:
                beam = radio_beam.Beam.from_fits_header(header)
                jtok = beam.jtok(header['RESTFRQ']*u.Hz).value
            except:
                jtok = 222. # approximated 0.32x0.34 at 225 GHz

            m0 = m0 * jtok
        else:
            cube = SpectralCube.read(fn)[:,yslice,xslice]
            bm = cube.beams[0]
            #jtok = bm.jtok(cube.wcs.wcs.restfrq*u.Hz)
            cube = cube.to(u.K, bm.jtok_equiv(cube.wcs.wcs.restfrq*u.Hz))

            slab = cube.spectral_slab(*vrange)
            cube.beam_threshold = 1
            #contguess = cube.spectral_slab(0*u.km/u.s, 40*u.km/u.s).percentile(50, axis=0)
            #contguess = cube.spectral_slab(70*u.km/u.s, 100*u.km/u.s).percentile(50, axis=0)
            mask = (cube.spectral_axis<40*u.km/u.s) | (cube.spectral_axis > 75*u.km/u.s)
            contguess = cube.with_mask(mask[:,None,None]).percentile(30, axis=0)
            slabsub = (slab-contguess)
            slabsub.beam_threshold = 0.15
            m0 = slabsub.moment0()
            header = m0.hdu.header

        label = linere.search(fn).groups()[0]
        frq = name_to_freq[label]

        closest_ind = np.argmin(np.abs(frqs - frq))
        closest_key = list(rt.keys())[closest_ind]
        closest_rt = rt[closest_key]
        upperstate = ch3oh.data['States'][closest_rt.UpperStateRef]
        upperen = u.Quantity(float(upperstate.StateEnergyValue),
                             unit=upperstate.StateEnergyUnit)

        maps[label] = m0
        energies[label] = upperen
        degeneracies[label] = int(upperstate.TotalStatisticalWeight)
        indices[label] = closest_ind
        frequencies[label] = frq

    # make sure the dict indices don't change order
    energy_to_key = {v:k for k,v in energies.items()}
    order = sorted(energy_to_key.keys())
    keys = [energy_to_key[k] for k in order]

    cube = np.empty((len(maps),)+maps[label].shape)
    xaxis = u.Quantity([energies[k] for k in keys])
    xaxis = xaxis.to(u.erg, u.spectral()).to(u.K, u.temperature_energy())
    for ii,key in enumerate(keys):
        # divide by degeneracy
        cube[ii,:,:] = maps[key]

    frequencies = u.Quantity([frequencies[k] for k in keys])
    indices = [indices[k] for k in keys]
    degeneracies = [degeneracies[k] for k in keys]

    assert xaxis.size == cube.shape[0]

    return xaxis,cube,maps,energies,frequencies,indices,degeneracies,header


#def Nu_thin_hightex(flux, line_strength, freq, fillingfactor=1.0, tau=None):
#    """
#    Optically-thin-ish approximation for the column density of the upper state
#    of a given line assuming T_ex >> T_bg and T_ex >> h nu
#
#    Parameters
#    ----------
#    flux : quantity
#        flux density in K km/s
#    line_strength : quantity
#        The strength of the line S_ij in 10^-18 esu cm (i.e., in debye)
#    freq : quantity
#        The frequency in Hz-equivalent
#
#
#    See eqn 29 of Mangum, rearranged to use eqn 11...
#    """
#    assert flux.unit.is_equivalent(u.K*u.km/u.s)
#    assert line_strength.unit.is_equivalent(u.esu*u.cm)
#    k = constants.k_B
#    term1 = (3*k/(8*np.pi**2 * freq * line_strength**2))
#    term5 = flux.to(u.K*u.km/u.s) / fillingfactor
#    term6 = 1 if tau is None else tau/(1-np.exp(-tau))
#    return (term1*term5*term6).to(u.cm**-2)
#
#
#def Nu_thin(flux, line_strength, tex, freq, Tbg=2.7315*u.K, fillingfactor=1.0, tau=None):
#    """
#    Derived from Eqn 30 + Eqn 80 of Mangum 2015
#    """
#    assert flux.unit.is_equivalent(u.K*u.km/u.s)
#    assert line_strength.unit.is_equivalent(u.esu*u.cm)
#    hnu = constants.h * freq
#    h = constants.h
#    kbt = constants.k_B * tex
#    term1 = (3*h/(8*np.pi**2 * line_strength**2))
#    term3 = 1. / (np.exp(hnu/kbt) - 1)
#    term4 = 1./(tex-Tbg)
#    term5 = flux.to(u.K*u.km/u.s) / fillingfactor
#    term6 = 1 if tau is None else tau/(1-np.exp(-tau))
#    return (term1*term3*term4*term5*term6).to(u.cm**-2)

def nupper_of_kkms(kkms, freq, Aul, degeneracies, Tex=100*u.K):
    """ Derived directly from pyspeckit eqns..."""
    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    kkms = u.Quantity(kkms, u.K*u.km/u.s)
    #nline = 1.95e3 * freq**2 / Aul * kkms
    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2
    # term2 = np.exp(-constants.h*freq/(constants.k_B*Tex)) -1
    # term2 -> kt / hnu
    # kelvin-hertz
    Khz = (kkms * (freq/constants.c)).to(u.K * u.MHz)
    return (nline * Khz / degeneracies).to(u.cm**-2)
    return nline.value / degeneracies *u.cm**-2 # because... something wrong.
    return nline.value * u.cm**-2

#def nupper_of_kkms(kkms, freq, Aul, degeneracies):
#    """ eqn 15.28 of Wilson 2009, with degeneracy dropped """
#    """ This is definitely wrong, which is sad. """
#    freq = u.Quantity(freq, u.GHz)
#    Aul = u.Quantity(Aul, u.Hz)
#    kkms = u.Quantity(kkms, u.K*u.km/u.s)
#    nline = 1.95e3 * freq**2 / Aul * kkms
#    return nline.value / degeneracies *u.cm**-2 # because... something wrong.
#    return nline.value * u.cm**-2
#
def fit_tex(eupper, nupperoverg, verbose=False, plot=False):
    """
    Fit the Boltzmann diagram
    """
    model = modeling.models.Linear1D()
    #fitter = modeling.fitting.LevMarLSQFitter()
    fitter = modeling.fitting.LinearLSQFitter()
    # ignore negatives
    good = nupperoverg > 0
    if good.sum() < len(nupperoverg)/2.:
        return 0*u.cm**-2, 0*u.K, 0, 0
    result = fitter(model, eupper[good], np.log(nupperoverg[good]))
    tex = -1./result.slope*u.K

    partition_func = specmodel.calculate_partitionfunction(ch3oh.data['States'],
                                                           temperature=tex.value)
    assert len(partition_func) == 1
    Q_rot = tuple(partition_func.values())[0]

    Ntot = np.exp(result.intercept + np.log(Q_rot)) * u.cm**-2

    if verbose:
        print(("Tex={0}, Ntot={1}, Q_rot={2}".format(tex, Ntot, Q_rot)))

    if plot:
        import pylab as pl
        L, = pl.plot(eupper, np.log10(nupperoverg), 'o')
        xax = np.array([0, eupper.max().value])
        line = (xax*result.slope.value +
                result.intercept.value)
        pl.plot(xax, np.log10(np.exp(line)), '-', color=L.get_color(),
                label='$T={0:0.1f} \log(N)={1:0.1f}$'.format(tex, np.log10(Ntot.value)))

    return Ntot, tex, result.slope, result.intercept

def pyspeckitfit(eupper, kkms, frequencies, degeneracies, einsteinAs,
                 verbose=False, plot=False, guess=(150,1e19)):
    """
    Fit the Boltzmann diagram (but do it right, with no approximations, just
    direct forward modeling)

    This doesn't seem to work, though: it can't reproduce its own output
    """

    bandwidth = (1*u.km/u.s/constants.c)*(frequencies)

    def model(tex, col, einsteinAs=einsteinAs, eupper=eupper):
        tex = u.Quantity(tex, u.K)
        col = u.Quantity(col, u.cm**-2)
        eupper = eupper.to(u.erg, u.temperature_energy())
        einsteinAs = u.Quantity(einsteinAs, u.Hz)

        partition_func = specmodel.calculate_partitionfunction(ch3oh.data['States'],
                                                               temperature=tex.value)
        assert len(partition_func) == 1
        Q_rot = tuple(partition_func.values())[0]
        return lte_molecule.line_brightness(tex, bandwidth, frequencies,
                                            total_column=col,
                                            partition_function=Q_rot,
                                            degeneracy=degeneracies,
                                            energy_upper=eupper,
                                            einstein_A=einsteinAs)

    def minmdl(args):
        tex, col = args
        return ((model(tex, col).value - kkms)**2).sum()

    from scipy import optimize
    result = optimize.minimize(minmdl, guess, method='Nelder-Mead')
    tex, Ntot = result.x

    if plot:
        import pylab as pl
        pl.subplot(2,1,1)
        pl.plot(eupper, kkms, 'o')
        order = np.argsort(eupper)
        pl.plot(eupper[order], model(tex, Ntot)[order])
        pl.subplot(2,1,2)
        xax = np.array([0, eupper.max().value])
        pl.plot(eupper, np.log(nupper_of_kkms(kkms, frequencies, einsteinAs, degeneracies).value), 'o')
        partition_func = specmodel.calculate_partitionfunction(ch3oh.data['States'],
                                                               temperature=tex)
        assert len(partition_func) == 1
        Q_rot = tuple(partition_func.values())[0]
        intercept = np.log(Ntot) - np.log(Q_rot)
        pl.plot(xax, np.log(xax*tex + intercept), '-',
                label='$T={0:0.1f} \log(N)={1:0.1f}$'.format(tex, np.log10(Ntot)))

    return Ntot, tex, result

def test_roundtrip(cubefrequencies=[218.44005, 234.68345, 220.07849, 234.69847, 231.28115]*u.GHz,
                   degeneracies=[9, 9, 17, 11, 21],
                   xaxis=[45.45959683,  60.92357159,  96.61387286, 122.72191958, 165.34856457]*u.K,
                   indices=[3503, 1504, 2500, 116, 3322],
                  ):

    # integrated line over 1 km/s (see dnu)
    onekms = 1*u.km/u.s / constants.c
    kkms = lte_molecule.line_brightness(tex=100*u.K,
                                        total_column=1e15*u.cm**-2,
                                        partition_function=1185,
                                        degeneracy=degeneracies,
                                        frequency=cubefrequencies,
                                        energy_upper=xaxis.to(u.erg,
                                                              u.temperature_energy()),
                                        einstein_A=einsteinAij[indices],
                                        dnu=onekms*cubefrequencies) * u.km/u.s
    col, tem, slope, intcpt = fit_tex(xaxis, nupper_of_kkms(kkms,
                                                            cubefrequencies,
                                                            einsteinAij[indices],
                                                            degeneracies).value,
                                      plot=True)
    print("temperature = {0} (input was 100)".format(tem))
    print("column = {0} (input was 1e15)".format(np.log10(col.value)))


def fit_all_tex(xaxis, cube, cubefrequencies, indices, degeneracies,
                replace_bad=False):
    """
    Parameters
    ----------
    replace_bad : bool
        Attempt to replace bad (negative) values with their upper limits?
    """

    tmap = np.empty(cube.shape[1:])
    Nmap = np.empty(cube.shape[1:])

    yy,xx = np.indices(cube.shape[1:])
    pb = ProgressBar(xx.size)
    count=0

    for ii,jj in (zip(yy.flat, xx.flat)):
        if any(np.isnan(cube[:,ii,jj])):
            tmap[ii,jj] = np.nan
        else:
            if replace_bad:
                neg = cube[:,ii,jj] <= 0
                cube[neg,ii,jj] = replace_bad
            nuppers = nupper_of_kkms(cube[:,ii,jj], cubefrequencies,
                                     einsteinAij[indices], degeneracies)
            fit_result = fit_tex(xaxis, nuppers.value)
            tmap[ii,jj] = fit_result[1].value
            Nmap[ii,jj] = fit_result[0].value
        pb.update(count)
        count+=1

    return tmap,Nmap

if __name__ == "__main__":
    import pylab as pl
    pl.matplotlib.rc_file('pubfiguresrc')

    detection_threshold_jykms = 0.065
    approximate_jytok = 221

    sources = {'e2': coordinates.SkyCoord('19:23:43.963', '+14:30:34.53',
                                          frame='fk5', unit=(u.hour, u.deg)),
               'e8': coordinates.SkyCoord('19:23:43.891', '+14:30:28.13',
                                          frame='fk5', unit=(u.hour, u.deg)),
               'ALMAmm14': coordinates.SkyCoord('19:23:38.571', '+14:30:41.80',
                                                frame='fk5', unit=(u.hour,
                                                                   u.deg)),
               'north': coordinates.SkyCoord('19:23:39.906', '+14:31:05.33',
                                             frame='fk5', unit=(u.hour,
                                                                u.deg)),
              }
    radii = {'e2':2.5*u.arcsec,
             'e8':3.0*u.arcsec,
             'ALMAmm14':2.1*u.arcsec,
             'north':4.0*u.arcsec,
            }

    # use precomputed moments
    for sourcename, region in (('e2','e2e8'), ('e8','e2e8'), ('north','north'),
                               ('ALMAmm14','ALMAmm14'),):
                                
        _ = cutout_id_chem_map(source=sources[sourcename],
                               radius=radii[sourcename],
                               sourcename=sourcename,
                               filelist=glob.glob(paths.dpath('12m/moments/*medsub_moment0.fits')),
                               chem_name='.CH3OH', # use dot to exclude 13CH3OH
                              )
        xaxis,cube,maps,energies,cubefrequencies,indices,degeneracies,header = _

        pl.figure(2, figsize=(12,12)).clf()
        sample_pos = np.linspace(0,1,7)[1:-1]
        nx = len(sample_pos)
        ny = len(sample_pos)
        for ii,(spy,spx) in enumerate(itertools.product(sample_pos,sample_pos)):
            rdx = int(spx*cube.shape[2])
            rdy = int(spy*cube.shape[1])
            plotnum = (nx*ny-(2*(ii//ny)*ny)+ii-ny)+1
            pl.subplot(nx,ny,plotnum)
            Ntot, tex, slope, intcpt = fit_tex(xaxis,
                                               nupper_of_kkms(cube[:,rdy,rdx],
                                                              cubefrequencies,
                                                              einsteinAij[indices],
                                                              degeneracies).value,
                                               plot=True)
            pl.ylim(11, 15)
            pl.xlim(0, 850)
            #pl.annotate("{0:d},{1:d}".format(rdx,rdy), (0.5, 0.85), xycoords='axes fraction',
            #            horizontalalignment='center')
            pl.annotate("T={0:d}".format(int(tex.value)),
                        (0.65, 0.85), xycoords='axes fraction',
                        horizontalalignment='left', fontsize=12)
            pl.annotate("N={0:0.1f}".format(np.log10(Ntot.value)),
                        (0.65, 0.75), xycoords='axes fraction',
                        horizontalalignment='left', fontsize=12)
            pl.annotate("{0:0.2f},{1:0.2f}".format(spx,spy),
                        (0.05, 0.05), xycoords='axes fraction',
                        horizontalalignment='left', fontsize=12)

            # show upper limits
            pl.plot(xaxis,
                    np.log10(nupper_of_kkms(detection_threshold_jykms*approximate_jytok,
                                            cubefrequencies,
                                            einsteinAij[indices],
                                            degeneracies,).value),
                    linestyle='none', marker='_', color='k',
                    markeredgewidth=2, alpha=0.5)

            if (plotnum-1) % ny == 0:
                pl.ylabel("log($N_u / g_u$)")
                if (plotnum-1) != (ny*(nx-1)):
                    ticks = pl.gca().get_yaxis().get_ticklocs()
                    pl.gca().get_yaxis().set_ticks(ticks[1:])
            else:
                pl.gca().get_yaxis().set_ticklabels([])
            if (plotnum-1) >= (ny*(nx-1)):
                pl.xlabel("$E_u$ [K]")
                tl = pl.gca().get_yaxis().get_ticklabels()
                xax = pl.gca().get_xaxis()
                if (plotnum-1) == (nx*ny-1):
                    xax.set_ticks((0,200,400,600,800))
                else:
                    xax.set_ticks((0,200,400,600))
                xax.set_tick_params(labelsize=14)
            else:
                pl.gca().get_xaxis().set_ticklabels([])
            #pl.legend(loc='best', fontsize='small')
        pl.subplots_adjust(hspace=0, wspace=0)
        pl.savefig(paths.fpath("chemistry/ch3oh_rotation_diagrams_{0}.png".format(sourcename)))

        tmap,Nmap = fit_all_tex(xaxis, cube, cubefrequencies, indices, degeneracies,
                                replace_bad=detection_threshold_jykms*approximate_jytok)

        pl.figure(1).clf()
        pl.imshow(tmap, vmin=0, vmax=600, cmap='hot')
        cb = pl.colorbar()
        cb.set_label("Temperature (K)")
        pl.savefig(paths.fpath("chemistry/ch3oh_temperature_map_{0}.png".format(sourcename)))
        pl.figure(3).clf()
        pl.imshow(np.log10(Nmap), vmin=16, vmax=19, cmap='viridis')
        cb = pl.colorbar()
        cb.set_label("log N(CH$_3$OH)")
        pl.savefig(paths.fpath("chemistry/ch3oh_column_map_{0}.png".format(sourcename)))

        hdu = fits.PrimaryHDU(data=tmap, header=header)
        hdu.writeto(paths.dpath('12m/moments/CH3OH_{0}_cutout_temperaturemap.fits'.format(sourcename)), clobber=True)

        hdu = fits.PrimaryHDU(data=Nmap, header=header)
        hdu.writeto(paths.dpath('12m/moments/CH3OH_{0}_cutout_columnmap.fits'.format(sourcename)), clobber=True)




    for sourcename, region, xslice, yslice, vrange, rdposns in (
        #('e2','e2e8',(114,214),(367,467),(51,60),[(10,10),(60,84),]),
        #('e8','e2e8',(119,239),(227,347),(52,63),[(10,60),(65,45),]),
        #('north','north',(152,350),(31,231),(54,64),[(100,80),(75,80),]),
        #('ALMAmm14','ALMAmm14',(80,180),(50,150),(58,67),[(65,40),(45,40),]),
    ):

        _ = cutout_id_chem_map(yslice=slice(*yslice), xslice=slice(*xslice),
                               vrange=vrange*u.km/u.s, sourcename=sourcename,
                               filelist=glob.glob(paths.dpath('12m/cutouts/*{0}*fits').format(region)),
                               chem_name='CH3OH',
                              )
        xaxis,cube,maps,energies,cubefrequencies,indices,degeneracies,header = _


        pl.figure(2).clf()
        for rdx,rdy in rdposns:
            fit_tex(xaxis, nupper_of_kkms(cube[:,rdy,rdx], cubefrequencies,
                                          einsteinAij[indices],
                                          degeneracies).value, plot=True)
        pl.ylabel("log($N_u / g_u$)")
        pl.xlabel("$E_u$ [K]")
        pl.legend(loc='best')
        pl.savefig(paths.fpath("chemistry/ch3oh_rotation_diagrams_{0}.png".format(sourcename)))

        tmap,Nmap = fit_all_tex(xaxis, cube, cubefrequencies, indices, degeneracies)

        pl.figure(1).clf()
        pl.imshow(tmap, vmin=0, vmax=1000, cmap='hot')
        cb = pl.colorbar()
        cb.set_label("Temperature (K)")
        pl.savefig(paths.fpath("chemistry/ch3oh_temperature_map_{0}.png".format(sourcename)))
        pl.figure(3).clf()
        pl.imshow(np.log10(Nmap), vmin=15, vmax=21, cmap='viridis')
        cb = pl.colorbar()
        cb.set_label("log N(CH$_3$OH)")
        pl.savefig(paths.fpath("chemistry/ch3oh_column_map_{0}.png".format(sourcename)))

        hdu = fits.PrimaryHDU(data=tmap, header=header)
        hdu.writeto(paths.dpath('12m/cutouts/CH3OH_{0}_cutout_temperaturemap.fits'.format(sourcename)), clobber=True)

        hdu = fits.PrimaryHDU(data=Nmap, header=header)
        hdu.writeto(paths.dpath('12m/cutouts/CH3OH_{0}_cutout_columnmap.fits'.format(sourcename)), clobber=True)
