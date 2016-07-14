"""

"""
import numpy as np
import pyspeckit
from pyspeckit.spectrum.models import model
from pyspeckit.spectrum.models import lte_molecule
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
from astropy.io import fits
from astroquery.splatalogue import Splatalogue, utils
from astropy import modeling

try:
    from .rotational_diagram_maps import nupper_of_kkms
except SystemError:
    from rotational_diagram_maps import nupper_of_kkms

from vamdclib import nodes
from vamdclib import request
from vamdclib import specmodel

tbl = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name=' CH3OCHO ',
                              energy_max=2500, energy_type='eu_k')
freqs = np.unique(tbl['Freq-GHz'])
#vdiff = (np.array((freqs-freqs[0])/freqs[0])*constants.c).to(u.km/u.s)

# SLAIM is missing an important (10,5) transition
slaim = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name=' CH3OCHO ',
                                energy_max=2500, energy_type='eu_k',
                                line_lists=['SLAIM'],
                                show_upper_degeneracy=True)
freqs = np.array(slaim['Freq-GHz'])*u.GHz
aij = slaim['Log<sub>10</sub> (A<sub>ij</sub>)']
deg = slaim['Upper State Degeneracy']
EU = (np.array(slaim['E_U (K)'])*u.K*constants.k_B).to(u.erg).value

# CDMS doesn't have CH3OCHO!!! It has CH2(OH)CHO though
# cdmsplat_ = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name='Glycolaldehyde',
#                                     energy_max=2500, energy_type='eu_k',
#                                     line_lists=['CDMS'],
#                                     show_upper_degeneracy=True)
# cdmsplat = utils.minimize_table(cdmsplat_)
# freqs = np.array(cdmsplat['Freq'])*u.GHz
# aij = cdmsplat['log10_Aij']
# deg = cdmsplat_['Upper State Degeneracy']
# EU = (np.array(cdmsplat['EU_K'])*u.K*constants.k_B).to(u.erg).value
#ref_freq = 220.74726*u.GHz
#vdiff = (np.array(-(freqs-ref_freq)/ref_freq)*constants.c).to(u.km/u.s).value



# nl = nodes.Nodelist()
# nl.findnode('cdms')
# cdms = nl.findnode('cdms')
# 
# request = request.Request(node=cdms)
# 
# 
# # Retrieve all species from CDMS
# result = request.getspecies()
# molecules = result.data['Molecules']
# 
# ch3ocho = [x for x in molecules.values()
#          if #hasattr(x,'MolecularWeight') and
#          #(x.ChemicalName == 'Isocyanic acid') and
#          #(x.StoichiometricFormula)==('CHNO') and
#          (x.OrdinaryStructuralFormula=='CH3OCHO')
#          #x.MolecularWeight=='32'
#         ][0]
# 
# ch3ocho_inchikey = ch3ocho.InChIKey
# 
# # query everything for ch3ocho
# query_string = "SELECT ALL WHERE VAMDCSpeciesID='%s'" % ch3ocho.VAMDCSpeciesID
# request.setquery(query_string)
# result = request.dorequest()



def ch3ocho_model(xarr, vcen, width, tex, column, background=None, tbg=2.73):

    if hasattr(tex,'unit'):
        tex = tex.value
    if hasattr(tbg,'unit'):
        tbg = tbg.value
    if hasattr(column, 'unit'):
        column = column.value
    if column < 25:
        column = 10**column
    if hasattr(vcen, 'unit'):
        vcen = vcen.value
    if hasattr(width, 'unit'):
        width = width.value

    ckms = constants.c.to(u.km/u.s).value

    # assume equal-width channels
    #kwargs = dict(rest=ref_freq)
    #equiv = u.doppler_radio(**kwargs)
    channelwidth = np.abs(xarr[1].to(u.Hz, ) - xarr[0].to(u.Hz, )).value
    #velo = xarr.to(u.km/u.s, equiv).value
    freq = xarr.to(u.Hz).value # same unit as nu below
    model = np.zeros_like(xarr).value

    freqs_ = freqs.to(u.Hz).value

    #Q = specmodel.calculate_partitionfunction(result.data['States'],
    #                                          temperature=tex)[ch3ocho.Id]

    Q = (deg * np.exp(-EU*u.erg / (constants.k_B * tex*u.K))).sum()

    for A, g, nu, eu in zip(aij, deg, freqs_, EU):
        taudnu = lte_molecule.line_tau_cgs(tex,
                                           column,
                                           Q,
                                           g,
                                           nu,
                                           eu,
                                           10**A)
        width_dnu = width / ckms * nu
        effective_linewidth_dnu = (2 * np.pi)**0.5 * width_dnu
        fcen = (1 - vcen/ckms) * nu
        tauspec = (np.exp(-(freq - fcen)**2 / (2 * width_dnu**2)) *
                   taudnu/effective_linewidth_dnu)
        jnu = (lte_molecule.Jnu_cgs(nu, tex)-lte_molecule.Jnu_cgs(nu, tbg))

        model = model + jnu*(1-np.exp(-tauspec))

    if background is not None:
        return background-model
    return model


#def fit_tex(eupper, nupperoverg, verbose=False, plot=False):
#    """
#    Fit the Boltzmann diagram
#    """
#    model = modeling.models.Linear1D()
#    #fitter = modeling.fitting.LevMarLSQFitter()
#    fitter = modeling.fitting.LinearLSQFitter()
#    ok_to_fit = nupperoverg > 0
#    result = fitter(model, eupper[ok_to_fit], np.log(nupperoverg[ok_to_fit]))
#    tex = -1./result.slope*u.K
#
#    Q_rot = (deg * np.exp(-EU*u.erg / (constants.k_B * tex))).sum()
#
#    Ntot = np.exp(result.intercept + np.log(Q_rot)) * u.cm**-2
#
#    if verbose:
#        print(("Tex={0}, Ntot={1}, Q_rot={2}".format(tex, Ntot, Q_rot)))
#
#    if plot:
#        import pylab as pl
#        L, = pl.plot(eupper, np.log10(nupperoverg), 'o')
#        xax = np.array([0, eupper.max().value])
#        line = (xax*result.slope.value +
#                result.intercept.value)
#        pl.plot(xax, np.log10(np.exp(line)), '-', color=L.get_color(),
#                label='$T={0:0.1f} \log(N)={1:0.1f}$'.format(tex, np.log10(Ntot.value)))
#
#    return Ntot, tex, result.slope, result.intercept
def fit_tex(eupper, nupperoverg, verbose=False, plot=False, uplims=None):
    """
    Fit the Boltzmann diagram
    """
    model = modeling.models.Linear1D()
    #fitter = modeling.fitting.LevMarLSQFitter()
    fitter = modeling.fitting.LinearLSQFitter()

    nupperoverg_tofit = nupperoverg.copy()

    if uplims is not None:
        upperlim_mask = nupperoverg < uplims
        if upperlim_mask.sum() > len(nupperoverg)/2.:
            # too many upper limits = bad idea to fit.
            return 0*u.cm**-2, 0*u.K, 0, 0
        nupperoverg_tofit[upperlim_mask] = uplims[upperlim_mask]

    # always ignore negatives
    good = nupperoverg_tofit > 0
    # skip any fits that have fewer than 50% good values
    if good.sum() < len(nupperoverg_tofit)/2.:
        return 0*u.cm**-2, 0*u.K, 0, 0

    result = fitter(model, eupper[good], np.log(nupperoverg_tofit[good]))
    tex = -1./result.slope*u.K

    #partition_func = specmodel.calculate_partitionfunction(ch3oh.data['States'],
    #                                                       temperature=tex.value)
    #assert len(partition_func) == 1
    #Q_rot = tuple(partition_func.values())[0]
    Q_rot = (deg * np.exp(-EU*u.erg / (constants.k_B * tex))).sum()

    Ntot = np.exp(result.intercept + np.log(Q_rot)) * u.cm**-2

    if verbose:
        print(("Tex={0}, Ntot={1}, Q_rot={2}, nuplim={3}".format(tex, Ntot, Q_rot, upperlim_mask.sum())))

    if plot:
        import pylab as pl
        L, = pl.plot(eupper, np.log10(nupperoverg_tofit), 'ro',
                     markeredgecolor='none', alpha=0.5)
        L, = pl.plot(eupper, np.log10(nupperoverg), 'o')
        xax = np.array([0, eupper.max().value])
        line = (xax*result.slope.value +
                result.intercept.value)
        pl.plot(xax, np.log10(np.exp(line)), '-', color=L.get_color(),
                label='$T={0:0.1f} \log(N)={1:0.1f}$'.format(tex, np.log10(Ntot.value)))

        if uplims is not None:
            pl.plot(eupper, np.log10(uplims), marker='_', alpha=0.5,
                    linestyle='none', color='k')

    return Ntot, tex, result.slope, result.intercept


def ch3ocho_fitter():
    """
    Generator for ch3ocho fitter class
    """

    myclass = model.SpectralModel(ch3ocho_model, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
    myclass.__name__ = "ch3ocho"
    
    return myclass

if 'ch3ocho' in pyspeckit.spectrum.fitters.default_Registry.multifitters:
    del pyspeckit.spectrum.fitters.default_Registry.multifitters['ch3ocho']
pyspeckit.spectrum.fitters.default_Registry.add_fitter('ch3ocho',ch3ocho_fitter(),4)

def ch3ocho_absorption_fitter():
    """
    Generator for ch3ocho absorption fitter class
    """

    myclass = model.SpectralModel(ch3ocho_model, 5,
            parnames=['shift','width','tex','column','background'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','T_{BG}'),
            centroid_par='shift',
            )
    myclass.__name__ = "ch3ocho_absorption"
    
    return myclass

pyspeckit.spectrum.fitters.default_Registry.add_fitter('ch3ocho_absorption',ch3ocho_absorption_fitter(),5)


if __name__ == "__main__":

    import glob
    import pyspeckit
    import paths
    import radio_beam
    import pylab as pl

    target = 'e2e'

    spectra = pyspeckit.Spectra(glob.glob(paths.spath("*{0}_spw*fits".format(target))))
    beam = radio_beam.Beam.from_fits_header(spectra.header)
    # "baseline"
    spectra.data -= np.nanpercentile(spectra.data, 10)
    spectra.data *= beam.jtok(spectra.xarr)
    spectra.unit = u.K

    spectra.xarr.convert_to_unit(u.GHz)
    spectra.plotter()

    spectra.specfit.Registry.add_fitter('ch3ocho', ch3ocho_fitter(), 4)
    spectra.specfit(fittype='ch3ocho', guesses=[55.6, 3.8, 200, 5e17],
                    limitedmin=[True]*4, limitedmax=[True]*4,
                    limits=[(50,70), (1,4), (20, 1000), (1e13, 1e18)])

    # advanced craziness: free-fit each CH3OH line
    sp = spectra[3][:1700]
    sp.xarr.convert_to_unit(u.GHz)
    sp.data -= np.nanpercentile(sp.data, 10)
    sp.data *= beam.jtok(sp.xarr)
    sp.unit = u.K
    okfreqs = np.array([sp.xarr.in_range(nu) and
                        np.isfinite(sp.data[sp.xarr.x_to_pix(nu)])
                        for nu in freqs], dtype='bool')
    okfreqs &= aij > -5
    guesses = [x for nu in freqs[okfreqs]
               for x in (20, nu.value*(1-55.626/constants.c.to(u.km/u.s).value),
                         2.79/constants.c.to(u.km/u.s).value*nu.value)]
    tied = ['','','']+[x for nu in freqs[okfreqs][1:] for x in
                       ('',
                        'p[1]+{0}'.format(nu.value-freqs[okfreqs][0].value),
                        'p[2]')]
    fixed = [False,True,True] * int((len(guesses)/3))
    limited=[(True,True)]*len(guesses)
    limits=[(0,1000), (200, 250), (0,0.05)]*int(len(guesses)/3)
    assert len(fixed) == len(guesses) == len(tied)
    sp.plotter(figure=pl.figure(3))
    sp.specfit(fittype='gaussian', guesses=guesses, tied=tied, fixed=fixed,
               limited=limited, limits=limits,
               annotate=False, verbose=True, renormalize=False)
    sp.plotter.line_ids(line_names=[str(x) for x in slaim['Resolved QNs']],
                        line_xvals=freqs, velocity_offset=55.626*u.km/u.s)

    qn_to_amp = {}

    for (amp, freq, width) in zip(sp.specfit.parinfo[::3],
                                  sp.specfit.parinfo[1::3],
                                  sp.specfit.parinfo[2::3]):
        rfreq = freq * (1+55.626/constants.c.to(u.km/u.s).value)
        closestind = np.argmin(np.abs(rfreq*u.GHz-freqs))
        closest = slaim[closestind]
        qn_to_amp[closest['Resolved QNs']] = (amp, width, closest['Freq-GHz'],
                                              closest['Log<sub>10</sub> (A<sub>ij</sub>)'],
                                              closest['E_U (K)'],
                                              closest['Upper State Degeneracy'],
                                             )

    pl.figure(4).clf()
    amps = u.Quantity([x[0].value for x in qn_to_amp.values()], u.K)
    widths = u.Quantity([x[1]/x[2]*constants.c for x in qn_to_amp.values()])
    thesefreqs = u.Quantity([x[2] for x in qn_to_amp.values()], u.GHz)
    these_aij = u.Quantity([10**x[3] for x in qn_to_amp.values()], u.s**-1)
    xaxis = u.Quantity([x[4] for x in qn_to_amp.values()], u.K)
    mydeg = [x[5] for x in qn_to_amp.values()]
    nupper = nupper_of_kkms(amps*widths*np.sqrt(2*np.pi), thesefreqs, these_aij, mydeg).value
    fit_tex(xaxis, nupper, plot=True, uplims=np.ones_like(nupper)*1e11)

    pl.savefig(paths.fpath('CH3OH_rotational_fit_notsogood.png'))
