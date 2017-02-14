"""
copied from Sgr B2

then from longbaseline/ch3cn_fits.py
"""
import paths
import numpy as np
import pyspeckit
from astropy.utils.console import ProgressBar
from pyspeckit.spectrum.models import model
from pyspeckit.spectrum.models import lte_molecule
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
from astropy import log
from astropy.io import fits
from astroquery.splatalogue import Splatalogue
from astropy import modeling
from line_to_image_list import line_to_image_list

from vamdclib import nodes
from vamdclib import request as r
from vamdclib import specmodel as m
from vamdclib import specmodel

tbl = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name='13CH3CN',
                              energy_max=1840, energy_type='eu_k')
freqs = np.unique(tbl['Freq-GHz'])
vdiff = (np.array((freqs-freqs[0])/freqs[0])*constants.c).to(u.km/u.s)
slaim = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name='13CH3CN',
                                energy_max=1840, energy_type='eu_k',
                                line_lists=['SLAIM'],
                                show_upper_degeneracy=True)
freqs = np.array(slaim['Freq-GHz'])*u.GHz
aij = slaim['Log<sub>10</sub> (A<sub>ij</sub>)']
deg = slaim['Upper State Degeneracy']
EU = (np.array(slaim['E_U (K)'])*u.K*constants.k_B).to(u.erg).value
ref_freq = 232.23419*u.GHz
vdiff = (np.array(-(freqs-ref_freq)/ref_freq)*constants.c).to(u.km/u.s).value



nl = nodes.Nodelist()
nl.findnode('cdms')
cdms = nl.findnode('cdms')

request = r.Request(node=cdms)


# Retrieve all species from CDMS
result = request.getspecies()
molecules = result.data['Molecules']

ch3cn = [x for x in molecules.values()
         if hasattr(x,'MolecularWeight') and
         (x.StoichiometricFormula)==('C2H3N')
         and x.MolecularWeight=='42'][0]

ch3cn_inchikey = ch3cn.InChIKey

# query everything for ch3cn
query_string = "SELECT ALL WHERE VAMDCSpeciesID='%s'" % ch3cn.VAMDCSpeciesID
request.setquery(query_string)
result = request.dorequest()
vamdc_result = result



def ch3cn_model(xarr, vcen, width, tex, column, background=None, tbg=2.73):

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

    # assume equal-width channels
    kwargs = dict(rest=ref_freq)
    equiv = u.doppler_radio(**kwargs)
    channelwidth = np.abs(xarr[1].to(u.Hz, equiv) - xarr[0].to(u.Hz, equiv)).value
    velo = xarr.to(u.km/u.s, equiv).value
    model = np.zeros_like(xarr).value

    freqs_ = freqs.to(u.Hz).value

    Q = m.calculate_partitionfunction(result.data['States'],
                                      temperature=tex)[ch3cn.Id]

    for voff, A, g, nu, eu in zip(vdiff, aij, deg, freqs_, EU):
        tau_per_dnu = lte_molecule.line_tau_cgs(tex,
                                                column,
                                                Q,
                                                g,
                                                nu,
                                                eu,
                                                10**A)
        s = np.exp(-(velo-vcen-voff)**2/(2*width**2))*tau_per_dnu/channelwidth
        jnu = (lte_molecule.Jnu_cgs(nu, tex)-lte_molecule.Jnu_cgs(nu, tbg))

        model = model + jnu*(1-np.exp(-s))

    if background is not None:
        return background-model
    return model

def nupper_of_kkms(kkms, freq, Aul, degeneracies):
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

def fit_tex(eupper, nupperoverg, errors=None, verbose=False, plot=False):
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
    if errors is not None:
        weights = 1/errors**2
    else:
        weights=None
    result = fitter(model, eupper[good], np.log(nupperoverg[good]),
                    weights=weights[good])
    tex = -1./result.slope*u.K

    assert not np.isnan(tex)

    partition_func = specmodel.calculate_partitionfunction(vamdc_result.data['States'],
                                                           temperature=tex.value)[ch3cn.Id]
    Q_rot = partition_func

    Ntot = np.exp(result.intercept + np.log(Q_rot)) * u.cm**-2

    if verbose:
        print(("Tex={0}, Ntot={1}, Q_rot={2}".format(tex, Ntot, Q_rot)))

    if plot:
        import pylab as pl
        if errors is not None:
            L,_,_ = pl.errorbar(x=eupper.value, y=np.log10(nupperoverg),
                                marker='o', linestyle='none',
                                yerr=np.array([np.log10(nupperoverg)-np.log10(nupperoverg-errors),
                                               np.log10(nupperoverg+errors)-np.log10(nupperoverg)]))
        else:
            L, = pl.plot(eupper, np.log10(nupperoverg), 'o')
        xax = np.array([0, eupper.max().value])
        line = (xax*result.slope.value +
                result.intercept.value)
        pl.plot(xax, np.log10(np.exp(line)), '-', color=L.get_color(),
                label='$T={0:0.1f} \log(N)={1:0.1f}$'.format(tex, np.log10(Ntot.value)))

    return Ntot, tex, result.slope, result.intercept

def fit_all_tex(xaxis, cube, cubefrequencies, degeneracies,
                einsteinAij,
                errorcube=None,
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
                                     einsteinAij, degeneracies)
            if errorcube is not None:
                enuppers = nupper_of_kkms(errorcube[:,ii,jj], cubefrequencies,
                                          einsteinAij, degeneracies)
            fit_result = fit_tex(xaxis, nuppers.value, errors=enuppers.value)
            tmap[ii,jj] = fit_result[1].value
            Nmap[ii,jj] = fit_result[0].value
        pb.update(count)
        count+=1

    return tmap,Nmap

line_name_dict = dict([(ln, u.Quantity(frq)) for ln,frq,_,_ in line_to_image_list
                       if '13CH3CN' in ln])
frequencies = sorted(u.Quantity(list(line_name_dict.values())))

line_eu = {line: EU[np.argmin(np.abs(freqs-line_name_dict[line]))]
           for line in line_name_dict if 'CH3CN' in line}
line_deg = {line: deg[np.argmin(np.abs(freqs-line_name_dict[line]))]
            for line in line_name_dict if 'CH3CN' in line}
line_aij = {line: aij[np.argmin(np.abs(freqs-line_name_dict[line]))]
            for line in line_name_dict if 'CH3CN' in line}

def multigaussian_model(xarr, vcen, width, amp1, amp2, amp3, amp4, amp5, amp6,
                        amp7, amp8, amp9,
                        background=0.0):

    xarr_frq = xarr.as_unit(u.GHz).value

    model = np.zeros(xarr.shape) + background
    for freq,amp in zip(frequencies, (amp1, amp2, amp3, amp4, amp5, amp6, amp7,
                                      amp8, amp9)):
        fcen = (freq*(1-u.Quantity(vcen, u.km/u.s) / constants.c).decompose().value).to(u.GHz).value
        fwidth = (u.Quantity(width,u.km/u.s)/constants.c).decompose().value * fcen
        model += amp * np.exp(-(xarr_frq-fcen)**2 / (2*fwidth**2))

    return model

def ch3cn_spw_fitter():
    """
    Generator for a multigaussian fitter covering the target spw
    """

    myclass = model.SpectralModel(multigaussian_model, 12,
            parnames=['shift','width', "13CH3CN 13_8", "13CH3CN 13_7",
                      "13CH3CN 13_6",  "13CH3CN 13_5", "13CH3CN 13_4",
                      "13CH3CN 13_3", "13CH3CN 13_2", "13CH3CN 13_1",
                      "13CH3CN 13_0", 'background'],
            parlimited=[(False,False),(True,False)]+[(False,False)]*9+[(True,False)],
            parlimits=[(0,0),]*12,
            shortvarnames=(r'\Delta x', r'\sigma',
                           "^{13}CH_3CN 13_8",
                           "^{13}CH_3CN 13_7", "^{13}CH_3CN 13_6",
                           "^{13}CH_3CN 13_5", "^{13}CH_3CN 13_4",
                           "^{13}CH_3CN 13_3", "^{13}CH_3CN 13_2",
                           "^{13}CH_3CN 13_1", "^{13}CH_3CN 13_0",
                           'T_{BG}'),
            centroid_par='shift',)
    myclass.__name__ = "multigauss"
    
    return myclass

pyspeckit.spectrum.fitters.default_Registry.add_fitter('13ch3cn_spw',ch3cn_spw_fitter(),11)

def ch3cn_fitter():
    """
    Generator for CH3CN fitter class
    """

    myclass = model.SpectralModel(ch3cn_model, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
    myclass.__name__ = "ch3cn"
    
    return myclass

pyspeckit.spectrum.fitters.default_Registry.add_fitter('ch3cn',ch3cn_fitter(),4)

def ch3cn_absorption_fitter():
    """
    Generator for CH3CN absorption fitter class
    """

    myclass = model.SpectralModel(ch3cn_model, 5,
            parnames=['shift','width','tex','column','background'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','T_{BG}'),
            centroid_par='shift',
            )
    myclass.__name__ = "ch3cn_absorption"
    
    return myclass

pyspeckit.spectrum.fitters.default_Registry.add_fitter('ch3cn_absorption',ch3cn_absorption_fitter(),5)


if __name__ == "__main__":

    sp = pyspeckit.Spectrum(paths.spath('e2e_radial_bin_0.00to0.38_spw2.fits'))
    sp.specfit.register_fitter('13ch3cn_spw',ch3cn_spw_fitter(),12)
    sp.xarr.convert_to_unit(u.GHz)
    sp.plotter(xmin=231.88, xmax=232.52)
    F = False
    T = True
    sp.baseline.basespec[:] = 0.380
    sp.baseline.baselinepars = [0,0.380]
    limits=[(54,59), (1.5,5),] + [(0.01,0.25)]*9 + [(0.0, 0.01)]
    fixed=[F,T,]+[T]*6+[F]*3+[T]
    sp.specfit(fittype='13ch3cn_spw', guesses=[55.7, 2.713, 0.01, 0.01, 0.01, 0.01, 0.01, 0.1, 0.2, 0.11, 0.11, 0.0], limited=[(T,T)]*12, limits=limits, fixed=fixed)
