"""
copied from Sgr B2
"""
import numpy as np
import pyspeckit
from pyspeckit.spectrum.models import model
from pyspeckit.spectrum.models import lte_molecule
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
from astropy.io import fits
from astroquery.splatalogue import Splatalogue

from vamdclib import nodes
from vamdclib import request as r
from vamdclib import specmodel as m

tbl = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name='CH3CN',
                              energy_max=1840, energy_type='eu_k')
freqs = np.unique(tbl['Freq-GHz'])
vdiff = (np.array((freqs-freqs[0])/freqs[0])*constants.c).to(u.km/u.s)
slaim = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name='CH3CN',
                                energy_max=1840, energy_type='eu_k',
                                line_lists=['SLAIM'],
                                show_upper_degeneracy=True)
freqs = np.array(slaim['Freq-GHz'])*u.GHz
aij = slaim['Log<sub>10</sub> (A<sub>ij</sub>)']
deg = slaim['Upper State Degeneracy']
EU = (np.array(slaim['E_U (K)'])*u.K*constants.k_B).to(u.erg).value
ref_freq = 220.74726*u.GHz
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
         and x.MolecularWeight=='41'][0]

ch3cn_inchikey = ch3cn.InChIKey

# query everything for ch3cn
query_string = "SELECT ALL WHERE VAMDCSpeciesID='%s'" % ch3cn.VAMDCSpeciesID
request.setquery(query_string)
result = request.dorequest()



def ch3cn_model(xarr, vcen, width, tex, column, background=None, tbg=2.73):

    if hasattr(tex,'unit'):
        tex = tex.value
    if hasattr(tbg,'unit'):
        tbg = tbg.value
    if hasattr(column, 'unit'):
        column = column.value
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


if __name__ == "__main__":
    cube = SpectralCube.read('../FITS/longbaseline/W51e2e_CH3CN_cutout.fits')
    sp_ = cube[:,43,43]
    hdr = sp_.header
    hdr['BUNIT'] = 'K'
    sphdu = fits.PrimaryHDU(data=sp_.to(u.K, u.brightness_temperature(cube.beam,
                                                                      cube.wcs.wcs.restfrq*u.Hz)).value,
                            header=hdr)
    sp = pyspeckit.Spectrum.from_hdu(sphdu)
    sp.data -= np.percentile(sp.data, 25)
    sp.plotter()
    sp.specfit(fittype='ch3cn', guesses=[62, 2, 500, 1e16])

    # do it again, excluding k=0 and k=1
    sp2 = sp[:400]
    sp2.plotter()
    sp2.specfit(fittype='ch3cn', guesses=[62, 2, 500, 1e16])
