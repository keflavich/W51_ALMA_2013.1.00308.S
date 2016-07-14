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
from vamdclib import request
from vamdclib import specmodel

tbl = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name=' CH3OH',
                              energy_max=1840, energy_type='eu_k')
freqs = np.unique(tbl['Freq-GHz'])
#vdiff = (np.array((freqs-freqs[0])/freqs[0])*constants.c).to(u.km/u.s)
slaim = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name=' CH3OH',
                                energy_max=1840, energy_type='eu_k',
                                line_lists=['SLAIM'],
                                show_upper_degeneracy=True)
freqs = np.array(slaim['Freq-GHz'])*u.GHz
aij = slaim['Log<sub>10</sub> (A<sub>ij</sub>)']
deg = slaim['Upper State Degeneracy']
EU = (np.array(slaim['E_U (K)'])*u.K*constants.k_B).to(u.erg).value
#ref_freq = 220.74726*u.GHz
#vdiff = (np.array(-(freqs-ref_freq)/ref_freq)*constants.c).to(u.km/u.s).value



nl = nodes.Nodelist()
nl.findnode('cdms')
cdms = nl.findnode('cdms')

request = request.Request(node=cdms)


# Retrieve all species from CDMS
result = request.getspecies()
molecules = result.data['Molecules']

ch3oh = [x for x in molecules.values()
         if #hasattr(x,'MolecularWeight') and
         (x.ChemicalName == 'Methanol') and
         (x.StoichiometricFormula)==('CH4O') and
         (x.OrdinaryStructuralFormula=='CH3OH')
         #x.MolecularWeight=='32'
        ][0]

ch3oh_inchikey = ch3oh.InChIKey

# query everything for ch3oh
query_string = "SELECT ALL WHERE VAMDCSpeciesID='%s'" % ch3oh.VAMDCSpeciesID
request.setquery(query_string)
result = request.dorequest()



def ch3oh_model(xarr, vcen, width, tex, column, background=None, tbg=2.73):

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
    #channelwidth = np.abs(xarr[1].to(u.Hz, ) - xarr[0].to(u.Hz, )).value
    #channelwidth = xarr.as_unit(u.Hz).dxarr
    #velo = xarr.to(u.km/u.s, equiv).value
    freq = xarr.to(u.Hz).value # same unit as nu below
    model = np.zeros_like(xarr).value

    freqs_ = freqs.to(u.Hz).value

    Q = specmodel.calculate_partitionfunction(result.data['States'],
                                              temperature=tex)[ch3oh.Id]

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

def ch3oh_fitter():
    """
    Generator for ch3oh fitter class
    """

    myclass = model.SpectralModel(ch3oh_model, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
    myclass.__name__ = "ch3oh"
    
    return myclass

pyspeckit.spectrum.fitters.default_Registry.add_fitter('ch3oh',ch3oh_fitter(),4)

def ch3oh_absorption_fitter():
    """
    Generator for ch3oh absorption fitter class
    """

    myclass = model.SpectralModel(ch3oh_model, 5,
            parnames=['shift','width','tex','column','background'],
            parlimited=[(False,False),(True,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N','T_{BG}'),
            centroid_par='shift',
            )
    myclass.__name__ = "ch3oh_absorption"
    
    return myclass

pyspeckit.spectrum.fitters.default_Registry.add_fitter('ch3oh_absorption',ch3oh_absorption_fitter(),5)


if __name__ == "__main__":

    import glob
    import pyspeckit
    import paths
    import radio_beam

    target = 'e2e'

    spectra = pyspeckit.Spectra(glob.glob(paths.spath("*{0}_spw*fits".format(target))))
    beam = radio_beam.Beam.from_fits_header(spectra.header)
    # "baseline"
    spectra.data -= np.nanpercentile(spectra.data, 10)
    spectra.data *= beam.jtok(spectra.xarr)

    spectra.plotter()

    spectra.specfit.Registry.add_fitter('ch3oh', ch3oh_fitter(), 4)
    spectra.specfit(fittype='ch3oh', guesses=[55, 4, 200, 5e15],
                    limitedmin=[True]*4, limitedmax=[True]*4, limits=[(50,70),
                                                                      (1,4),
                                                                      (20,
                                                                       1000),
                                                                      (1e13,
                                                                       1e18)])

    ok = slaim['Species'] =='CH3OHvt=0'
    spectra.xarr.refX = 225*u.GHz
    spectra.plotter.line_ids(line_names=[str(x)
                                         for x in slaim['Resolved QNs'][ok]],
                             line_xvals=freqs[ok],
                             velocity_offset=55.626*u.km/u.s)

    import pylab as pl
    pl.matplotlib.rc_file('pubfiguresrc')
    import line_to_image_list
    target = 'SelectedPixel{0}'
    for selreg, tem, col, vel, width, ylim in (
        (1, 401, 5.1e18, 55.2, 5.3, (-5,100)),
        (2, 220, 1.3e18, 55.6, 5.3, (-5,100)),
        (3, 167, 5.7e17, 53.0, 5.3, (-5,100)),
        (4, 127, 1.7e17, 54.0, 5.3, (-5,60)),
       ):

        speclist = [pyspeckit.Spectrum(fn) for fn in
                    glob.glob(paths.spath("*{0}_spw*fits".format(target.format(selreg))))]
        for sp in speclist:
            sp.data -= np.nanpercentile(sp.data, 10)
        spectra = pyspeckit.Spectra(speclist)
        beam = radio_beam.Beam.from_fits_header(spectra.header)
        # "baseline"
        #spectra.data -= np.nanpercentile(spectra.data, 10)
        spectra.data *= beam.jtok(spectra.xarr)
        spectra.unit='K'

        pl.figure(selreg, figsize=(24,16)).clf()
        spectra.plotter(figure=selreg)

        if 'ch3oh' in spectra.specfit.Registry.multifitters:
            del spectra.specfit.Registry.multifitters['ch3oh']
        spectra.specfit.Registry.add_fitter('ch3oh', ch3oh_fitter(), 4)
        spectra.specfit.fitter = spectra.specfit.Registry.multifitters['ch3oh']

        spectra.specfit.plot_model([vel, width/2.35, tem, col])

        pl.figure(4+selreg).clf()
        spectra.xarr.convert_to_unit(u.GHz)
        linenamecen = [(x[0],float(x[1].strip('GHz'))) for x in line_to_image_list.line_to_image_list
                       if 'CH3OH' in x[0][:5]]

        nx,ny = 3,3
        for ii,(linename,linecen) in enumerate(linenamecen):
            plotnum=ii+1
            ax = pl.subplot(nx,ny,ii+1)
            spectra.xarr.convert_to_unit(u.GHz)
            spectra.xarr.convert_to_unit(u.km/u.s, refX=linecen*u.GHz)
            #fcen = (1-vel/constants.c.to(u.km/u.s).value)*linecen
            #dx = width/constants.c.to(u.km/u.s).value * linecen*3
            spectra.plotter(xmin=vel-3*width, xmax=vel+3*width,
                            axis=ax)
            spectra.specfit.plot_model([vel, width/2.35, tem, col])
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
                pl.xlabel("$V_{LSR}$ [km/s]")
                tl = pl.gca().get_yaxis().get_ticklabels()
                xax = pl.gca().get_xaxis()
                if (plotnum-1) == (nx*ny-1):
                    pass
                    #xax.set_ticks((0,200,400,600,800))
                else:
                    pass
                    #xax.set_ticks((0,200,400,600))
                xax.set_tick_params(labelsize=14)
                print("Xlabel -> labeled: {0}".format(plotnum))
            else:
                print("Xlabel -> blank: {0}".format(plotnum))
                pl.xlabel("")
                pl.gca().get_xaxis().set_ticklabels([])
            pl.subplots_adjust(hspace=0, wspace=0)

        pl.savefig(paths.fpath("chemistry/ch3oh_rotdiagram_fits_SelectedPixel{0}.png"
                               .format(selreg)))

        speclist = [pyspeckit.Spectrum(fn) for fn in
                    glob.glob(paths.merge_spath("*{0}_spw*7m12m_hires*fits".format(target.format(selreg))))]
        for sp in speclist:
            sp.data -= np.nanpercentile(sp.data, 10)
        spectra = pyspeckit.Spectra(speclist)
        beam = radio_beam.Beam.from_fits_header(spectra.header)
        # "baseline"
        #spectra.data -= np.nanpercentile(spectra.data, 10)
        spectra.data *= beam.jtok(spectra.xarr)
        spectra.unit='K'

        spectra.plotter(figure=selreg+8)

        del spectra.specfit.Registry.multifitters['ch3oh']
        spectra.specfit.Registry.add_fitter('ch3oh', ch3oh_fitter(), 4)
        spectra.specfit.fitter = spectra.specfit.Registry.multifitters['ch3oh']

        spectra.specfit.plot_model([vel, width/2.35, tem, col])
