"""
copied from Sgr B2
"""
import numpy as np
import pyspeckit
from pyspeckit.spectrum.models import model
from pyspeckit.spectrum.models import lte_molecule
from astropy import units as u
from astropy import constants
from astropy.io import fits
from astropy import log
from astroquery.splatalogue import Splatalogue

import glob
import paths
from astropy import wcs
import radio_beam


from vamdclib import nodes
from vamdclib import request
from vamdclib import specmodel

import line_to_image_list

import pylab as pl

filetype = 'pdf'

pl.style.use('classic')
pl.matplotlib.rc_file('pubfiguresrc')

if filetype == 'png':
    fontsize = 8
    pl.rcParams['font.size'] = 10
    pl.rcParams['axes.labelsize'] = 12
    pl.rcParams['axes.titlesize'] = 14
    major_labelsize = 10
    minor_labelsize = 8
elif filetype == 'pdf':
    fontsize = 14
    pl.rcParams['font.size'] = 18
    pl.rcParams['axes.labelsize'] = 20
    pl.rcParams['axes.titlesize'] = 22
    major_labelsize = 20
    minor_labelsize = 16

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


def show_modelfit(spectra, vel, width, tem, col, figsavename=None, fignum=1,
                  vscale=3, ylim=(-5,100), separation_tolerance=0.005*u.GHz,):
    """
    width=fwhm
    """

    assert spectra.unit == 'K'

    if 'ch3oh' in spectra.specfit.Registry.multifitters:
        del spectra.specfit.Registry.multifitters['ch3oh']
    spectra.specfit.Registry.add_fitter('ch3oh', ch3oh_fitter(), 4)
    spectra.specfit.fitter = spectra.specfit.Registry.multifitters['ch3oh']

    #spectra.specfit.plot_model([vel, width/2.35, tem, col])

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
        spectra.xarr.convert_to_unit(u.km/u.s, refX=linecen*u.GHz)
        #fcen = (1-vel/constants.c.to(u.km/u.s).value)*linecen
        #dx = width/constants.c.to(u.km/u.s).value * linecen*3
        try:
            spectra.plotter(xmin=vel-vscale*width, xmax=vel+vscale*width,
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
        if tem > 0: # check to avoid div-by-zero error
            spectra.specfit.plot_model([vel, width/2.35, tem, col])
        else:
            print("Model had T=0 for {0}".format(linename))
        ax.set_ylim(*ylim)
        ax.annotate(line_to_image_list.labeldict[linename],
                    (0.05, 0.85), horizontalalignment='left',
                    fontsize=fontsize,
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
            #tl = pl.gca().get_yaxis().get_ticklabels()
            xax = pl.gca().get_xaxis()
            xax.set_ticks([40,45,50,55,60,65,70])
            pl.gca().tick_params(axis='both', which='major', labelsize=major_labelsize)
            pl.gca().tick_params(axis='both', which='minor', labelsize=minor_labelsize)
            if (plotnum-1) == (nx*ny-1):
                pass
                #xax.set_ticks((0,200,400,600,800))
            else:
                pass
                #xax.set_ticks((0,200,400,600))
            xax.set_tick_params(labelsize=major_labelsize)
            log.debug("Xlabel -> labeled: {0}".format(plotnum))
        else:
            log.debug("Xlabel -> blank: {0}".format(plotnum))
            pl.xlabel("")
            pl.gca().get_xaxis().set_ticklabels([])
        pl.subplots_adjust(hspace=0, wspace=0)
        pl.gca().tick_params(axis='both', which='major', labelsize=major_labelsize)
        pl.gca().tick_params(axis='both', which='minor', labelsize=minor_labelsize)
        plotnum += 1

    if figsavename is not None:
        pl.savefig(figsavename, bbox_inches='tight')


def load_and_convert_spectra(globname):
    speclist = [pyspeckit.Spectrum(fn) for fn in
                glob.glob(paths.spath(globname))]
    for sp in speclist:
        sp.data -= np.nanpercentile(sp.data, 25)
    spectra = pyspeckit.Spectra(speclist)
    beam = radio_beam.Beam.from_fits_header(spectra.header)
    # "baseline"
    #spectra.data -= np.nanpercentile(spectra.data, 10)
    spectra.data *= beam.jtok(spectra.xarr)
    spectra.unit='K'

    spectra.xarr.refX = 220*u.GHz # hackalack

    return spectra



if __name__ == "__main__":

    import glob
    import pyspeckit
    from astropy import wcs
    from astropy import coordinates
    import pyregion
    import paths
    import radio_beam

    target = 'e2e'

    spectra = pyspeckit.Spectra(glob.glob(paths.spath("*{0}_spw[0-9]_mean.fits".format(target))))
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

    pixel_regions = pyregion.open(paths.rpath('three_e2_pixels.reg'))
    colmap = fits.open(paths.dpath('12m/moments/CH3OH_e2_cutout_columnmap.fits'))
    temmap = fits.open(paths.dpath('12m/moments/CH3OH_e2_cutout_temperaturemap.fits'))
    ch3ohwcs = wcs.WCS(colmap[0].header)

    target = 'SelectedPixel{0}'
    for selreg, tem, col, vel, width, ylim in (
        (1, 9999401, 99995.1e18, 55.2, 5.3, (-5,100)),
        (2, 9999220, 99991.3e18, 55.6, 5.3, (-5,100)),
        (3, 9999167, 99995.7e17, 53.0, 5.3, (-5,100)),
        (4, 9999127, 99991.7e17, 54.0, 5.3, (-5,60)),
        (5, 9999999, 9999999999, 55.8, 7.0, (-5,110)),
       ):

        # tem, col should be extractd from RTD fits
        reg = [rr for rr in pixel_regions if int(rr.attr[1]['text'][-1]) == selreg][0]
        coord = coordinates.SkyCoord(reg.coord_list[0], reg.coord_list[1],
                                     frame='fk5', unit=(u.deg, u.deg))
        xpix, ypix = ch3ohwcs.celestial.wcs_world2pix(coord.ra.deg,
                                                      coord.dec.deg,
                                                      0)
        tem = temmap[0].data[int(ypix), int(xpix)]
        col = colmap[0].data[int(ypix), int(xpix)]

        speclist = [pyspeckit.Spectrum(fn) for fn in
                    glob.glob(paths.spath("*{0}_spw[0-9]_mean.fits".format(target.format(selreg))))]
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

        show_modelfit(spectra, vel, width, tem, col,
                      figsavename=paths.fpath("chemistry/ch3oh_rotdiagram_fits_SelectedPixel{0}.{1}"
                                              .format(selreg, filetype)),
                      fignum=4+selreg, ylim=ylim)

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
