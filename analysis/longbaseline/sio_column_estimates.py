"""
copied from ch3cn
"""
import numpy as np
import pyspeckit
import paths
from astropy.utils.console import ProgressBar
from pyspeckit.spectrum.models import lte_molecule
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
from astropy import log
from astropy.io import fits
from astroquery.splatalogue import Splatalogue
from astropy import modeling
from astropy.convolution import Gaussian1DKernel
from astropy import wcs
import pyregion
from astropy import coordinates
import image_tools
import os


from vamdclib import nodes
from vamdclib import request as r
from vamdclib import specmodel as m
from vamdclib import specmodel

tbl = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name='SiO',
                              energy_max=1840, energy_type='eu_k')
freqs = np.unique(tbl['Freq-GHz'])
vdiff = (np.array((freqs-freqs[0])/freqs[0])*constants.c).to(u.km/u.s)
slaim = Splatalogue.query_lines(210*u.GHz, 235*u.GHz, chemical_name='SiO',
                                energy_max=1840, energy_type='eu_k',
                                line_lists=['SLAIM'],
                                noHFS=True, # there seems to be a problem where HFS
                                # for K >= 6 is included *incorrectly*
                                show_upper_degeneracy=True)
freqs = np.array(slaim['Freq-GHz'])*u.GHz
aij = slaim['Log<sub>10</sub> (A<sub>ij</sub>)']
deg = slaim['Upper State Degeneracy']
EU = (np.array(slaim['E_U (K)'])*u.K*constants.k_B).to(u.erg).value
ref_freq = 217.10498*u.GHz
vdiff = (np.array(-(freqs-ref_freq)/ref_freq)*constants.c).to(u.km/u.s).value



nl = nodes.Nodelist()
nl.findnode('cdms')
cdms = nl.findnode('cdms')

request = r.Request(node=cdms)


# Retrieve all species from CDMS
result = request.getspecies()
molecules = result.data['Molecules']

sio = [x for x in molecules.values()
       if (x.StoichiometricFormula)==('OSi')
       and (x.OrdinaryStructuralFormula == 'SiO')
      ][0]

sio_inchikey = sio.InChIKey

# query everything for sio
query_string = "SELECT ALL WHERE VAMDCSpeciesID='%s'" % sio.VAMDCSpeciesID
request.setquery(query_string)
result = request.dorequest()
vamdc_result = result



def sio_model(xarr, vcen, width, tex, column, background=None, tbg=2.73):

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

    if background is not None:
        tbg = background

    # assume equal-width channels
    kwargs = dict(rest=ref_freq)
    equiv = u.doppler_radio(**kwargs)
    channelwidth = np.abs(xarr[1].to(u.Hz, equiv) - xarr[0].to(u.Hz, equiv)).value
    velo = xarr.to(u.km/u.s, equiv).value
    mol_model = np.zeros_like(xarr).value

    freqs_ = freqs.to(u.Hz).value

    Q = m.calculate_partitionfunction(result.data['States'],
                                      temperature=tex)[sio.Id]

    jnu_bg = lte_molecule.Jnu_cgs(xarr.to(u.Hz).value, tbg)
    bg_model = np.ones_like(xarr).value * jnu_bg

    for voff, A, g, nu, eu in zip(vdiff, aij, deg, freqs_, EU):
        tau_per_dnu = lte_molecule.line_tau_cgs(tex,
                                                column,
                                                Q,
                                                g,
                                                nu,
                                                eu,
                                                10**A)
        s = np.exp(-(velo-vcen-voff)**2/(2*width**2))*tau_per_dnu/channelwidth
        jnu_mol = lte_molecule.Jnu_cgs(nu, tex)

        # the "emission" model is generally zero everywhere, so we can just
        # add to it as we move along
        mol_model = mol_model + jnu_mol*(1-np.exp(-s))

        # background is assumed to start as a uniform value, then each
        # absorption line multiplies to reduce it.  s is zero for most velocities,
        # so this is mostly bg_model *= 1
        bg_model *= np.exp(-s)

    if background:
        # subtract jnu_bg because we *must* rezero for the case of
        # having multiple components, otherwise the backgrounds add,
        # which is nonsense
        model = bg_model + mol_model - jnu_bg
    else:
        model = mol_model

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


def ntot_of_nupper(nupper, eupper, tex, degen=1):

    partition_func = specmodel.calculate_partitionfunction(vamdc_result.data['States'],
                                                           temperature=tex.value)[sio.Id]
    Q_rot = partition_func

    Ntot = nupper * (Q_rot/degen) * np.exp(eupper / (constants.k_B*tex))

    return Ntot

if __name__ == "__main__":

    # try estimating some SiO column densities

    # peak T_B ~ 250 K ~ 13 mJy/beam
    integrated_intensity = 14e3 * u.K * u.km/u.s
    nsio = ntot_of_nupper(nupper_of_kkms(integrated_intensity, ref_freq,
                                         10**aij.mean(), 1),
                          tbl[-1]['E_U (K)']*u.K*constants.k_B,
                          250*u.K)

    print("nsio = {0}, log={1}".format(nsio, np.log10(nsio.value)))

    cosmic_si_abundance = 650 / 739000

    print("nh2, assuming all Si in SiO: {0}".format(nsio/cosmic_si_abundance))

    # Leurini+ 2013 suggest X_sio ~ 1-5x10^-8, upper limit 2e-7
    xsio = 1e-7

    print("nh2, X(SiO) = {1}: {0}".format(nsio/xsio, xsio))

    beam_area = 1e31*u.cm**2

    print("Mass assuming 1 beam = {0} to {1}"
          .format((nsio/cosmic_si_abundance*beam_area*u.Da/0.739).to(u.M_sun),
                  (nsio/xsio*beam_area*2.8*u.Da).to(u.M_sun)))

    peak_velocity = 115*u.km/u.s
    # assume ~2 beams...
    distance = (0.07 * u.arcsec * 5.4*u.kpc).to(u.km, u.dimensionless_angles())

    timescale = distance / peak_velocity

    print("Mass rate = {0} to {1}"
          .format((nsio/cosmic_si_abundance*beam_area*u.Da/0.739).to(u.M_sun)/timescale.to(u.yr),
                  (nsio/xsio*beam_area*2.8*u.Da).to(u.M_sun)/timescale.to(u.yr)))

    north_ds = paths.dpath('longbaseline/W51north_siocube_downsampled.fits')
    if os.path.exists(north_ds):
        sm_sio_cube = SpectralCube.read(north_ds)
    else:
        siocube = (SpectralCube.read(
            paths.dpath('longbaseline/linked/W51northcax.SPW0_ALL_medsub_cutout.fits'))
            .with_spectral_unit(u.km/u.s, rest_value=ref_freq,
                                velocity_convention='radio')
            .spectral_slab(-140*u.km/u.s, 260*u.km/u.s)
        )
        fwhm_factor = np.sqrt(8*np.log(2))
        hanning_factor = 1129/977
        current_resolution = np.mean(np.diff(siocube.spectral_axis)) * hanning_factor
        target_resolution = 10.0 * u.km/u.s
        pixel_scale = current_resolution
        gaussian_width = ((target_resolution**2 - current_resolution**2)**0.5 /
                          pixel_scale / fwhm_factor)
        kernel = Gaussian1DKernel(gaussian_width)

        new_xaxis = np.arange(-140, 265, 5) * u.km/u.s
        sm_sio_cube = siocube.spectral_smooth(kernel).spectral_interpolate(new_xaxis)

        sm_sio_cube.write(north_ds)


    celhdr = sm_sio_cube.wcs.celestial.to_header()
    celhdr['NAXIS1'] = sm_sio_cube.shape[2]
    celhdr['NAXIS2'] = sm_sio_cube.shape[1]
    bluemask = pyregion.open(paths.rpath('sio_blue_boxmask_lb_north.reg')).as_imagecoord(celhdr).get_mask(sm_sio_cube[0,:,:].hdu)
    redmask = pyregion.open(paths.rpath('sio_red_boxmask_lb_north.reg')).as_imagecoord(celhdr).get_mask(sm_sio_cube[0,:,:].hdu)

    sio_m0_blue = sm_sio_cube.with_mask(sm_sio_cube > 3*u.mJy).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).moment0() / u.Jy * sm_sio_cube.beam.jtok(ref_freq)
    sio_m0_red = sm_sio_cube.with_mask(sm_sio_cube > 3*u.mJy).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).moment0() / u.Jy * sm_sio_cube.beam.jtok(ref_freq)
    sio_m1_blue = sm_sio_cube.with_mask(sm_sio_cube > 3*u.mJy).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).moment1()
    sio_m1_red = sm_sio_cube.with_mask(sm_sio_cube > 3*u.mJy).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).moment1()

    north_center = coordinates.SkyCoord(290.91688055555557*u.deg,
                                        14.51818888888889*u.deg,
                                        frame='icrs')
    nc_x, nc_y = sm_sio_cube.wcs.celestial.wcs_world2pix(north_center.ra.deg,
                                                         north_center.dec.deg,
                                                         0)
    yy,xx = np.indices(sio_m0_blue.shape)
    rr = ((yy-nc_y)**2 + (xx-nc_x)**2)**0.5

    nr,centers,profile = image_tools.azimuthalAverage(np.nan_to_num(sio_m0_blue),
                                                      mask=np.isfinite(sio_m0_blue),
                                                      center=(nc_x,nc_y),
                                                      binsize=2,
                                                      return_nr=True)

    vnr,vcenters,vprofile = image_tools.azimuthalAverage(np.nan_to_num(sio_m1_blue),
                                                         mask=np.isfinite(sio_m1_blue),
                                                         center=(nc_x,nc_y),
                                                         binsize=2,
                                                         return_nr=True)

    pixscale = (wcs.utils.proj_plane_pixel_scales(sm_sio_cube.wcs)[0]*u.deg * 5.4*u.kpc).to(u.pc, u.dimensionless_angles())
    inclination = 45*u.deg
    # hard-coded: 50 pixels was measured approximately from the radial graph
    # 42 km/s was also selected as v_lsr = 18, v_lsr(north)=60
    age = (pixscale*50 / (42*u.km/u.s / np.tan(inclination))).to(u.yr)

    ppbeam = (beam_area / pixscale**2).decompose()

    nsio_profile = ntot_of_nupper(nupper_of_kkms((profile*u.K*u.km/u.s),
                                                 ref_freq, 10**aij.mean(), 1),
                                  tbl[-1]['E_U (K)']*u.K*constants.k_B,
                                  250*u.K)

    m_sio_profile = (nsio_profile * nr * beam_area / xsio * 2.8*u.Da).to(u.M_sun) / ppbeam

    massloss_rate = m_sio_profile / age

    single_event_rate = np.nansum(massloss_rate)
