"""
copied from ch3cn
"""
import pylab as pl
import numpy as np
import pyspeckit
import paths
from astropy.utils.console import ProgressBar
from pyspeckit.spectrum.models import lte_molecule
from pyspeckit.spectrum.models.lte_molecule import nupper_of_kkms, ntot_of_nupper
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
from astropy import log
from astropy.io import fits
from astroquery.splatalogue import Splatalogue, utils
from astropy import modeling
from astropy.convolution import Gaussian1DKernel
from astropy import wcs
from astropy.visualization import quantity_support
import pyregion
from astropy import coordinates
import image_tools
import os

import warnings

import sys

quantity_support()

# from vamdclib import nodes
# from vamdclib import request as r
# from vamdclib import specmodel as m
# from vamdclib import specmodel

try:
    with open('sio_column_estimate_results.txt', 'w') as outfile:
        sys.stdout = outfile

        dw51 = 5.4*u.kpc

        warnings.filterwarnings('ignore', message='This function (.*) requires loading the entire cube', append=True)
        warnings.filterwarnings('ignore', message='All-NaN slice encountered')
        warnings.filterwarnings('ignore', message='invalid value encountered in true_divide')
        fmin,fmax = 217*u.GHz, 218*u.GHz

        tbl = Splatalogue.query_lines(fmin, fmax, chemical_name='SiO',
                                      energy_max=1840, energy_type='eu_k')
        tbl = utils.clean_column_headings(tbl)
        tbl = tbl[tbl['Species'] == 'SiOv=0']
        freqs = np.unique(tbl['FreqGHz'])
        vdiff = (np.array((freqs-freqs[0])/freqs[0])*constants.c).to(u.km/u.s)
        slaim = Splatalogue.query_lines(fmin, fmax, chemical_name=' SiO ',
                                        energy_max=1840, energy_type='eu_k',
                                        line_lists=['SLAIM'],
                                        noHFS=True, # there seems to be a problem where HFS
                                        # for K >= 6 is included *incorrectly*
                                        show_upper_degeneracy=True)
        slaim = utils.clean_column_headings(slaim)
        freqs = np.array(slaim['FreqGHz'])*u.GHz
        aij = slaim['log10_Aij']
        deg = slaim['Upper State Degeneracy'][0]
        EU = (np.array(slaim['EU_K'])*u.K*constants.k_B).to(u.erg).value
        ref_freq = 217.10498*u.GHz
        vdiff = (np.array(-(freqs-ref_freq)/ref_freq)*constants.c).to(u.km/u.s).value

        freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters(molecule_name='SiO', fmin=fmin,
                                              fmax=fmax)
        deg = deg.item()


        # nl = nodes.Nodelist()
        # nl.findnode('cdms')
        # cdms = nl.findnode('cdms')

        # request = r.Request(node=cdms)


        # # Retrieve all species from CDMS
        # result = request.getspecies()
        # molecules = result.data['Molecules']

        # sio = [x for x in molecules.values()
        #        if (x.StoichiometricFormula)==('OSi')
        #        and (x.OrdinaryStructuralFormula == 'SiO')
        #       ][0]

        # sio_inchikey = sio.InChIKey

        # # query everything for sio
        # query_string = "SELECT ALL WHERE VAMDCSpeciesID='%s'" % sio.VAMDCSpeciesID
        # request.setquery(query_string)
        # result = request.dorequest()
        # vamdc_result = result



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

            #Q = m.calculate_partitionfunction(result.data['States'],
            #                                  temperature=tex)[sio.Id]
            Q = partfunc(tex)

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




        if __name__ == "__main__":


            # SANITY CHECK: Comparison to eqn 3 of 
            # https://ui.adsabs.harvard.edu/abs/2019ApJ...878...29L/abstract
            tex = 250*u.K
            freq = freqs.item()
            ourvalue = ntot_of_nupper(nupper_of_kkms(1*u.K*u.km/u.s, ref_freq,
                                                     10**aij.item()), tex=tex,
                                      degeneracy=11, Q_rot=partfunc(tex),
                                      eupper=EU.item()*u.erg)
            theirvalue = (1.6e11 * u.cm**-2 *
                          (tex+0.35*u.K)*np.exp(31.26*u.K/tex) /
                          (np.exp(10.4*u.K/tex)-1) / (lte_molecule.Jnu(freq,
                                                                       tex) -
                                                      lte_molecule.Jnu(freq,
                                                                       2.73*u.K))).decompose()
            print(f"Sanity check: this should be ~1 = {(ourvalue/theirvalue).decompose()}")

            # try estimating some SiO column densities

            # peak T_B ~ 250 K ~ 13 mJy/beam
            tex = 250*u.K
            # this value should come from the peak integrated intensity below
            integrated_intensity = 14e3 * u.K * u.km/u.s
            nsio = ntot_of_nupper(nupper_of_kkms(integrated_intensity, ref_freq,
                                                 10**aij.mean()),
                                  tbl[-1]['EU_K']*u.K*constants.k_B,
                                  degeneracy=deg,
                                  Q_rot=partfunc(tem=tex),
                                  tex=tex)

            print("north nsio = {0:0.3g}, log={1:0.3g}".format(nsio, np.log10(nsio.value)))

            cosmic_si_abundance = 650 / 739000

            print(f"north nh2, assuming all Si in SiO (X={cosmic_si_abundance:0.3g}): {nsio/cosmic_si_abundance:0.3g}")

            # Leurini+ 2013 suggest X_sio ~ 1-5x10^-8, upper limit 2e-7
            xsio = 1e-7

            print("north nh2, X(SiO) = {1:0.3g}: {0:0.3g}".format(nsio/xsio, xsio))


            peak_velocity = 115*u.km/u.s
            # assume ~2 beams...
            distance = (0.07 * u.arcsec * dw51).to(u.km, u.dimensionless_angles())

            timescale = distance / peak_velocity


            north_ds = paths.dpath('longbaseline/W51north_siocube_downsampled.fits')
            if os.path.exists(north_ds):
                sm_sio_cube_north = sm_sio_cube = SpectralCube.read(north_ds)
            elif os.path.exists(north_ds+".gz"):
                sm_sio_cube_north = sm_sio_cube = SpectralCube.read(north_ds+".gz")
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

            sm_sio_cube = sm_sio_cube.to(u.K, sm_sio_cube.beam.jtok_equiv(ref_freq))

            beam_area = (sm_sio_cube.beam.sr * dw51**2).to(u.cm**2, u.dimensionless_angles())
            print(f"Beam area = {beam_area:0.3g}")

            print("north Mass assuming 1 beam = {0:0.3g} to {1:0.3g}"
                  .format((nsio/cosmic_si_abundance*beam_area*u.Da/0.739).to(u.M_sun),
                          (nsio/xsio*beam_area*2.8*u.Da).to(u.M_sun)))

            print("north Mass rate = {0:0.3g} (cosmic) to {1:0.3g} (X=1e-7)"
                  .format((nsio/cosmic_si_abundance*beam_area*u.Da/0.739).to(u.M_sun)/timescale.to(u.yr),
                          (nsio/xsio*beam_area*2.8*u.Da).to(u.M_sun)/timescale.to(u.yr)))

            celhdr = sm_sio_cube.wcs.celestial.to_header()
            celhdr['NAXIS1'] = sm_sio_cube.shape[2]
            celhdr['NAXIS2'] = sm_sio_cube.shape[1]
            bluemask = pyregion.open(paths.rpath('sio_blue_boxmask_lb_north.reg')).as_imagecoord(celhdr).get_mask(sm_sio_cube[0,:,:].hdu)
            redmask = pyregion.open(paths.rpath('sio_red_boxmask_lb_north.reg')).as_imagecoord(celhdr).get_mask(sm_sio_cube[0,:,:].hdu)

            # noise is about 0.6 mJy, so 3 mJy is 5-sigma
            thrsh = (sm_sio_cube.beam.jtok(ref_freq) / u.Jy * 3*u.mJy).to(u.K)

            sio_m0_blue_north = sm_sio_cube.with_mask(sm_sio_cube > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).moment0()
            sio_m0_red_north = sm_sio_cube.with_mask(sm_sio_cube > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).moment0()
            sio_m1_blue_north = sm_sio_cube.with_mask(sm_sio_cube > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).moment1()
            sio_m1_red_north = sm_sio_cube.with_mask(sm_sio_cube > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).moment1()
            sio_peak_blue_north = sm_sio_cube.with_mask(sm_sio_cube > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).max(axis=0)
            sio_peak_red_north = sm_sio_cube.with_mask(sm_sio_cube > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).max(axis=0)

            print(f"north max integrated intensity: blue={np.nanmax(sio_m0_blue_north):0.3g}, red={np.nanmax(sio_m0_red_north):0.3g}")

            pixscale = (wcs.utils.proj_plane_pixel_scales(sm_sio_cube.wcs)[0]*u.deg * dw51).to(u.pc, u.dimensionless_angles())

            kkms_to_mass = ntot_of_nupper(nupper_of_kkms(1*u.K*u.km/u.s, ref_freq,
                                                         10**aij.mean()),
                                          eupper=tbl[-1]['EU_K']*u.K*constants.k_B,
                                          degeneracy=deg,
                                          Q_rot=partfunc(tem=tex),
                                          tex=tex) * 2.8*u.Da * pixscale**2
            dv = np.mean(np.diff(sm_sio_cube.spectral_axis))

            mass_sm_sio_cube = ((sm_sio_cube * kkms_to_mass * dv) /
                                (u.K*u.km/u.s)).to(u.M_sun)
            print(f"kkms_to_mass = {kkms_to_mass.to(u.M_sun) / (u.K*u.km/u.s)}")
            print(f"north Total SiO mass: {mass_sm_sio_cube.sum()}, h2 for xsio={xsio} -> {mass_sm_sio_cube.sum()/xsio}")
            totalnorthmass_masked = mass_sm_sio_cube.with_mask(sm_sio_cube > thrsh).sum()
            print(f"north Total SiO mass: {totalnorthmass_masked}, h2 for xsio={xsio} -> {totalnorthmass_masked/xsio}")
            velcenter = 60*u.km/u.s
            velaxis = mass_sm_sio_cube.spectral_axis - velcenter
            sio_momentum_blue_north = ((mass_sm_sio_cube * velaxis[:,None,None])
                                       .with_mask(sm_sio_cube > thrsh)
                                       .with_mask(bluemask)
                                       .spectral_slab(-140*u.km/u.s, 60*u.km/u.s)
                                       .sum(axis=0))
            sio_momentum_red_north = ((mass_sm_sio_cube * velaxis[:,None,None])
                                      .with_mask(sm_sio_cube > thrsh)
                                      .with_mask(redmask)
                                      .spectral_slab(60*u.km/u.s, 260*u.km/u.s)
                                      .sum(axis=0))

            north_center = coordinates.SkyCoord(290.91688055555557*u.deg,
                                                14.51818888888889*u.deg,
                                                frame='icrs')
            nc_x, nc_y = sm_sio_cube.wcs.celestial.wcs_world2pix(north_center.ra.deg,
                                                                 north_center.dec.deg,
                                                                 0)
            yy,xx = np.indices(sio_m0_blue_north.shape)
            rr = ((yy-nc_y)**2 + (xx-nc_x)**2)**0.5

            nr,centers,profile = image_tools.azimuthalAverage(np.nan_to_num(sio_m0_blue_north),
                                                              mask=np.isfinite(sio_m0_blue_north),
                                                              center=(nc_x,nc_y),
                                                              binsize=2,
                                                              return_nr=True)

            vnr,vcenters,vprofile = image_tools.azimuthalAverage(np.nan_to_num(sio_m1_blue_north),
                                                                 mask=np.isfinite(sio_m1_blue_north),
                                                                 center=(nc_x,nc_y),
                                                                 binsize=2,
                                                                 return_nr=True)

            # inclination is probably close to 75 deg, but we want to put it in the paper
            # as "quantity tan(ii/45)"
            inclination = 75*u.deg
            inclination = 45*u.deg
            # hard-coded: 50 pixels was measured approximately from the radial graph
            # 42 km/s was also selected as v_lsr = 18, v_lsr(north)=60
            max_velocity = (60+25)*u.km/u.s
            max_velocity = 25*u.km/u.s
            max_velocity = 95 * u.km/u.s # we see blue emission starting at -35 km/s

            age = (pixscale*50 / max_velocity).to(u.yr) / np.tan(inclination)
            dmax = (0.285*u.arcsec*dw51).to(u.pc, u.dimensionless_angles())
            age = (dmax / max_velocity).to(u.yr) / np.tan(inclination)

            print("North inclination = {0:0.3g}".format(inclination))
            print("North: dmax/velocity = {0:0.3g}".format((dmax/max_velocity).to(u.yr)))
            print("north age estimate from max velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(max_velocity, dmax.to(u.au), age))

            # second age estimate
            # this time, we use the location of the highest value in the moment 0 map
            # to represent the "average" outflowing mass.
            # If we're looking at a single event, i.e., a hubble flow, the age should
            # be the same for the average and the max.
            # If we're looking at a constantly-driven outflow... we should see the same
            # velocities at all points and the age should be lower because we're looking closer
            # to the source

            # compute the positions of brightest integrated emission
            bymax,bxmax = np.unravel_index(np.nanargmax(sio_m0_blue_north), sio_m0_blue_north.shape)
            rymax,rxmax = np.unravel_index(np.nanargmax(sio_m0_red_north), sio_m0_red_north.shape)

            # find the offset from the central position, the velocity offset, and compute age
            m0pk_sep_blue_north = (((bxmax-nc_x)**2 + (bymax-nc_y)**2)**0.5 * pixscale).to(u.au)
            vatmax_blue_north = sio_m1_blue_north[bymax, bxmax] - velcenter
            avg_age_blue_north = (m0pk_sep_blue_north / vatmax_blue_north).to(u.yr) / np.tan(inclination)

            m0pk_sep_red_north = (((rxmax-nc_x)**2 + (rymax-nc_y)**2)**0.5 * pixscale).to(u.au)
            vatmax_red_north = sio_m1_red_north[rymax, rxmax] - velcenter
            avg_age_red_north = (m0pk_sep_red_north / vatmax_red_north).to(u.yr) / np.tan(inclination)

            print("north blue age estimate from avg velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(vatmax_blue_north, m0pk_sep_blue_north, avg_age_blue_north))
            print("north red age estimate from avg velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(vatmax_red_north, m0pk_sep_red_north, avg_age_red_north))



            ppbeam = (beam_area / pixscale**2).decompose()

            print("Velocity inclination factor = 1/cos({0:0.3g}) = {1:0.3g}"
                  .format(inclination, 1/np.cos(inclination)))
            print("north total momentum (inclination-corrected): blue={0:0.3g}, red={1:0.3g}"
                  .format(np.nansum(sio_momentum_blue_north)/xsio/np.cos(inclination),
                          np.nansum(sio_momentum_red_north)/xsio/np.cos(inclination)))
            windspeed = 500*u.km/u.s
            print("north momentum -> mass loss: blue={0:0.3g}, red={1:0.3g}"
                  .format(np.nansum(sio_momentum_blue_north)/xsio/age/windspeed/np.cos(inclination),
                          np.nansum(sio_momentum_red_north)/xsio/age/windspeed/np.cos(inclination)))

            # this is nonsense velmom1_north = np.nansum(velaxis * profile)/np.nansum(profile[np.isfinite(velaxis)])
            # this is nonsense print("north Average Velocity (moment 1): {0:0.3g}".format(velmom1_north))

            nsio_profile = ntot_of_nupper(nupper_of_kkms(u.Quantity(profile, u.K*u.km/u.s),
                                                         ref_freq, 10**aij.mean()),
                                          tbl[-1]['EU_K']*u.K*constants.k_B,
                                          degeneracy=deg,
                                          Q_rot=partfunc(tem=tex),
                                          tex=tex)

            m_sio_profile = (nsio_profile * nr * beam_area / xsio * 2.8*u.Da).to(u.M_sun) / ppbeam
            print(f"north Total mass = {np.nansum(m_sio_profile).to(u.M_sun)}")

            massloss_rate = m_sio_profile / age

            single_event_rate = np.nansum(massloss_rate)

            single_event_mass = np.nansum(m_sio_profile)
            print("north Mass accreted = {0:0.3g} in {1:0.3g} years "
                  "gives rate {2:0.3g} (assuming 100% accreted in single event)"
                  .format(single_event_mass, age, single_event_rate,))

            distances = centers*pixscale
            age_profile = (distances / max_velocity).to(u.yr) / np.tan(inclination)
            massloss_profile = m_sio_profile / age_profile

            integintens_at_max = profile[np.argmin(np.abs(distances-dmax))]
            surfdens_at_max = nsio_profile[np.argmin(np.abs(distances-dmax))]
            print("north dmax = {0:0.3g}".format(dmax.to(u.au)))
            print("north integ intens, Surface density at maxdist = {1:0.3g}, {0:0.3g}".format(surfdens_at_max, integintens_at_max))
            north_blob_area = (np.pi*(0.04*u.arcsec)**2*(dw51)**2).to(u.cm**2, u.dimensionless_angles())
            north_blob_mass = (north_blob_area * surfdens_at_max / xsio * 2.8*u.Da).to(u.M_sun)
            print("Another estimate that is probably not useful, since it assumes "
                  "that a single blob was ejected over its dynamical age (which is "
                  "not self-consistent):")
            print("  north blob area, mass = {0:0.3g}, {1:0.3g}"
                  .format(north_blob_area, north_blob_mass))
            print("  north accr. rate (single blob at age {1:0.3g}) = {0:0.3g}"
                  .format(2*north_blob_mass / age, age))
            print()

            fig1 = pl.figure(1)
            fig1.clf()
            ax = fig1.gca()
            ax.plot(age_profile, massloss_profile)
            ax.set_title('north')
            ax.set_xlabel("Age (yr)")
            ax.set_ylabel("Mass Ejection Rate ($M_\odot$ yr$^{-1}$)")
            fig1.savefig(paths.fpath('longbaseline/sio_masslossrate_history_north.png'))

            fig3 = pl.figure(3)
            fig3.clf()
            ax = fig3.gca()
            ax.set_title('north')
            ax.plot((centers*pixscale).to(u.au), m_sio_profile.to(u.M_sun))
            ax.set_xlabel("Distance from Central Source (AU)")
            ax.set_ylabel("Mass in Outflow ($M_\odot$)")
            fig3.savefig(paths.fpath('longbaseline/sio_massvsdistance_north.png'))

            mxnorth = sm_sio_cube_north.max(axis=0)
            pl.figure(7)
            pl.clf()
            pl.title('north')
            pl.imshow(mxnorth.value, cmap='gray', interpolation='none', origin='lower')
            pl.contour(sio_m0_blue_north.value, colors=['b']*10)
            pl.contour(sio_m0_red_north.value, colors=['r']*10)
            pl.plot([bxmax], [bymax], 'bx')
            pl.plot([rxmax], [rymax], 'rx')
            pl.plot([nc_x], [nc_y], 'wo')
            pl.axis([197, 320, 216, 302])
            pl.savefig(paths.fpath('longbaseline/sio_contour_figure_north.png'))



            e2_ds = paths.dpath('longbaseline/W51e2_siocube_downsampled.fits')
            if os.path.exists(e2_ds):
                sm_sio_cube_e2 = SpectralCube.read(e2_ds)
            elif os.path.exists(e2_ds+".gz"):
                sm_sio_cube_e2 = SpectralCube.read(e2_ds+".gz")
            else:
                siocube = (SpectralCube.read(
                    paths.dpath('longbaseline/linked/W51e2cax.SPW0_ALL_medsub_cutout.fits'))
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
                sm_sio_cube_e2 = siocube.spectral_smooth(kernel).spectral_interpolate(new_xaxis)

                sm_sio_cube_e2.write(e2_ds)

            sm_sio_cube_e2 = sm_sio_cube_e2.to(u.K, sm_sio_cube_e2.beam.jtok_equiv(ref_freq))

            celhdr = sm_sio_cube_e2.wcs.celestial.to_header()
            celhdr['NAXIS1'] = sm_sio_cube_e2.shape[2]
            celhdr['NAXIS2'] = sm_sio_cube_e2.shape[1]
            bluemask = pyregion.open(paths.rpath('sio_blue_boxmask_lb_e2.reg')).as_imagecoord(celhdr).get_mask(sm_sio_cube_e2[0,:,:].hdu)
            redmask = pyregion.open(paths.rpath('sio_red_boxmask_lb_e2.reg')).as_imagecoord(celhdr).get_mask(sm_sio_cube_e2[0,:,:].hdu)

            thrsh = (sm_sio_cube_e2.beam.jtok(ref_freq) / u.Jy * 3*u.mJy).to(u.K)

            sio_m0_blue_e2 = sm_sio_cube_e2.with_mask(sm_sio_cube_e2 > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).moment0()
            sio_m0_red_e2 = sm_sio_cube_e2.with_mask(sm_sio_cube_e2 > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).moment0()
            sio_m1_blue_e2 = sm_sio_cube_e2.with_mask(sm_sio_cube_e2 > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).moment1()
            sio_m1_red_e2 = sm_sio_cube_e2.with_mask(sm_sio_cube_e2 > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).moment1()
            sio_peak_blue_e2 = sm_sio_cube_e2.with_mask(sm_sio_cube_e2 > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).max(axis=0)
            sio_peak_red_e2 = sm_sio_cube_e2.with_mask(sm_sio_cube_e2 > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).max(axis=0)

            print(f"e2 max integrated intensity: blue={np.nanmax(sio_m0_blue_e2):0.3g}, red={np.nanmax(sio_m0_red_e2):0.3g}")

            # DEBUG - trying to understand where the various peaks are
            # pl.contour(sio_peak_red_e2, levels=[5,10,25,50,100,150,200,250,300,350], colors=['r']*20)
            # pl.contour(sio_peak_blue_e2, levels=[5,10,25,50,100,150,200,250,300,350], colors=['b']*20)

            #e2_center = coordinates.SkyCoord('19:23:43.76', '14:30:34.500',
            e2_center = coordinates.SkyCoord('19:23:43.9665', '14:30:34.505',
                                             unit=(u.hour, u.deg), frame='icrs')
            nc_x, nc_y = sm_sio_cube_e2.wcs.celestial.wcs_world2pix(e2_center.ra.deg,
                                                                    e2_center.dec.deg,
                                                                    0)
            yy,xx = np.indices(sio_m0_blue_e2.shape)
            rr = ((yy-nc_y)**2 + (xx-nc_x)**2)**0.5

            nr,centers,profile = image_tools.azimuthalAverage(np.nan_to_num(sio_m0_blue_e2),
                                                              mask=np.isfinite(sio_m0_blue_e2),
                                                              center=(nc_x,nc_y),
                                                              binsize=2,
                                                              return_nr=True)

            vnr,vcenters,vprofile = image_tools.azimuthalAverage(np.nan_to_num(sio_m1_blue_e2),
                                                                 mask=np.isfinite(sio_m1_blue_e2),
                                                                 center=(nc_x,nc_y),
                                                                 binsize=2,
                                                                 return_nr=True)

            pixscale = (wcs.utils.proj_plane_pixel_scales(sm_sio_cube_e2.wcs)[0]*u.deg * dw51).to(u.pc, u.dimensionless_angles())
            inclination = 45*u.deg

            max_velocity = 105*u.km/u.s
            dmax = (0.474*u.arcsec*dw51).to(u.pc, u.dimensionless_angles())
            age = (dmax / max_velocity).to(u.yr) / np.tan(inclination)


            kkms_to_nupper = nupper_of_kkms(1*u.K*u.km/u.s, ref_freq,
                                            10**aij.mean())
            print(f"Conversion from integrated intensity to N_upper is 1 K km/s = {kkms_to_nupper:0.3g}")

            # just in case pixscale changed...
            kkms_to_column = ntot_of_nupper(kkms_to_nupper,
                                            eupper=tbl[-1]['EU_K']*u.K*constants.k_B,
                                            degeneracy=deg,
                                            Q_rot=partfunc(tem=tex), tex=tex)
            # note that this isn't a real mass of anything; M_sio would use 44 Da....
            kkms_to_mass = kkms_to_column * 2.8*u.Da * pixscale**2
            print(f"kkms_to_mass = {kkms_to_mass.to(u.M_sun) / (u.K*u.km/u.s)}")
            dv = np.mean(np.diff(sm_sio_cube_e2.spectral_axis))

            print(f"Conversion from N_upper to N_tot for 1 K km/s at tex={tex} K = {kkms_to_column:0.3g}")

            mass_sm_sio_cube_e2 = ((sm_sio_cube_e2 * kkms_to_mass *
                                    dv) / (u.K*u.km/u.s)).to(u.M_sun)
            e2_outflow_total_mass = mass_sm_sio_cube_e2.sum()/xsio
            print(f"e2 mass for xsio={xsio}: {e2_outflow_total_mass}")
            print(f"e2 peak column: peak kkms={np.nanmax(sio_m0_blue_e2)} -> column={np.nanmax(sio_m0_blue_e2)*kkms_to_column/(u.K*u.km/u.s)}")
            velcenter = 60*u.km/u.s
            velaxis = mass_sm_sio_cube_e2.spectral_axis - velcenter
            sio_momentum_blue_e2 = ((mass_sm_sio_cube_e2 * velaxis[:,None,None])
                                    .with_mask(sm_sio_cube_e2 > thrsh)
                                    .with_mask(bluemask)
                                    .spectral_slab(-140*u.km/u.s, 60*u.km/u.s)
                                    .sum(axis=0))
            sio_momentum_red_e2 = ((mass_sm_sio_cube_e2 * velaxis[:,None,None])
                                   .with_mask(sm_sio_cube_e2 > thrsh)
                                   .with_mask(redmask)
                                   .spectral_slab(60*u.km/u.s, 260*u.km/u.s)
                                   .sum(axis=0))
            print("Velocity inclination factor = 1/cos({0:0.3g}) = {1:0.3g}"
                  .format(inclination, 1/np.cos(inclination)))
            print("e2 total momentum (inclination corrected): blue={0:0.3g}, red={1:0.3g}"
                  .format(np.nansum(sio_momentum_blue_e2)/xsio/np.cos(inclination),
                          np.nansum(sio_momentum_red_e2)/xsio/np.cos(inclination)))
            windspeed = 500*u.km/u.s
            print("e2 momentum -> mass loss: blue={0:0.3g}, red={1:0.3g}  [windspeed={windspeed:0.3g}, age={age:0.3g}]"
                  .format(np.nansum(sio_momentum_blue_e2)/xsio/age/windspeed/np.cos(inclination),
                          np.nansum(sio_momentum_red_e2)/xsio/age/windspeed/np.cos(inclination),
                          windspeed=windspeed, age=age
                         ))

            # this is nonsense velmom1_e2 = np.nansum(velaxis * profile)/np.nansum(profile[np.isfinite(velaxis)])
            # this is nonsense print("e2 Average Velocity (moment 1): {0:0.3g}".format(velmom1_e2))


            print("e2 age estimate from max velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(max_velocity, dmax.to(u.au), age))

            # second age estimate
            # this time, we use the location of the highest value in the moment 0 map
            # to represent the "average" outflowing mass.
            # If we're looking at a single event, i.e., a hubble flow, the age should
            # be the same for the average and the max.
            # If we're looking at a constantly-driven outflow... we should see the same
            # velocities at all points and the age should be lower because we're looking closer
            # to the source

            # compute the positions of brightest integrated emission
            bymax,bxmax = np.unravel_index(np.nanargmax(sio_m0_blue_e2), sio_m0_blue_e2.shape)
            rymax,rxmax = np.unravel_index(np.nanargmax(sio_m0_red_e2), sio_m0_red_e2.shape)

            # find the offset from the central position, the velocity offset, and compute age
            m0pk_sep_blue_e2 = (((bxmax-nc_x)**2 + (bymax-nc_y)**2)**0.5 * pixscale).to(u.au)
            vatmax_blue_e2 = sio_m1_blue_e2[bymax, bxmax] - velcenter
            avg_age_blue_e2 = (m0pk_sep_blue_e2 / vatmax_blue_e2).to(u.yr) / np.tan(inclination)

            m0pk_sep_red_e2 = (((rxmax-nc_x)**2 + (rymax-nc_y)**2)**0.5 * pixscale).to(u.au)
            vatmax_red_e2 = sio_m1_red_e2[rymax, rxmax] - velcenter
            avg_age_red_e2 = (m0pk_sep_red_e2 / vatmax_red_e2).to(u.yr) / np.tan(inclination)

            print("e2 blue age estimate from avg velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(vatmax_blue_e2, m0pk_sep_blue_e2, avg_age_blue_e2))
            print("e2 red age estimate from avg velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(vatmax_red_e2, m0pk_sep_red_e2, avg_age_red_e2))



            ppbeam = (beam_area / pixscale**2).decompose()

            # we assume T=250 K, which is close to the peak brightness temperature
            # The total N(SiO) is approximately linear with temperature in this regime
            # (Q~T), so the errors resulting from this estimate are small
            nsio_profile = ntot_of_nupper(nupper_of_kkms(u.Quantity(profile,u.K*u.km/u.s),
                                                         ref_freq, 10**aij.mean()),
                                          tbl[-1]['EU_K']*u.K*constants.k_B,
                                          degeneracy=deg,
                                          Q_rot=partfunc(tem=tex),
                                          tex=250*u.K)
            print(f"e2 integrated intensity = {np.nanmax(profile)}")

            m_sio_profile = (nsio_profile * nr * beam_area / xsio * 2.8*u.Da).to(u.M_sun) / ppbeam

            # this assumes a hubble flow: all mass sent out at same age
            massloss_rate = m_sio_profile / age

            single_event_rate = np.nansum(massloss_rate)

            single_event_mass = np.nansum(m_sio_profile)
            print("e2 Mass accreted = {0:0.3g} in {1:0.3g} years gives rate {2:0.3g} from the integrated profile of SiO column profile / age"
                  .format(single_event_mass, age, single_event_rate,))

            distances = centers*pixscale
            age_profile = (distances / max_velocity).to(u.yr) / np.tan(inclination)
            massloss_profile = m_sio_profile / age_profile

            integintens_at_max = profile[np.argmin(np.abs(distances-dmax))]
            surfdens_at_max = nsio_profile[np.argmin(np.abs(distances-dmax))]
            print("e2 dmax = {0:0.3g}".format(dmax.to(u.au)))
            print("e2 integ intens, Surface density at maxdist = {1:0.3g}, {0:0.3g}".format(surfdens_at_max, integintens_at_max))
            e2_blob_area = (np.pi*0.058*u.arcsec*0.031*u.arcsec*(dw51)**2).to(u.cm**2,
                                                                                   u.dimensionless_angles())
            e2_blob_mass = (e2_blob_area * surfdens_at_max / xsio * 2.8*u.Da).to(u.M_sun)
            print("Another estimate that is probably not useful, since it assumes "
                  "that a single blob was ejected over its dynamical age (which is "
                  "not self-consistent):")
            print("  e2 blob area, mass = {0:0.3g}, {1:0.3g}".format(e2_blob_area, e2_blob_mass))
            print("  e2 accr. rate from this single blob at age {1:0.3g} = {0:0.3g}".format(2*e2_blob_mass / age, age))
            print()

            fig2 = pl.figure(2)
            fig2.clf()
            ax = fig2.gca()
            ax.set_title('e2')
            ax.plot(age_profile, massloss_profile)
            ax.set_xlabel("Age (yr)")
            ax.set_ylabel("Mass Ejection Rate ($M_\odot$ yr$^{-1}$)")
            ax.set_xlim(0, 200)
            fig2.savefig(paths.fpath('longbaseline/sio_masslossrate_history_e2.png'))


            fig4 = pl.figure(4)
            fig4.clf()
            ax = fig4.gca()
            ax.set_title('e2')
            ax.plot((centers*pixscale).to(u.au), m_sio_profile.to(u.M_sun))
            ax.set_xlabel("Distance from Central Source (AU)")
            ax.set_ylabel("Mass in Outflow ($M_\odot$)")
            ax.set_xlim(0, 3000)
            fig4.savefig(paths.fpath('longbaseline/sio_massvsdistance_e2.png'))

            mxe2 = sm_sio_cube_e2.max(axis=0)
            pl.figure(5)
            pl.clf()
            pl.title('e2')
            pl.imshow(mxe2.value, cmap='gray', interpolation='none', origin='lower')
            pl.contour(sio_m0_blue_e2.value, colors=['b']*10)
            pl.contour(sio_m0_red_e2.value, colors=['r']*10)
            pl.plot([bxmax], [bymax], 'bx')
            pl.plot([rxmax], [rymax], 'rx')
            pl.plot([nc_x], [nc_y], 'wo')
            pl.axis([268, 439, 305, 428])
            pl.savefig(paths.fpath('longbaseline/sio_contour_figure_e2.png'))




            e8_ds = paths.dpath('longbaseline/W51e8_siocube_downsampled.fits')
            if os.path.exists(e8_ds):
                sm_sio_cube_e8 = SpectralCube.read(e8_ds)
            elif os.path.exists(e8_ds+".gz"):
                sm_sio_cube_e8 = SpectralCube.read(e8_ds+".gz")
            else:
                siocube = (SpectralCube.read(
                    paths.dpath('longbaseline/linked/W51e8cax.SPW0_ALL_medsub_cutout.fits'))
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
                sm_sio_cube_e8 = siocube.spectral_smooth(kernel).spectral_interpolate(new_xaxis)

                sm_sio_cube_e8.write(e8_ds)


            sm_sio_cube_e8 = sm_sio_cube_e8.to(u.K, sm_sio_cube_e8.beam.jtok_equiv(ref_freq))

            celhdr = sm_sio_cube_e8.wcs.celestial.to_header()
            celhdr['NAXIS1'] = sm_sio_cube_e8.shape[2]
            celhdr['NAXIS2'] = sm_sio_cube_e8.shape[1]
            bluemask = pyregion.open(paths.rpath('sio_blue_boxmask_lb_e8.reg')).as_imagecoord(celhdr).get_mask(sm_sio_cube_e8[0,:,:].hdu)
            redmask = pyregion.open(paths.rpath('sio_red_boxmask_lb_e8.reg')).as_imagecoord(celhdr).get_mask(sm_sio_cube_e8[0,:,:].hdu)

            thrsh = (sm_sio_cube_e8.beam.jtok(ref_freq) / u.Jy * 3*u.mJy).to(u.K)
            sio_m0_blue_e8 = sm_sio_cube_e8.with_mask(sm_sio_cube_e8 > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).moment0()
            sio_m0_red_e8 = sm_sio_cube_e8.with_mask(sm_sio_cube_e8 > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).moment0()
            sio_m1_blue_e8 = sm_sio_cube_e8.with_mask(sm_sio_cube_e8 > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).moment1()
            sio_m1_red_e8 = sm_sio_cube_e8.with_mask(sm_sio_cube_e8 > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).moment1()
            sio_peak_blue_e8 = sm_sio_cube_e8.with_mask(sm_sio_cube_e8 > thrsh).with_mask(bluemask).spectral_slab(-140*u.km/u.s, 60*u.km/u.s).max(axis=0)
            sio_peak_red_e8 = sm_sio_cube_e8.with_mask(sm_sio_cube_e8 > thrsh).with_mask(redmask).spectral_slab(60*u.km/u.s, 260*u.km/u.s).max(axis=0)

            print(f"e8 max integrated intensity: blue={np.nanmax(sio_m0_blue_e8):0.3g}, red={np.nanmax(sio_m0_red_e8):0.3g}")

            e8_center = coordinates.SkyCoord('19:23:43.9053', '14:30:28.2385',
                                             unit=(u.hour, u.deg), frame='icrs')
            nc_x, nc_y = sm_sio_cube_e8.wcs.celestial.wcs_world2pix(e8_center.ra.deg,
                                                                    e8_center.dec.deg,
                                                                    0)
            yy,xx = np.indices(sio_m0_blue_e8.shape)
            rr = ((yy-nc_y)**2 + (xx-nc_x)**2)**0.5

            nr,centers,profile = image_tools.azimuthalAverage(np.nan_to_num(sio_m0_blue_e8),
                                                              mask=np.isfinite(sio_m0_blue_e8),
                                                              center=(nc_x,nc_y),
                                                              binsize=2,
                                                              return_nr=True)

            vnr,vcenters,vprofile = image_tools.azimuthalAverage(np.nan_to_num(sio_m1_blue_e8),
                                                                 mask=np.isfinite(sio_m1_blue_e8),
                                                                 center=(nc_x,nc_y),
                                                                 binsize=2,
                                                                 return_nr=True)

            pixscale = (wcs.utils.proj_plane_pixel_scales(sm_sio_cube_e8.wcs)[0]*u.deg * dw51).to(u.pc, u.dimensionless_angles())
            # For e8, it's probably better to assume something closer to 75 deg?
            # this is supposed to be nearly edge-on
            inclination = 75*u.deg
            # but we're putting it in the paper as quantity*tan(ii/45deg)
            inclination = 45*u.deg


            # 0.667 contains *all* of the sio (i.e., is the most distant SiO)
            dmax = (0.667*u.arcsec*dw51).to(u.pc, u.dimensionless_angles())
            # this is the extent of the "bright" sio
            dmax = (0.195*u.arcsec*dw51).to(u.pc, u.dimensionless_angles())
            # the approximate distance to the highest-velocity thing
            dmax = (0.072*u.arcsec*dw51).to(u.pc, u.dimensionless_angles())
            max_velocity = 42*u.km/u.s
            age = (dmax / max_velocity).to(u.yr) / np.tan(inclination)


            # just in case pixscale changed...
            kkms_to_mass = ntot_of_nupper(nupper_of_kkms(1*u.K*u.km/u.s, ref_freq,
                                                         10**aij.mean()),
                                          eupper=tbl[-1]['EU_K']*u.K*constants.k_B,
                                          degeneracy=deg,
                                          Q_rot=partfunc(tem=tex),
                                          tex=tex) * 2.8*u.Da * pixscale**2
            print(f"kkms_to_mass = {kkms_to_mass.to(u.M_sun) / (u.K*u.km/u.s)}")
            dv = np.mean(np.diff(sm_sio_cube_e8.spectral_axis))

            mass_sm_sio_cube_e8 = ((sm_sio_cube_e8 * kkms_to_mass *
                                    dv) / (u.K*u.km/u.s)).to(u.M_sun)
            print(f"e8 mass (xsio={xsio}): {mass_sm_sio_cube_e8.sum()/xsio}")
            velcenter = 60*u.km/u.s
            velaxis = mass_sm_sio_cube_e8.spectral_axis - velcenter
            sio_momentum_blue_e8 = ((mass_sm_sio_cube_e8 * velaxis[:,None,None])
                                    .with_mask(sm_sio_cube_e8 > thrsh)
                                    .with_mask(bluemask)
                                    .spectral_slab(-140*u.km/u.s, 60*u.km/u.s)
                                    .sum(axis=0))
            sio_momentum_red_e8 = ((mass_sm_sio_cube_e8 * velaxis[:,None,None])
                                   .with_mask(sm_sio_cube_e8 > thrsh)
                                   .with_mask(redmask)
                                   .spectral_slab(60*u.km/u.s, 260*u.km/u.s)
                                   .sum(axis=0))
            print("Velocity inclination factor = 1/cos({0:0.3g}) = {1:0.3g}"
                  .format(inclination, 1/np.cos(inclination)))
            print("e8 total momentum: blue={0:0.3g}, red={1:0.3g}"
                  .format(np.nansum(sio_momentum_blue_e8)/xsio/np.cos(inclination),
                          np.nansum(sio_momentum_red_e8)/xsio/np.cos(inclination)))
            windspeed = 500*u.km/u.s
            print("e8 momentum -> mass loss: blue={0:0.3g}, red={1:0.3g}"
                  .format(np.nansum(sio_momentum_blue_e8)/xsio/age/windspeed/np.cos(inclination),
                          np.nansum(sio_momentum_red_e8)/xsio/age/windspeed/np.cos(inclination)))

            # this is nonsense velmom1_e8 = np.nansum(velaxis * profile)/np.nansum(profile[np.isfinite(velaxis)])
            # this is nonsense print("e8 Average Velocity (moment 1): {0:0.3g}".format(velmom1_e8))

            print('e8 inclination: {0:0.3g}'.format(inclination))
            print("e8 age estimate from max velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(max_velocity, dmax.to(u.au), age))

            # second age estimate
            # this time, we use the location of the highest value in the moment 0 map
            # to represent the "average" outflowing mass.
            # If we're looking at a single event, i.e., a hubble flow, the age should
            # be the same for the average and the max.
            # If we're looking at a constantly-driven outflow... we should see the same
            # velocities at all points and the age should be lower because we're looking closer
            # to the source

            # compute the positions of brightest integrated emission
            bymax,bxmax = np.unravel_index(np.nanargmax(sio_m0_blue_e8), sio_m0_blue_e8.shape)
            rymax,rxmax = np.unravel_index(np.nanargmax(sio_m0_red_e8), sio_m0_red_e8.shape)

            mxe8 = sm_sio_cube_e8.max(axis=0)
            pl.figure(6)
            pl.clf()
            pl.title('e8')
            pl.imshow(mxe8.value, cmap='gray', interpolation='none', origin='lower')
            pl.contour(sio_m0_blue_e8.value, colors=['b']*10)
            pl.contour(sio_m0_red_e8.value, colors=['r']*10)
            pl.plot([bxmax], [bymax], 'bx')
            pl.plot([rxmax], [rymax], 'rx')
            pl.plot([nc_x], [nc_y], 'wo')
            pl.axis([334, 474, 420, 495])
            pl.savefig(paths.fpath('longbaseline/sio_contour_figure_e8.png'))

            # find the offset from the central position, the velocity offset, and compute age
            m0pk_sep_blue_e8 = (((bxmax-nc_x)**2 + (bymax-nc_y)**2)**0.5 * pixscale).to(u.au)
            vatmax_blue_e8 = sio_m1_blue_e8[bymax, bxmax] - velcenter
            avg_age_blue_e8 = (m0pk_sep_blue_e8 / vatmax_blue_e8).to(u.yr) / np.tan(inclination)

            m0pk_sep_red_e8 = (((rxmax-nc_x)**2 + (rymax-nc_y)**2)**0.5 * pixscale).to(u.au)
            vatmax_red_e8 = sio_m1_red_e8[rymax, rxmax] - velcenter
            avg_age_red_e8 = (m0pk_sep_red_e8 / vatmax_red_e8).to(u.yr) / np.tan(inclination)

            print("e8 blue age estimate from avg velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(vatmax_blue_e8, m0pk_sep_blue_e8, avg_age_blue_e8))
            print("e8 red age estimate from avg velocity={0:0.3g} separation={1:0.3g} age={2:0.3g}"
                  .format(vatmax_red_e8, m0pk_sep_red_e8, avg_age_red_e8))




            ppbeam = (beam_area / pixscale**2).decompose()

            # we assume T=250 K, which is close to the peak brightness temperature
            # The total N(SiO) is approximately linear with temperature in this regime
            # (Q~T), so the errors resulting from this estimate are small
            nsio_profile = ntot_of_nupper(nupper_of_kkms(u.Quantity(profile,u.K*u.km/u.s),
                                                         ref_freq, 10**aij.mean()),
                                          tbl[-1]['EU_K']*u.K*constants.k_B,
                                          degeneracy=deg,
                                          Q_rot=partfunc(tem=tex),
                                          tex=tex)

            m_sio_profile = (nsio_profile * nr * beam_area / xsio * 2.8*u.Da).to(u.M_sun) / ppbeam

            massloss_rate = m_sio_profile / age

            single_event_rate = np.nansum(massloss_rate)

            single_event_mass = np.nansum(m_sio_profile)
            print("e8 Mass accreted = {0:0.3g} in {1:0.3g} years gives rate {2:0.3g}"
                  .format(single_event_mass, age, single_event_rate,))

finally:
    sys.stdout = sys.__stdout__
