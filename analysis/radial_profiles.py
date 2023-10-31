import matplotlib
matplotlib.use('Qt5Agg')
import os
import pylab as pl
assert matplotlib.get_backend() == 'Qt5Agg'
import numpy as np
import radio_beam
import pyregion
import paths
#import image_tools
from image_tools.radialprofile import azimuthalAverage
from astropy import wcs
from astropy import constants
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import units as u
from astropy import coordinates
import itertools
import masscalc
import dust_emissivity
import reproject
from label_lines import labelLine

pl.matplotlib.rc_file('pubfiguresrc')
pl.rcParams['backend'] = 'Qt5Agg'
pl.rcParams['axes.titlesize'] = 12
pl.rcParams['axes.labelsize'] = 12
pl.rcParams['font.size'] = 12
pl.rcParams['figure.dpi'] = 150
figsize=(10,10)
pl.rcParams['figure.figsize'] = figsize

assert matplotlib.get_backend() == 'Qt5Agg'

ffiles = """
selfcal_allspw_mfs.image.pbcor.fits
selfcal_allspw_selfcal_3_mfs_deeper.image.pbcor.fits
selfcal_spw3_selfcal_4ampphase_mfs_tclean.model.fits
selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_10mJy.image.pbcor.fits
w51_cont_spw0_hires.image.fits
w51_cont_spw1_hires.image.fits
w51_cont_spw2_hires.image.fits
w51_cont_spw3.image.fits
w51_cont_spw3_hires.image.fits
w51_spw2_continuum_noflag.image.fits
w51_spw3_continuum.image.fits
w51_spw3_continuum_7m12m.image.fits
w51_spw3_continuum_flagged_tclean.image.fits
w51_spw3_continuum_flagged_uniform_tclean.image.fits
w51_spw3_continuum_noflag.image.fits
w51_spw3_continuum_r0.image.fits
w51_spw3_continuum_r0_mulstiscale.image.fits
""".split()



def make_rprof(regions, ploteach=False):
    names = [r.attr[1]['text'] for r in regions]
    center_positions = coordinates.SkyCoord([r.coord_list
                                             for r in regions],
                                            unit=(u.deg, u.deg),
                                            frame='fk5')

    size = u.Quantity([1.25,1.25], u.arcsec)

    if ploteach:
        nplots = len(names)
        for ii in range(nplots):
            pl.figure(ii, figsize=figsize).clf()
            pl.figure(nplots+ii, figsize=figsize).clf()
            pl.figure(nplots*2+ii, figsize=figsize).clf()

        linestyles = {name: itertools.cycle(['-'] + ['--'] + [':'] + ['-.'])
                      for name in names}

        for fn in ffiles:
            if os.path.exists(paths.dpath("12m/continuum/"+fn+".gz")):
                fh = fits.open(paths.dpath("12m/continuum/"+fn+".gz"))
            else:
                fh = fits.open(paths.dpath("12m/continuum/"+fn))
            mywcs = wcs.WCS(fh[0].header)

            if 'BMAJ' not in fh[0].header:
                #print("File {0} does not have BMAJ".format(fn))
                continue
            try:
                beam = radio_beam.Beam.from_fits_header(fh[0].header)
            except KeyError:
                #print("File {0} doesn't have beam info in the header".format(fn))
                continue

            pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5
            ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam
            #print("fn  {0} ppbeam={1:0.2f}".format(fn, ppbeam))

            for jj,(name,position) in enumerate(zip(names, center_positions)):
                cutout = Cutout2D(fh[0].data, position, size, wcs=mywcs)

                nr, bins, rprof = azimuthalAverage(cutout.data, binsize=1.0,
                                                   return_nr=True)

                linestyle = next(linestyles[name])

                pl.figure(jj, figsize=figsize)
                pl.title(name)
                pl.plot(bins*pixscale*3600., rprof/ppbeam,
                        label=fn.split(".")[0], linestyle=linestyle)
                pl.ylabel("Azimuthally Averaged Flux (Jy)")
                pl.xlabel("Radius [arcsec]")

                cumul_rprof = np.nan_to_num(rprof*nr/ppbeam).cumsum()

                pl.figure(nplots+jj, figsize=figsize)
                pl.title(name)
                pl.plot(bins*pixscale*3600., cumul_rprof,
                        label=fn.split(".")[0], linestyle=linestyle)
                if jj == 0:
                    pl.fill_between([0, beam.major.to(u.arcsec).value], [100,100], [0,0], zorder=-5, alpha=0.1, color='k')
                pl.ylabel("Cumulative Flux (Jy)")
                pl.xlabel("Radius [arcsec]")
                if jj == 0:
                    ax = pl.gca()
                    ax2 = ax.twiny()
                    ax3 = ax.twinx()
                    def tick_function(old_x):
                        newx = (old_x*u.arcsec*masscalc.distance).to(u.pc, u.dimensionless_angles()).value
                        return ["%.1f" % z for z in newx]
                    new_tick_locations = [0.005,0.01,0.015,0.02,0.025]*u.pc
                    new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
                    ax2.set_xlim(ax.get_xlim())
                    ax2.set_xticks(new_tick_locs_as.value)
                    ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
                    ax2.set_xlabel(r"Radius (pc)")
                    ax3.set_ylim(ax.get_ylim())
                    yticks_mass = np.arange(0,6000,1000)
                    yticks_Jy = yticks_mass/masscalc.mass_conversion_factor().value
                    ax3.set_yticks(yticks_Jy)
                    ax3.set_yticklabels(yticks_mass)
                    ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=20$ K)")

                pl.figure(nplots*2+jj, figsize=figsize)
                pl.title(name)
                pl.plot(((bins*pixscale*u.deg)*masscalc.distance).to(u.pc,
                                                                     u.dimensionless_angles()),
                        cumul_rprof * masscalc.mass_conversion_factor(),
                        label=fn.split(".")[0], linestyle=linestyle)
                if jj == 0:
                    pl.fill_between([0, beam.major.to(u.arcsec).value], [100,100], [0,0], zorder=-5, alpha=0.1, color='k')
                pl.ylabel("Cumulative Mass (M$_\\odot$, $T=20$ K)")
                pl.xlabel("Radius (pc)")

        for ii in range(nplots):
            for xtra in (0,nplots*2):
                ax = pl.figure(ii+xtra, figsize=figsize).gca()
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

                # Put a legend to the right of the current axis
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        nplots = 0

    for jj in range(nplots*3+1, nplots*3+15+1):
        pl.figure(jj, figsize=figsize).clf()
    # $ find ~/work/w51/alma/FITS/ -samefile ~/work/w51/alma/FITS/W51_te_continuum_best.fits
    # /Users/adam/work/w51/alma/FITS//12m/continuum/selfcal_allspw_selfcal_3_mfs_deeper.image.pbcor.fits
    # /Users/adam/work/w51/alma/FITS//W51_te_continuum_best.fits
    fn = "selfcal_allspw_selfcal_3_mfs_deeper.image.pbcor.fits"
    fh = fits.open(paths.dpath("12m/continuum/"+fn))
    mywcs = wcs.WCS(fh[0].header)
    beam = radio_beam.Beam.from_fits_header(fh[0].header)
    pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5
    ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam
    for ii,(name,position) in enumerate(zip(names, center_positions)):
        cutout = Cutout2D(fh[0].data, position, size, wcs=mywcs)

        nr, bins, rprof = azimuthalAverage(cutout.data, binsize=1.0,
                                           return_nr=True)

        pl.figure(nplots*3+1, figsize=figsize)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(bins*pixscale*3600., rprof/ppbeam,
                label=name)
        pl.ylabel("Azimuthally Averaged Flux (Jy)")
        pl.xlabel("Radius [arcsec]")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        cumul_rprof = np.nan_to_num(rprof*nr/ppbeam).cumsum()

        pl.figure(nplots*3+2, figsize=figsize)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(bins*pixscale*3600., cumul_rprof,
                label=name)
        pl.ylabel("Cumulative Flux (Jy)")
        pl.xlabel("Radius [arcsec]")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if ii == len(names) - 1:
            ax = pl.gca()
            ax2 = ax.twiny()
            ax3 = ax.twinx()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")
            ax3.set_ylim(ax.get_ylim())
            yticks_mass = np.arange(0,6000,1000)
            yticks_Jy = yticks_mass/masscalc.mass_conversion_factor(TK=40).value
            ax3.set_yticks(yticks_Jy)
            ax3.set_yticklabels(yticks_mass)
            ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=40$ K)")



        radii = ((bins*pixscale*u.deg)*masscalc.distance).to(u.pc,
                                                             u.dimensionless_angles())
        mass_40k_profile = (cumul_rprof * masscalc.mass_conversion_factor(TK=40) / u.beam).to(u.M_sun)

        pl.figure(nplots*3+3, figsize=figsize)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(radii, mass_40k_profile, label=name)
        pl.ylabel("Cumulative Mass (M$_\\odot$, $T=40$ K)")
        pl.xlabel("Radius (pc)")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


        density_40k_profile = (mass_40k_profile / (4/3.*np.pi*radii**3) / (2.8*u.Da)).to(u.cm**-3)
        #print(density_40k_profile)
        #print(radii)

        pl.figure(nplots*3+4, figsize=figsize)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.semilogy(radii, density_40k_profile, label=name)
        pl.ylabel("Cumulative Density [n(H$_2$), $T=40$ K]")
        pl.xlabel("Radius (pc)")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        angular_radii = bins*pixscale*3600.
        azimuthal_average_flux = rprof/ppbeam
        azimuthal_average_mass_40K = azimuthal_average_flux * masscalc.mass_conversion_factor(TK=40) / u.beam
        bindiff_as = np.diff(bins).mean() * pixscale * 3600.
        bindiff_cm = bindiff_as * masscalc.distance / 206265.
        bins_cm = (angular_radii * masscalc.distance / 206265.).to(u.cm)
        sqdiffbins_cm = u.Quantity([bins_cm[0].to(u.cm).value**2] +
                                   (bins_cm.to(u.cm)[1:]**2 -
                                    bins_cm.to(u.cm)[:-1]**2).value.tolist(),
                                   u.cm**2)
        cudiffbins_cm = u.Quantity([bins_cm[0].to(u.cm).value**3] +
                                   (bins_cm.to(u.cm)[1:]**3 -
                                    bins_cm.to(u.cm)[:-1]**3).value.tolist(),
                                   u.cm**3)
        sqdiffbins_pix = [bins[0]**2] + (bins[1:]**2 - bins[:-1]**2).tolist()
        azimuthal_average_density_40K = (azimuthal_average_mass_40K * sqdiffbins_pix /
                                        (2.8*u.Da) /
                                        (4/3.*np.pi*(cudiffbins_cm))).to(u.cm**-3)

        pl.figure(nplots*3+5, figsize=figsize)
        pl.semilogy(angular_radii, azimuthal_average_density_40K, label=name)
        pl.ylabel("Azimuthally Averaged Density\nat 40K [cm$^{-3}$]")
        pl.xlabel("Radius [arcsec]")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if ii == len(names) - 1:
            ax = pl.gca()
            ax.set_ylim(1e7, 5e10)
            ax2 = ax.twiny()
            #ax3 = ax.twinx()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")
            #ax3.set_ylim(ax.get_ylim())
            #yticks_mass = np.arange(0,6000,1000)
            #yticks_Jy = yticks_mass/masscalc.mass_conversion_factor(TK=40).value
            #ax3.set_yticks(yticks_Jy)
            #ax3.set_yticklabels(yticks_mass)
            #ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=40$ K)")



        c_s_40k = ((constants.k_B * 40*u.K / (2.4*u.Da))**0.5).to(u.km/u.s)
        azimuthal_average_MJ_40k = (np.pi/6. * c_s_40k**3 /
                                    (constants.G**1.5 *
                                     (2.8*u.Da*azimuthal_average_density_40K)**0.5)).to(u.M_sun)

        pl.figure(nplots*3+6, figsize=figsize)
        pl.semilogy(angular_radii, azimuthal_average_MJ_40k, label=name)
        pl.ylabel("Azimuthally Averaged $M_J$\nat $T=40$K")
        pl.xlabel("Radius [arcsec]")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if ii == len(names) - 1:
            ax = pl.gca()
            #ax.set_ylim(1e7, 5e10)
            ax2 = ax.twiny()
            #ax3 = ax.twinx()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")

        source = 'e2' if name == 'e2e' else name
        temperature_map_fn = paths.dpath('12m/moments/CH3OH_{0}_cutout_temperaturemap.fits'.format(source))
        temperature_map_fh = fits.open(temperature_map_fn)

        # this whole section is copied from overlay_contours_on_ch3oh
        ch3ohN_hdul = fits.open(paths.dpath('12m/moments/CH3OH_{0}_cutout_columnmap.fits'.format(source)))
        ch3ohT_hdul = fits.open(paths.dpath('12m/moments/CH3OH_{0}_cutout_temperaturemap.fits'.format(source)))
        bigwcs = wcs.WCS(ch3ohT_hdul[0].header)
        bigpixscale = (bigwcs.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg
        ch3ohN = ch3ohN_hdul[0].data
        ch3ohT = ch3ohT_hdul[0].data
        dust_brightness,wts = reproject.reproject_interp(fits.open(paths.dpath('W51_te_continuum_best.fits')),
                                                         ch3ohN_hdul[0].header)
        contbm = radio_beam.Beam.from_fits_header(paths.dpath("W51_te_continuum_best.fits"))

        yy,xx = np.indices(ch3ohN.shape)
        if source == 'north':
            center = [84.,38.]
        else:
            center = [ch3ohN.shape[0]/2., ch3ohN.shape[1]/2.]
        yyc = (yy-center[0])
        xxc = (xx-center[1])
        rr = (yyc**2 + xxc**2)**0.5
        rr_as = (rr*bigpixscale).to(u.arcsec)
        theta = np.arctan2(yyc,xxc)*u.rad

        dust_column = dust_emissivity.dust.colofsnu(nu=225*u.GHz, snu_per_beam=dust_brightness*u.Jy/contbm.sr,
                                                    #beamomega=contbm,
                                                    temperature=ch3ohT*u.K)
        ch3oh_abundance = ch3ohN / dust_column.value
        mask = (ch3oh_abundance > 1e-10) & (ch3oh_abundance < 1e-5)
        if source == 'e2':
            mask = mask & (((theta > 15*u.deg) & (theta < 345*u.deg)) | (theta < -15*u.deg))
        mask = mask & np.isfinite(ch3oh_abundance)
        # exclude high-abundance, low-column regions: likely to be div-by-zero zones
        mask = mask & (~((ch3ohN < 1e18) & (ch3oh_abundance > 5e-6)))
        mask = mask & (~((dust_brightness<1e-2) & (ch3ohT > 500) & (ch3oh_abundance > 1e-6)))

        fits.PrimaryHDU(data=dust_column.to(u.cm**-2).value,
                        header=ch3ohN_hdul[0].header).writeto(f'dust_column_density_with_ch3oh_tem_{source}.fits', overwrite=True)

        #mask = mask & (~((ch3ohT > 250) &
        #                 (ch3ohN < 1e18) &
        #                 (rr_as>1.5*u.arcsec))
        #               )# these are low-column,

        temwcs = wcs.WCS(temperature_map_fh[0].header)

        temcutout = Cutout2D(temperature_map_fh[0].data, position, size, wcs=temwcs)
        maskcutout = Cutout2D(mask.astype('float'), position, size, wcs=bigwcs)

        tem_nr, tem_bins, tem_rprof = azimuthalAverage(temcutout.data,
                                                       weights=(maskcutout.data>0).astype('float'),
                                                       binsize=1.0,
                                                       return_nr=True,
                                                       interpnan=True,
                                                      )
        azimuthal_average_temperature = tem_rprof

        mass_map = (cutout.data * masscalc.mass_conversion_factor(TK=temcutout.data) ).to(u.M_sun)
        yy,xx = np.indices(cutout.data.shape)
        rr_map = ((xx-cutout.data.shape[1]/2.)**2 + (yy-cutout.data.shape[0]/2.)**2)
        mass_profile = (cumul_rprof *
                        masscalc.mass_conversion_factor(TK=tem_rprof) /
                        u.beam).to(u.M_sun)
        density_profile = (mass_profile / (4/3.*np.pi*radii**3) /
                           (2.8*u.Da)).to(u.cm**-3)

        azimuthal_average_mass = (azimuthal_average_flux *
                                  masscalc.mass_conversion_factor(TK=tem_rprof)
                                  / u.beam)

        # how accurate is this?  We are measuring mass in cylindrical annuli,
        # but we are assuming the underlying source structure is spherical
        # See Cores.ipynb (in notes/, not in this repo) for a detailed
        # analysis...
        azimuthal_average_density = (azimuthal_average_mass * sqdiffbins_pix /
                                     (2.8*u.Da) /
                                     (4/3.*np.pi*(cudiffbins_cm))).to(u.cm**-3)


        # this is not a good approximation for the spherial density profile...
        # see cores.ipynb.
        # Really should be this number +1 (so see the +1 below)
        # Also, negatived...
        density_alpha_ = 1-((np.log(density_profile[1:].value) -
                             np.log(density_profile[:-1].value)) /
                            (np.log(angular_radii[1:]) -
                             np.log(angular_radii[:-1]))
                           )
        # from cores.ipynb, this actually gets (close to) the right answer
        density_alpha = -((np.log(azimuthal_average_density[1:].value) -
                           np.log(azimuthal_average_density[:-1].value)) /
                          (np.log(angular_radii[1:]) -
                           np.log(angular_radii[:-1]))
                         )

        c_s = ((constants.k_B * tem_rprof*u.K / (2.4*u.Da))**0.5).to(u.km/u.s)
        azimuthal_average_MJ = (np.pi/6. * c_s**3 /
                                (constants.G**1.5 *
                                 (2.8*u.Da*azimuthal_average_density)**0.5)
                               ).to(u.M_sun)

        cumulative_mass_weighted_temperature = ((azimuthal_average_mass *
                                                 azimuthal_average_temperature).cumsum()
                                                / mass_profile)
        cumulative_c_s = ((constants.k_B *
                           cumulative_mass_weighted_temperature*u.K /
                           (2.4*u.Da))**0.5).to(u.km/u.s)
        cumulative_profile_MJ = (np.pi/6. * cumulative_c_s**3 /
                                (constants.G**1.5 *
                                 (2.8*u.Da*density_profile)**0.5)
                               ).to(u.M_sun)



        pl.figure(nplots*3+7, figsize=figsize)
        line, = pl.plot(angular_radii, azimuthal_average_MJ, label=name)
        # WRONG pl.plot(angular_radii, azimuthal_average_mass, linestyle='--',
        # WRONG         alpha=0.75, color=line.get_color())
        #pl.plot(rr_map*pixscale*3600., mass_map, 'k.', alpha=0.25)
        pl.axis([0,1.2,0.0,5.0])
        pl.ylabel("Azimuthally Averaged $M_J$\nat $T(\\mathrm{CH}_3\\mathrm{OH})$ [$M_\\odot$]")
        pl.xlabel("Radius [arcsec]")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if ii == len(names) - 1:
            ax = pl.gca()
            #ax.set_ylim(1e7, 5e10)
            ax2 = ax.twiny()
            #ax3 = ax.twinx()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")


        pl.figure(nplots*3+8, figsize=figsize)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(bins*pixscale*3600., mass_profile,
                label=name)
        if ii == 0:
            pl.fill_between([0, beam.major.to(u.arcsec).value], [100,100], [0,0], zorder=-5, alpha=0.1, color='k')
        pl.ylabel("Cumulative Mass at T(CH$_3$OH) [M$_\\odot$]")
        pl.xlabel("Radius [arcsec]")
        if len(names) < 5:
            pl.legend(loc='upper left')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if ii == len(names) - 1:
            ax = pl.gca()
            ax2 = ax.twiny()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")

        pl.figure(nplots*3+9, figsize=figsize)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(bins*pixscale*3600., tem_rprof,
                label=name)
        pl.legend(loc='best')
        pl.ylabel("CH$_3$OH temperature")
        pl.xlabel("Radius [as]")


        # Jeans length (radius = length/2)
        azimuthal_average_RJ = (c_s**1 /
                                (constants.G**0.5 *
                                 (2.8*u.Da*azimuthal_average_density)**0.5)
                               ).to(u.au) / 2.

        pl.figure(nplots*3+10, figsize=figsize)
        pl.plot(angular_radii, azimuthal_average_RJ, label=name)
        pl.plot([0,(7000*u.au/masscalc.distance).to(u.arcsec,
                                                    u.dimensionless_angles()).value],
                [0,7000], 'k--', alpha=0.5)
        pl.plot([beam.major.to(u.arcsec).value,
                 (7000*u.au/masscalc.distance).to(u.arcsec,
                                                  u.dimensionless_angles()).value],
                [1000,1000], 'k:', alpha=0.5)
        if ii == 0:
            pl.fill_between([0, beam.major.to(u.arcsec).value], [8000,8000], zorder=-5, alpha=0.1, color='k')
        pl.ylabel("Azimuthally Averaged $R_\\mathrm{J}$\nat $T(\\mathrm{CH}_3\\mathrm{OH})$ [au]")
        pl.xlabel("Radius [arcsec]")
        pl.ylim(0,8000)
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if ii == len(names) - 1:
            ax = pl.gca()
            #ax.set_ylim(1e7, 5e10)
            ax2 = ax.twiny()
            #ax3 = ax.twinx()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")


        pl.figure(nplots*3+11, figsize=figsize)
        line, = pl.plot((angular_radii[1:]+angular_radii[:-1])/2., density_alpha)
        pl.plot((angular_radii[1:]+angular_radii[:-1])/2., density_alpha_, '--',
                color=line.get_color(), alpha=0.5)
        pl.xlabel("Radius [as]")
        pl.ylabel("Density Profile Exponent $\\kappa_\\rho$")


        pl.figure(nplots*3+12, figsize=figsize)
        pl.semilogy(angular_radii, azimuthal_average_density, label=name)
        pl.ylabel("Azimuthally Averaged Density\nat T(CH$_3$OH) [cm$^{-3}$]")
        pl.xlabel("Radius [arcsec]")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if ii == len(names) - 1:
            ax = pl.gca()
            ax.set_ylim(1e6, 2e9)
            ax.set_xlim(0,1.2)
            ax2 = ax.twiny()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,7000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")

            ax3 = ax.twinx()
            def tick_function(old_x):
                newx = ((constants.G * 2.8*u.Da * (old_x*u.cm**-3))**-0.5).to(u.kyr).value
                return ["%.1f" % z for z in newx]
            new_tick_locations = [1.3,1.5,2,3,5]*u.kyr
            new_tick_locs_dens = ((new_tick_locations**2 * constants.G *
                                   2.8*u.Da)**-1).to(u.cm**-3)
            ax3.set_ylim(ax.get_ylim())
            ax3.set_yticks(new_tick_locs_dens.value)
            ax3.set_yticklabels(tick_function(new_tick_locs_dens.value))
            ax3.set_ylabel(r"$t_{ff}$ [kyr]")

            # start from 0.05 arcsec
            xlim = u.Quantity([0.05, ax.get_xlim()[1]], u.arcsec)
            xax_as = np.linspace(xlim[0], xlim[1], 100)
            xax_m = (xax_as *
                     masscalc.distance).to(u.m, u.dimensionless_angles())
            yax_s = ((ax.get_ylim() * u.cm**-3 * 2.8 * u.Da * constants.G)**-0.5).to(u.s)

            # plot 1 km/s line
            line_1 = (((xax_m / (7.5*u.km/u.s))**2 * constants.G *
                       2.8*u.Da)**-1).to(u.cm**-3)
            line_2 = (((xax_m / (5*u.km/u.s))**2 * constants.G *
                       2.8*u.Da)**-1).to(u.cm**-3)
            line_3 = (((xax_m / (10*u.km/u.s))**2 * constants.G *
                       2.8*u.Da)**-1).to(u.cm**-3)

            line, = pl.plot(xax_as.value, line_1.value, 'k:',
                            label='7.5 km s$^{-1}$',
                            zorder=-10)
            labelLine(line, 0.45, "7.5 km s$^{-1}$")
            line, = pl.plot(xax_as.value, line_2.value, 'k:',
                            label='5 km s$^{-1}$',
                            zorder=-10)
            labelLine(line, 0.6, "5 km s$^{-1}$")
            line, = pl.plot(xax_as.value, line_3.value, 'k-.',
                            label='10 km s$^{-1}$',
                            zorder=-10)
            labelLine(line, 0.7, "10 km s$^{-1}$")

            #ax.set_xlim(0,1.2)


            #ax3.set_ylim(ax.get_ylim())
            #yticks_mass = np.arange(0,6000,1000)
            #yticks_Jy = yticks_mass/masscalc.mass_conversion_factor(TK=40).value
            #ax3.set_yticks(yticks_Jy)
            #ax3.set_yticklabels(yticks_mass)
            #ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=40$ K)")

        pl.figure(nplots*3+13, figsize=figsize)
        pl.plot(angular_radii, cumulative_profile_MJ, label=name)
        pl.ylabel("$M_J(r<R)$ at T(CH$_3$OH) [$M_\odot$]")
        pl.xlabel("Radius [arcsec]")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if ii == len(names) - 1:
            ax = pl.gca()
            ax.set_ylim(0, 2)
            ax2 = ax.twiny()
            #ax3 = ax.twinx()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")
            #ax3.set_ylim(ax.get_ylim())
            #yticks_mass = np.arange(0,6000,1000)
            #yticks_Jy = yticks_mass/masscalc.mass_conversion_factor(TK=40).value
            #ax3.set_yticks(yticks_Jy)
            #ax3.set_yticklabels(yticks_mass)
            #ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=40$ K)")

        pl.figure(nplots*3+14)
        line, = pl.plot(angular_radii, azimuthal_average_MJ, label=name)
        pl.plot(angular_radii, azimuthal_average_mass,
                linestyle='--', color=line.get_color())
        if ii == 0:
            pl.fill_between([0, beam.major.to(u.arcsec).value], [100,100], [0,0], zorder=-5, alpha=0.1, color='k')
        pl.ylabel("Mass at T(CH$_3$OH) [$M_\odot$]")
        pl.xlabel("Radius [arcsec]")

        handles, labels = ax.get_legend_handles_labels()
        if '$\\bar{M}$' not in labels:
            handles.append(pl.matplotlib.lines.Line2D([],[],
                                                      color='k', linestyle='--'))
            labels.append('$\\bar{M}$')
        if '$M_J$' not in labels:
            handles.append(pl.matplotlib.lines.Line2D([],[],
                                                      color='k', linestyle='-'))
            labels.append('$M_J$')


        if ii == len(names) - 1:
            ax = pl.gca()
            ax.set_ylim(0, 5)
            ax2 = ax.twiny()
            #ax3 = ax.twinx()
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%i" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius [au]")
            #ax3.set_ylim(ax.get_ylim())
            #yticks_mass = np.arange(0,6000,1000)
            #yticks_Jy = yticks_mass/masscalc.mass_conversion_factor(TK=40).value
            #ax3.set_yticks(yticks_Jy)
            #ax3.set_yticklabels(yticks_mass)
            #ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=40$ K)")

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])
            ax2.set_position([box.x0, box.y0, box.width * 1.0, box.height])

            # Put a legend to the right of the current axis
            ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))


        pl.figure(nplots*3+15, figsize=figsize)
        pl.subplot(2,2,1)
        pl.semilogy(angular_radii, azimuthal_average_density, label=name)
        pl.ylabel("Azimuthally Averaged Density\nat T(CH$_3$OH) [cm$^{-3}$]")
        pl.subplot(2,2,3)
        pl.semilogy(angular_radii, density_profile, label=name)
        pl.ylabel("Cumulative Average Density\nat T(CH$_3$OH) [cm$^{-3}$]")
        pl.xlabel("Radius [arcsec]")
        pl.subplot(2,2,2)
        pl.plot(azimuthal_average_density, density_profile, label=name)
        pl.plot(sorted(azimuthal_average_density.value),
                sorted(azimuthal_average_density.value), 'k--')
        pl.xlabel("Azimuthally Averaged Density\nat T(CH$_3$OH) [cm$^{-3}$]")
        pl.ylabel("Cumulative Average Density\nat T(CH$_3$OH) [cm$^{-3}$]")
        pl.subplot(2,2,4)
        pl.ylabel("Azimuthal / Cumulative")
        pl.plot(angular_radii, azimuthal_average_density / density_profile,
                label=name)
        #pl.plot(angular_radii, [1.0] * len(angular_radii), 'k--')
        pl.xlabel("Radius [arcsec]")

    return locals()

regions = pyregion.open(paths.rpath("hmcore_centroids.reg"))
variables = make_rprof(regions, ploteach=True)
for k,v in variables.items():
    globals()[k] = v
nplots = len(regions)
for suffix in ("png","pdf"):
    pl.figure(nplots*3+2).savefig(paths.fpath("cumulative_radial_flux_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+3).savefig(paths.fpath("cumulative_radial_mass_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+4).savefig(paths.fpath("cumulative_density_40K_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+5).savefig(paths.fpath("azimuthalaverage_density_40K_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+12).savefig(paths.fpath("azimuthalaverage_density_of_TCH3OH_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+6).savefig(paths.fpath("azimuthalaverage_radial_mj_40K_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+7).savefig(paths.fpath("azimuthalaverage_radial_mj_of_TCH3OH_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+8).savefig(paths.fpath("cumulative_radial_mass_of_TCH3OH_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+9).savefig(paths.fpath("azimuthalaverage_radial_TCH3OH_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+10).savefig(paths.fpath("azimuthalaverage_radial_rj_of_TCH3OH_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+11).savefig(paths.fpath("radialprofileexponent_of_TCH3OH_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+13).savefig(paths.fpath("cumulative_MJ_of_TCH3OH_massivecores.{0}".format(suffix)), bbox_inches='tight')
    #assert all(pl.figure(nplots*3+14).get_size_inches() == (12,8))
    pl.figure(nplots*3+14).savefig(paths.fpath("azimuthalaverage_radial_mj_and_mbar_of_TCH3OH_massivecores.{0}".format(suffix)), bbox_inches='tight')
    pl.figure(nplots*3+15).savefig(paths.fpath("compare_cumulative_to_average_density.{0}".format(suffix)), bbox_inches='tight')

fig = pl.figure(nplots*3+8, figsize=figsize)
ax = fig.axes[0]
xx = np.linspace(0,1.2)
# power of R is +3 for volume conversion
ax.plot(xx, 300*xx**1.0, 'k:', label='$\\rho\\propto R^{-2}$', alpha=0.5)
ax.plot(xx, 300*xx**1.5, 'k--', label='$\\rho\\propto R^{-3/2}$', alpha=0.5)
ax.plot(xx, 300*xx**2.0, 'k-', label='$\\rho\\propto R^{-1}$', alpha=0.5)
pl.legend(loc='upper left')
fig.savefig(paths.fpath("cumulative_radial_mass_of_TCH3OH_massivecores_withmodel.png"), bbox_inches='tight')

# g31.41 overplot
r_g31 = np.array([0.11,0.33,0.55,0.77,0.99,1.21,1.43,1.65]) * 7.9/5.4
t_g31 = [544,436,343,233,182,136,110,114]
fig = pl.figure(nplots*3+9)
ax = fig.gca()
ax.plot(r_g31, t_g31, color=(0.2,0.9,0.2), label='G31.41')
pl.legend(loc='best')
fig.savefig(paths.fpath("azimuthalaverage_radial_TCH3OH_massivecores_withG31.png"), bbox_inches='tight')

#regions = pyregion.open(paths.rpath("cores.reg"))
#make_rprof(regions, ploteach=False)
#nplots = 3
#pl.figure(nplots*3+2).savefig(paths.fpath("cumulative_radial_flux_allcores.png"))
#pl.figure(nplots*3+3).savefig(paths.fpath("cumulative_radial_mass_allcores.png"))

pl.draw()
pl.show()
