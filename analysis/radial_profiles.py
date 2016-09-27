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
import pylab as pl
import itertools
import masscalc
import dust_emissivity
import reproject

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
            pl.figure(ii).clf()
            pl.figure(nplots+ii).clf()
            pl.figure(nplots*2+ii).clf()

        linestyles = {name: itertools.cycle(['-'] + ['--'] + [':'] + ['-.'])
                      for name in names}

        for fn in ffiles:
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
            
            for ii,(name,position) in enumerate(zip(names, center_positions)):
                cutout = Cutout2D(fh[0].data, position, size, wcs=mywcs)

                nr, bins, rprof = azimuthalAverage(cutout.data, binsize=1.0,
                                                   return_nr=True)

                linestyle = next(linestyles[name])

                pl.figure(ii)
                pl.title(name)
                pl.plot(bins*pixscale*3600., rprof/ppbeam,
                        label=fn.split(".")[0], linestyle=linestyle)
                pl.ylabel("Azimuthally Averaged Flux (Jy)")
                pl.xlabel("Radius (arcsec)")

                cumul_rprof = np.nan_to_num(rprof*nr/ppbeam).cumsum()

                pl.figure(nplots+ii)
                pl.title(name)
                pl.plot(bins*pixscale*3600., cumul_rprof,
                        label=fn.split(".")[0], linestyle=linestyle)
                pl.ylabel("Cumulative Flux (Jy)")
                pl.xlabel("Radius (arcsec)")
                if ii == 0:
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

                pl.figure(nplots*2+ii)
                pl.title(name)
                pl.plot(((bins*pixscale*u.deg)*masscalc.distance).to(u.pc,
                                                                     u.dimensionless_angles()),
                        cumul_rprof * masscalc.mass_conversion_factor(),
                        label=fn.split(".")[0], linestyle=linestyle)
                pl.ylabel("Cumulative Mass (M$_\\odot$, $T=20$ K)")
                pl.xlabel("Radius (pc)")

        for ii in range(nplots):
            for xtra in (0,nplots*2):
                ax = pl.figure(ii+xtra).gca()
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

                # Put a legend to the right of the current axis
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        nplots = 0

    pl.matplotlib.rc_file('pubfiguresrc')
    for jj in range(nplots*3+1, nplots*3+9):
        pl.figure(jj).clf()
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

        pl.figure(nplots*3+1)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(bins*pixscale*3600., rprof/ppbeam,
                label=name)
        pl.ylabel("Azimuthally Averaged Flux (Jy)")
        pl.xlabel("Radius (arcsec)")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        cumul_rprof = np.nan_to_num(rprof*nr/ppbeam).cumsum()

        pl.figure(nplots*3+2)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(bins*pixscale*3600., cumul_rprof,
                label=name)
        pl.ylabel("Cumulative Flux (Jy)")
        pl.xlabel("Radius (arcsec)")
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
                return ["%.1f" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius (au)")
            ax3.set_ylim(ax.get_ylim())
            yticks_mass = np.arange(0,6000,1000)
            yticks_Jy = yticks_mass/masscalc.mass_conversion_factor(TK=40).value
            ax3.set_yticks(yticks_Jy)
            ax3.set_yticklabels(yticks_mass)
            ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=40$ K)")



        radii = ((bins*pixscale*u.deg)*masscalc.distance).to(u.pc,
                                                             u.dimensionless_angles())
        mass_40k_profile = (cumul_rprof * masscalc.mass_conversion_factor(TK=40) / u.beam).to(u.M_sun)

        pl.figure(nplots*3+3)
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

        pl.figure(nplots*3+4)
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
        azimuthal_average_mass = azimuthal_average_flux * masscalc.mass_conversion_factor(TK=40) / u.beam
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
        azimuthal_average_density = (azimuthal_average_mass * sqdiffbins_pix /
                                     (2.8*u.Da) /
                                     (4/3.*np.pi*(cudiffbins_cm))).to(u.cm**-3)

        pl.figure(nplots*3+5)
        pl.semilogy(angular_radii, azimuthal_average_density, label=name)
        pl.ylabel("Azimuthally Averaged Density [cm$^{-3}$]")
        pl.xlabel("Radius (arcsec)")
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
                return ["%.1f" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius (au)")
            #ax3.set_ylim(ax.get_ylim())
            #yticks_mass = np.arange(0,6000,1000)
            #yticks_Jy = yticks_mass/masscalc.mass_conversion_factor(TK=40).value
            #ax3.set_yticks(yticks_Jy)
            #ax3.set_yticklabels(yticks_mass)
            #ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=40$ K)")


        c_s_40k = ((constants.k_B * 40*u.K / (2.4*u.Da))**0.5).to(u.km/u.s)
        azimuthal_average_MJ_40k = (np.pi/6. * c_s_40k**3 /
                                    (constants.G**1.5 *
                                     (2.8*u.Da*azimuthal_average_density)**0.5)).to(u.M_sun)

        pl.figure(nplots*3+6)
        pl.semilogy(angular_radii, azimuthal_average_MJ_40k, label=name)
        pl.ylabel("Azimuthally Averaged $M_J$ at $T=40$K")
        pl.xlabel("Radius (arcsec)")
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
                return ["%.1f" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius (au)")

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
        bm = radio_beam.Beam.from_fits_header(paths.dpath("W51_te_continuum_best.fits"))

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

        dust_column = dust_emissivity.dust.colofsnu(225*u.GHz, dust_brightness*u.Jy,
                                                    beamomega=bm,
                                                    temperature=ch3ohT*u.K)
        ch3oh_abundance = ch3ohN / dust_column.value
        mask = (ch3oh_abundance > 1e-10) & (ch3oh_abundance < 1e-5)
        if source == 'e2':
            mask = mask & (((theta > 15*u.deg) & (theta < 345*u.deg)) | (theta < -15*u.deg))
        mask = mask & np.isfinite(ch3oh_abundance)
        # exclude high-abundance, low-column regions: likely to be div-by-zero zones
        mask = mask & (~((ch3ohN < 1e18) & (ch3oh_abundance > 5e-6)))
        mask = mask & (~((dust_brightness<1e-2) & (ch3ohT > 500) & (ch3oh_abundance > 1e-6)))

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

        mass_profile = (cumul_rprof * masscalc.mass_conversion_factor(TK=tem_rprof) / u.beam).to(u.M_sun)
        density_profile = (mass_profile / (4/3.*np.pi*radii**3) / (2.8*u.Da)).to(u.cm**-3)

        c_s = ((constants.k_B * tem_rprof*u.K / (2.4*u.Da))**0.5).to(u.km/u.s)
        azimuthal_average_MJ = (np.pi/6. * c_s**3 /
                                (constants.G**1.5 *
                                 (2.8*u.Da*azimuthal_average_density)**0.5)
                               ).to(u.M_sun)


        pl.figure(nplots*3+7)
        pl.plot(angular_radii, azimuthal_average_MJ, label=name)
        pl.ylabel("Azimuthally Averaged $M_J$ at $T(\\mathrm{CH}_3\\mathrm{OH})$")
        pl.xlabel("Radius (arcsec)")
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
                return ["%.1f" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius (au)")


        pl.figure(nplots*3+8)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(bins*pixscale*3600., mass_profile,
                label=name)
        pl.ylabel("Cumulative Mass at T(CH$_3$OH)")
        pl.xlabel("Radius (arcsec)")
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
            def tick_function(old_x):
                newx = (old_x*u.arcsec*masscalc.distance).to(u.au, u.dimensionless_angles()).value
                return ["%.1f" % z for z in newx]
            new_tick_locations = np.arange(0,8000,1000)*u.au
            new_tick_locs_as = (new_tick_locations/masscalc.distance).to(u.arcsec, u.dimensionless_angles())
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(new_tick_locs_as.value)
            ax2.set_xticklabels(tick_function(new_tick_locs_as.value))
            ax2.set_xlabel(r"Radius (au)")

        pl.figure(nplots*3+9)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(bins*pixscale*3600., tem_rprof,
                label=name)
        pl.legend(loc='best')
        pl.ylabel("CH$_3$OH temperature")
        pl.xlabel("Radius (as)")



regions = pyregion.open(paths.rpath("hmcore_centroids.reg"))
make_rprof(regions, ploteach=True)
nplots = len(regions)
pl.figure(nplots*3+2).savefig(paths.fpath("cumulative_radial_flux_massivecores.png"))
pl.figure(nplots*3+3).savefig(paths.fpath("cumulative_radial_mass_massivecores.png"))
pl.figure(nplots*3+4).savefig(paths.fpath("cumulative_density_40K_massivecores.png"))
pl.figure(nplots*3+5).savefig(paths.fpath("azimuthalaverage_density_40K_massivecores.png"))
pl.figure(nplots*3+6).savefig(paths.fpath("azimuthalaverage_radial_mj_40K_massivecores.png"))
pl.figure(nplots*3+7).savefig(paths.fpath("azimuthalaverage_radial_mj_of_TCH3OH_massivecores.png"))
pl.figure(nplots*3+8).savefig(paths.fpath("cumulative_radial_mass_of_TCH3OH_massivecores.png"))

#regions = pyregion.open(paths.rpath("cores.reg"))
#make_rprof(regions, ploteach=False)
#nplots = 3
#pl.figure(nplots*3+2).savefig(paths.fpath("cumulative_radial_flux_allcores.png"))
#pl.figure(nplots*3+3).savefig(paths.fpath("cumulative_radial_mass_allcores.png"))

pl.draw()
pl.show()
