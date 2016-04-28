import numpy as np
import radio_beam
import pyregion
import paths
import image_tools
from astropy import wcs
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import units as u
from astropy import coordinates
import pylab as pl
import itertools
import masscalc

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
            print("fn  {0} ppbeam={1:0.2f}".format(fn, ppbeam))
            
            for ii,(name,position) in enumerate(zip(names, center_positions)):
                cutout = Cutout2D(fh[0].data, position, size, wcs=mywcs)

                nr, bins, rprof = image_tools.radialprofile.azimuthalAverage(cutout.data,
                                                                             binsize=1.0,
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
    pl.figure(nplots*3+1).clf()
    pl.figure(nplots*3+2).clf()
    pl.figure(nplots*3+3).clf()
    fn = "selfcal_allspw_selfcal_3_mfs_deeper.image.pbcor.fits"
    fh = fits.open(paths.dpath("12m/continuum/"+fn))
    mywcs = wcs.WCS(fh[0].header)
    beam = radio_beam.Beam.from_fits_header(fh[0].header)
    pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5
    ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam
    for ii,(name,position) in enumerate(zip(names, center_positions)):
        cutout = Cutout2D(fh[0].data, position, size, wcs=mywcs)

        nr, bins, rprof = image_tools.radialprofile.azimuthalAverage(cutout.data,
                                                                     binsize=1.0,
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
            yticks_Jy = yticks_mass/masscalc.mass_conversion_factor().value
            ax3.set_yticks(yticks_Jy)
            ax3.set_yticklabels(yticks_mass)
            ax3.set_ylabel("Cumulative Mass (M$_\\odot$, $T=20$ K)")

        pl.figure(nplots*3+3)
        #pl.title(fn.replace(".image.pbcor.fits",""))
        pl.plot(((bins*pixscale*u.deg)*masscalc.distance).to(u.pc,
                                                             u.dimensionless_angles()),
                cumul_rprof * masscalc.mass_conversion_factor(),
                label=name)
        pl.ylabel("Cumulative Mass (M$_\\odot$, $T=20$ K)")
        pl.xlabel("Radius (pc)")
        if len(names) < 5:
            pl.legend(loc='best')
        elif ii==len(names)-1:
            ax = pl.gca()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

            # Put a legend to the right of the current axis
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


regions = pyregion.open(paths.rpath("hmcore_centroids.reg"))
make_rprof(regions, ploteach=True)
nplots = len(regions)
pl.figure(nplots*3+2).savefig(paths.fpath("cumulative_radial_flux_massivecores.png"))
pl.figure(nplots*3+3).savefig(paths.fpath("cumulative_radial_mass_massivecores.png"))

regions = pyregion.open(paths.rpath("cores.reg"))
make_rprof(regions, ploteach=False)
nplots = 3
pl.figure(nplots*3+2).savefig(paths.fpath("cumulative_radial_flux_allcores.png"))
pl.figure(nplots*3+3).savefig(paths.fpath("cumulative_radial_mass_allcores.png"))

pl.draw()
pl.show()
