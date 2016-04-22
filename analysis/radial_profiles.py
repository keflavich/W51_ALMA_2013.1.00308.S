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


regions = pyregion.open(paths.rpath("hmcore_centroids.reg"))
names = [r.attr[1]['text'] for r in regions]
center_positions = coordinates.SkyCoord([r.coord_list
                                         for r in regions],
                                        unit=(u.deg, u.deg),
                                        frame='fk5')
nplots = len(names)
for ii in range(nplots):
    pl.figure(ii).clf()
    pl.figure(nplots+ii).clf()

size = u.Quantity([1.25,1.25], u.arcsec)

linestyles = itertools.cycle(['-']*8 + ['--']*8 + [':']*8)

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

        linestyle = next(linestyles)

        pl.figure(ii)
        pl.title(name)
        pl.plot(bins*pixscale*3600., rprof*nr/ppbeam,
                label=fn.split(".")[0], linestyle=linestyle)
        pl.ylabel("Azimuthally Summed Flux")
        pl.xlabel("Radius (arcsec)")

        cumul_rprof = np.nan_to_num(rprof*nr/ppbeam).cumsum()

        pl.figure(nplots+ii)
        pl.title(name)
        pl.plot(bins*pixscale*3600., cumul_rprof,
                label=fn.split(".")[0], linestyle=linestyle)
        pl.ylabel("Cumulative Flux (Jy)")
        pl.xlabel("Radius (arcsec)")

for ii in range(len(names)):
    for xtra in (0,nplots):
        ax = pl.figure(ii+xtra).gca()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

pl.draw()
pl.show()
