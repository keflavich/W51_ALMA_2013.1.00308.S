import numpy as np
import pylab as pl
import paths
from spectral_cube import SpectralCube
from astropy import coordinates
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy import log
import aplpy
from constants import distance

e1e2 = coordinates.SkyCoord(290.93268,14.508363,unit=('deg','deg'), frame='icrs')

def set_tight_ticks(F):
    F.tick_labels.set_yformat('dd:mm:ss.ss')
    F.tick_labels.set_xformat('hh:mm:ss.ss')
    F.ticks.set_xspacing(0.001)
    F.ticks.set_yspacing(0.001)
    #F.tick_labels.set_x_full_label_side('left')

#aplpy.make_rgb_cube( ('W51-CBAND-feathered.fits','W51-X-ABCD-S1.VTESS.VTC.DAVID-MEH.fits','W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'), 'W51_CXU_rgb' )

pl.close(1)
figure = pl.figure(1)
figure.clf()

# clean the header of junk axes
hdu = fits.open(paths.vpath('data/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'))
hdu[0].data = hdu[0].data.squeeze()
hdu[0].header = wcs.WCS(hdu[0].header).sub([wcs.WCSSUB_CELESTIAL]).to_header()

F = aplpy.FITSFigure(hdu,convention='calabretta',figure=figure)
#F = aplpy.FITSFigure(dpath+'W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.rot45.fits',convention='calabretta',figure=figure)
F.tick_labels.set_xformat('dd.dd')
F.tick_labels.set_yformat('dd.dd')
F.tick_labels.set_font(size=20)
F.axis_labels.set_font(size=20)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.05,invert=True)
#e1 = coordinates.ICRS(290.93263,14.50745,unit=('deg','deg'))
#F.recenter(e1.ra.value,e1.dec.value,width=1/60.,height=1/60.)
#F.recenter(290.92633,14.514769,radius=1.4/60.)
# IRS2:
F.recenter(290.91644,14.518939,radius=0.3/60.)


#F.add_scalebar(length=((0.5 * u.pc)/distance*u.radian).to(u.degree).value)
#F.scalebar.set_label('0.5 pc')
#F.scalebar.set_color('black')
#F.scalebar.set_linewidth(3)
#F.scalebar.set_font_size(20)
#
#pl.draw()
#pl.show()
F.add_scalebar(length=((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
#F.scalebar.set_length(((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.1 pc')
#F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)

F.show_contour(paths.dpath("w51_te_continuum_best.fits"),
               levels=[0.02, 0.04, 0.08, 0.16],
               colors=['r']*10,
               layer='almate_cont_ours')
set_tight_ticks(F)
F.save(paths.fpath("almate_contours_on_w51ku.png"))
F.recenter(290.91644,14.518939,radius=0.15/60.)
F.save(paths.fpath("almate_contours_on_w51ku_zoomier.png"))

F.remove_layer('almate_cont_ours')

F.recenter(290.91644,14.518939,radius=0.3/60.)
hdu = fits.open(paths.dpath("longbaseline/W51ncax.cont.image.pbcor.fits"))[0]
crop_hdu = fits.PrimaryHDU(data=hdu.data.squeeze()[1500:-1500, 1500:-1500],
                           header=wcs.WCS(hdu.header).celestial[1500:-1500, 1500:-1500].to_header())

F.show_contour(crop_hdu,
               levels=[0.001, 0.004, 0.008, 0.016],
               colors=['r']*10,
               layer='almalb_cont_ours')
F.save(paths.fpath("almalb_contours_on_w51ku.png"))

F.recenter(290.91644,14.518939,radius=0.15/60.)
F.save(paths.fpath("almalb_contours_on_w51ku_zoomier.png"))
F.recenter(290.9163,14.5182,width=0.13/60., height=0.08/60.)
F.save(paths.fpath("almalb_contours_on_w51ku_more_zoomier.png"))
