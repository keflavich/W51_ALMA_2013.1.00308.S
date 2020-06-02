import numpy as np
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling, time; from astropy.io import fits, ascii
fh = fits.open('/Users/adam/work/students/JoshMachado2019/W51/data/par_maps.fits')
import pylab as pl
pl.rcParams['font.size'] = 16

ww = wcs.WCS(fh[0].header)
w51e2=coordinates.SkyCoord.from_name('W51 e2')
xpix, ypix = ww.celestial.wcs_world2pix(w51e2.fk5.ra, w51e2.fk5.dec, 0)
yy,xx = np.indices(fh[0].data.shape[1:])
rr2 = (xx-xpix)**2 + (yy-ypix)**2
pixscale = wcs.utils.proj_plane_pixel_scales(ww.celestial)[0]*u.deg
rr = ((rr2**0.5)*pixscale).to(u.arcsec)
tem = fh[0].data[0,:,:]



ch3oht = fits.open('/Users/adam/work/w51/alma/FITS/12m/moments/CH3OH_e2_cutout_temperaturemap.fits')

wwc = wcs.WCS(ch3oht[0].header)
xpixc, ypixc = wwc.celestial.wcs_world2pix(w51e2.fk5.ra, w51e2.fk5.dec, 0)
yyc,xxc = np.indices(ch3oht[0].data.shape)
rr2c = (xxc-xpixc)**2 + (yyc-ypixc)**2
pixscalec = wcs.utils.proj_plane_pixel_scales(wwc.celestial)[0]*u.deg
rrc = ((rr2c**0.5)*pixscalec).to(u.arcsec)
okc = np.isfinite(ch3oht[0].data) & (ch3oht[0].data > 5)

pl.clf()

ok = np.isfinite(tem) & (rr<30*u.arcsec) & (tem > 5)
pl.plot(rr[ok] * 5400/206265*u.pc/u.arcsec, tem[ok], '.', label='NH$_3$')

pl.plot(rrc[okc] * 5400/206265*u.pc/u.arcsec, ch3oht[0].data[okc], '.', label='CH$_3$OH')
pl.ylim(1, 600)
pl.xlim(0, 0.7)

xax = np.linspace(0.02, 1)
pl.plot(xax, 500/(xax/xax[0]), label='1/r', color='k', linestyle='--')
#pl.plot(xax, 500/(xax/xax[0])**2, label='1/r$^2$')
pl.xlabel("Radius (pc)")
pl.ylabel("Temperature (K)")
pl.legend(loc='best')
