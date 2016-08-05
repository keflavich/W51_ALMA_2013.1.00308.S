import numpy as np
import aplpy
import paths
import matplotlib
import reproject
import radio_beam
import dust_emissivity
from astropy.io import fits
import pylab as pl
from astropy import units as u
matplotlib.rc_file('pubfiguresrc')


matplotlib.pyplot.figure(1).clf()
F = aplpy.FITSFigure(paths.dpath('chemslices/chemical_m0_slabs_e2_CH3OH1029-936_merge.fits'),
                     figure=matplotlib.pyplot.figure(1))
F.show_grayscale(invert=True, vmax=1100, vmin=-20)
F.show_contour(paths.dpath('W51_te_continuum_best.fits'),
               colors=['r']*11,
               levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
                       0.282361, 0.4, ])
F.save(paths.fpath("continuum_contours_on_ch3oh1029.png"))

matplotlib.pyplot.figure(2).clf()
F = aplpy.FITSFigure(paths.dpath('W51_te_continuum_best.fits'),
                     figure=matplotlib.pyplot.figure(2))
F.recenter(290.93315, 14.5097, 2.8/3600.)
F.show_grayscale(invert=True, vmax=0.43, vmin=-0.01, stretch='arcsinh')
F.show_contour(paths.dpath('chemslices/chemical_m0_slabs_e2_CH3OH1029-936_merge.fits'),
               colors=['b']*11,
               levels=np.linspace(200,1100,6), layer='CH3OH')
#F.ticks.set_xspacing(0.0005)
F.save(paths.fpath("ch3oh1029_contours_on_continuum.png"))

matplotlib.pyplot.figure(3).clf()
F = aplpy.FITSFigure(paths.dpath('12m/moments/CH3OH_e2_cutout_temperaturemap.fits'),
                     figure=matplotlib.pyplot.figure(3))
F.show_colorscale(vmax=600, vmin=0, cmap='hot')
F.show_contour(paths.dpath('W51_te_continuum_best.fits'),
               colors=['b']*11,
               levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
                       0.282361, 0.4, ])
F.save(paths.fpath("chemistry/continuum_contours_on_ch3oh_temperature.png"))

matplotlib.pyplot.figure(4).clf()
F = aplpy.FITSFigure(paths.dpath('12m/moments/CH3OH_e2_cutout_columnmap.fits'),
                     figure=matplotlib.pyplot.figure(4))
F.show_colorscale(vmax=1e19, vmin=1e16, cmap='gray_r', stretch='arcsinh')
#F.show_colorscale(vmax=1e21, vmin=1e15, cmap='viridis', stretch='log')
F.show_contour(paths.dpath('W51_te_continuum_best.fits'),
               colors=['r']*11,
               levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
                       0.282361, 0.4, ])
F.save(paths.fpath("chemistry/continuum_contours_on_ch3oh_temperature.png"))



# compare CH3OH LTE column to dust emission
ch3ohN_hdul = fits.open(paths.dpath('12m/moments/CH3OH_e2_cutout_columnmap.fits'))
ch3ohT_hdul = fits.open(paths.dpath('12m/moments/CH3OH_e2_cutout_temperaturemap.fits'))
dust_brightness,wts = reproject.reproject_interp(fits.open(paths.dpath('W51_te_continuum_best.fits')),
                                                 ch3ohN_hdul[0].header)

w = wcs.WCS(ch3ohT_hdul[0].header)
pixscale = (w.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg

pl.figure(5).clf()
ch3ohN = ch3ohN_hdul[0].data
ch3ohT = ch3ohT_hdul[0].data
pl.scatter(dust_brightness, ch3ohN, c=ch3ohT, vmax=600, vmin=0,
           edgecolor='none', alpha=0.9)
pl.axis((1e-3, 0.2, 1e17, 1e19))
pl.loglog()
pl.xlabel("Dust Brightness (Jy/beam)")
pl.ylabel("N(CH$_3$OH) [cm$^{-2}$]")
cb = pl.colorbar()
cb.set_label('CH$_3$OH-derived Temperature')
pl.savefig(paths.fpath('chemistry/CH3OH_LTE_vs_dust_brightness.png'))

bm = radio_beam.Beam.from_fits_header(paths.dpath("W51_te_continuum_best.fits"))

dust_column = dust_emissivity.dust.colofsnu(225*u.GHz, dust_brightness*u.Jy,
                                            beamomega=bm,
                                            temperature=ch3ohT*u.K)

pl.figure(6).clf()
pl.plot([1e23,1e25], [1e16, 1e18], 'k-', label='$X=10^{-7}$', zorder=-1)
pl.plot([1e23,1e25], [1e17, 1e19], 'k--', label='$X=10^{-6}$', zorder=-2)
pl.scatter(dust_column, ch3ohN, c=ch3ohT, vmax=600, vmin=0, edgecolor='none',
           alpha=0.9)
pl.loglog()
pl.axis((1e23, 1e25, 1e17, 1e19))
pl.xlabel("Dust-derived N(H$_2$) [cm$^{-2}$]")
pl.ylabel("N(CH$_3$OH) [cm$^{-2}$]")
cb = pl.colorbar()
cb.set_label('CH$_3$OH-derived Temperature')
pl.legend(loc='lower right')
pl.savefig(paths.fpath('chemistry/CH3OH_LTE_vs_dust_column.png'))

pl.figure(7).clf()
ch3oh_abundance = ch3ohN / dust_column.value 
pl.imshow(ch3oh_abundance, vmin=1e-7, vmax=1e-5, cmap='gray_r', norm=matplotlib.colors.LogNorm())
cb = pl.colorbar()
cb.set_label("Dust-derived N(H$_2$) [cm$^{-2}$]")
pl.xticks([])
pl.yticks([])
pl.plot([80,80+(1*u.arcsec/pixscale)], [5, 5], 'k-')
pl.annotate("1\" = 5400 AU", ((160+(1*u.arcsec/pixscale))/2. - 5, 7.5),
            horizontalalignment='center')
pl.savefig(paths.fpath('chemistry/CH3OH_LTE_abundance_map.png'))

pl.figure(8).clf()
yy,xx = np.indices(ch3ohN.shape)
yyc = (yy-ch3ohN.shape[0]/2.)
xxc = (xx-ch3ohN.shape[1]/2.)
rr = (yyc**2 + xxc**2)**0.5
theta = np.arctan2(yyc,xxc)*u.rad
mask = ((theta > 15*u.deg) & (theta < 345*u.deg)) | (theta < -15*u.deg)

pl.scatter((rr*pixscale).to(u.arcsec).value[mask],
           ch3oh_abundance[mask],
           c=ch3ohT[mask], vmax=600, vmin=0, edgecolor='none', alpha=0.9)
#pl.semilogy()
pl.axis((0,(rr.max()*pixscale).to(u.arcsec).value,
         1e-7, 5e-6))
pl.xlabel("Separation from e2e (\")")
pl.ylabel("$X$(CH$_3$OH)")
cb = pl.colorbar()
cb.set_label('CH$_3$OH-derived Temperature')
pl.savefig(paths.fpath('chemistry/CH3OH_LTE_abundance_radial_profile.png'))
