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
from astropy import wcs
import image_tools
matplotlib.rc_file('pubfiguresrc')

for source in ('e8','north','e2',):

    matplotlib.pyplot.figure(1).clf()
    F1 = aplpy.FITSFigure(paths.dpath('chemslices/chemical_m0_slabs_{0}_CH3OH1029-936_merge.fits'.format(source)),
                          figure=matplotlib.pyplot.figure(1))
    F1.show_grayscale(invert=True, vmax=1100, vmin=-20)
    F1.show_contour(paths.dpath('W51_te_continuum_best.fits'),
                    colors=['r']*11,
                    levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
                            0.282361, 0.4, ])

    F1.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
    F1.scalebar.set_label('5000 au / 0.025 pc')
    F1.scalebar.set_color('k')

    F1.save(paths.fpath("{0}_continuum_contours_on_ch3oh1029.png".format(source)))

    matplotlib.pyplot.figure(2).clf()
    F = aplpy.FITSFigure(paths.dpath('W51_te_continuum_best.fits'),
                         figure=matplotlib.pyplot.figure(2))
    if source == 'e2': # TODO: mapping for each source
        F.recenter(290.93315, 14.5097, 2.8/3600.)
    F.show_grayscale(invert=True, vmax=0.43, vmin=-0.01, stretch='arcsinh')
    F.show_contour(paths.dpath('chemslices/chemical_m0_slabs_{0}_CH3OH1029-936_merge.fits'.format(source)),
                   colors=['b']*11,
                   levels=np.linspace(200,1100,6), layer='CH3OH')
    #F.ticks.set_xspacing(0.0005)

    F.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
    F.scalebar.set_label('5000 au / 0.025 pc')
    F.scalebar.set_color('k')

    F.save(paths.fpath("{0}_ch3oh1029_contours_on_continuum.png".format(source)))

    matplotlib.pyplot.figure(3).clf()
    F = aplpy.FITSFigure(paths.dpath('12m/moments/CH3OH_{0}_cutout_temperaturemap.fits'.format(source)),
                         figure=matplotlib.pyplot.figure(3))
    F.show_colorscale(vmax=600, vmin=50, cmap='hot')
    F.show_contour(paths.dpath('W51_te_continuum_best.fits'),
                   colors=['b']*11,
                   levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
                           0.282361, 0.4, ])

    F.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
    F.scalebar.set_label('5000 au / 0.025 pc')
    F.scalebar.set_color('k')

    F.save(paths.fpath("chemistry/{0}_continuum_contours_on_ch3oh_temperature.png".format(source)))

    matplotlib.pyplot.figure(4).clf()
    F = aplpy.FITSFigure(paths.dpath('12m/moments/CH3OH_{0}_cutout_columnmap.fits'.format(source)),
                         figure=matplotlib.pyplot.figure(4))
    F.show_colorscale(vmax=1e19, vmin=1e16, cmap='gray_r', stretch='arcsinh')
    #F.show_colorscale(vmax=1e21, vmin=1e15, cmap='viridis', stretch='log')
    F.show_contour(paths.dpath('W51_te_continuum_best.fits'),
                   colors=['r']*11,
                   levels=[0.015, 0.0256944, 0.0577778, 0.11125, 0.186111,
                           0.282361, 0.4, ])

    F.add_scalebar((0.025*u.pc / (5400*u.pc)).to(u.deg,u.dimensionless_angles()))
    F.scalebar.set_label('5000 au / 0.025 pc')
    F.scalebar.set_color('k')

    F.save(paths.fpath("chemistry/{0}_continuum_contours_on_ch3oh_temperature.png".format(source)))



    # compare CH3OH LTE column to dust emission
    ch3ohN_hdul = fits.open(paths.dpath('12m/moments/CH3OH_{0}_cutout_columnmap.fits'.format(source)))
    ch3ohT_hdul = fits.open(paths.dpath('12m/moments/CH3OH_{0}_cutout_temperaturemap.fits'.format(source)))
    dust_brightness,wts = reproject.reproject_interp(fits.open(paths.dpath('W51_te_continuum_best.fits')),
                                                     ch3ohN_hdul[0].header)

    w = wcs.WCS(ch3ohT_hdul[0].header)
    pixscale = (w.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg

    pl.figure(5).clf()
    ch3ohN = ch3ohN_hdul[0].data
    ch3ohT = ch3ohT_hdul[0].data
    pl.scatter(dust_brightness, ch3ohN, c=ch3ohT, vmax=600, vmin=50,
               edgecolor='none', alpha=0.9)
    pl.axis((1e-3, 0.2, 1e17, 1e19))
    pl.loglog()
    pl.xlabel("Dust Brightness (Jy/beam)")
    pl.ylabel("N(CH$_3$OH) [cm$^{-2}$]")
    cb = pl.colorbar()
    cb.set_label('CH$_3$OH-derived Temperature')
    pl.savefig(paths.fpath('chemistry/{0}_CH3OH_LTE_vs_dust_brightness.png'.format(source)),
               bbox_inches='tight')

    bm = radio_beam.Beam.from_fits_header(paths.dpath("W51_te_continuum_best.fits"))

    dust_column = dust_emissivity.dust.colofsnu(225*u.GHz, dust_brightness*u.Jy,
                                                beamomega=bm,
                                                temperature=ch3ohT*u.K)

    pl.figure(6).clf()
    pl.plot([1e23,1e25], [1e16, 1e18], 'k-', label='$X=10^{-7}$', zorder=-1)
    pl.plot([1e23,1e25], [1e17, 1e19], 'k--', label='$X=10^{-6}$', zorder=-2)
    pl.scatter(dust_column, ch3ohN, c=ch3ohT, vmax=600, vmin=50, edgecolor='none',
               alpha=0.6)
    pl.loglog()
    pl.axis((1e23, 1e25, 1e17, 1e19))
    pl.xlabel("Dust-derived N(H$_2$) [cm$^{-2}$]")
    pl.ylabel("N(CH$_3$OH) [cm$^{-2}$]")
    cb = pl.colorbar()
    cb.set_label('CH$_3$OH-derived Temperature')
    pl.legend(loc='upper left')
    pl.savefig(paths.fpath('chemistry/{0}_CH3OH_LTE_vs_dust_column.png'.format(source)),
               bbox_inches='tight')

    pl.figure(7).clf()
    ch3oh_abundance = ch3ohN / dust_column.value
    pl.imshow(ch3oh_abundance, vmin=5e-8, vmax=1e-5, cmap='gray_r', norm=matplotlib.colors.LogNorm())
    cb = pl.colorbar()
    cb.set_label("Dust-derived N(H$_2$) [cm$^{-2}$]")
    pl.xticks([])
    pl.yticks([])
    pl.plot([80,80+(1*u.arcsec/pixscale)], [5, 5], 'k-')
    pl.annotate("1\" = 5400 AU", ((160+(1*u.arcsec/pixscale))/2. - 5, 7.5),
                horizontalalignment='center')
    pl.savefig(paths.fpath('chemistry/{0}_CH3OH_LTE_abundance_map.png'.format(source)), bbox_inches='tight')

    yy,xx = np.indices(ch3ohN.shape)
    if source == 'north':
        center = [84.,38.]
    else:
        center = [ch3ohN.shape[0]/2., ch3ohN.shape[1]/2.]
    yyc = (yy-center[0])
    xxc = (xx-center[1])
    rr = (yyc**2 + xxc**2)**0.5
    rr_as = (rr*pixscale).to(u.arcsec)
    theta = np.arctan2(yyc,xxc)*u.rad

    # mask is an INCLUDE mask
    # we don't see anything at all for abundances below 1e-10
    mask = (ch3oh_abundance > 1e-10) & (ch3oh_abundance < 1e-5)
    if source == 'e2':
        mask = mask & (((theta > 15*u.deg) & (theta < 345*u.deg)) | (theta < -15*u.deg))
    mask = mask & (rr_as < 3*u.arcsec) # there's not really any valid data out of this radius
    mask = mask & np.isfinite(ch3oh_abundance)
    # exclude high-abundance, low-column regions: likely to be div-by-zero zones
    mask = mask & (~((ch3ohN < 1e18) & (ch3oh_abundance > 5e-6)))
    mask = mask & (~((dust_brightness<1e-2) & (ch3ohT > 500) & (ch3oh_abundance > 1e-6)))

    mask = mask & (~((ch3ohT > 250) &
                     (ch3ohN < 1e18) &
                     (rr_as>1.5*u.arcsec))
                   )# these are low-column,
    # high-temperature: they're very likely to be bad fits, since there is no
    # high-excitation data at these positions.  Selecting on large radii also helps
    # protect against parameter-based-selection bias

    mask_hdu = fits.PrimaryHDU(data=mask.astype('float'), header=ch3ohN_hdul[0].header)
    F1.hide_layer('contour_set_1')
    F1.show_contour(mask_hdu, levels=[0.5,1.5], colors=[(1,0,0,0.1)]*2, filled=True,
                    layer='mask_contours')
    F1.save(paths.fpath("chemistry/{0}_CH3OH_LTE_mask_contours_on_ch3oh1029.png".format(source)))


    #inds = np.argsort(rr.ravel())
    #bad = ~np.isfinite(ch3oh_abundance)
    #cumul_mean = (np.cumsum((np.nan_to_num(ch3oh_abundance)*mask).ravel()[inds]) /
    #              np.cumsum(np.arange(inds.size)*(~bad).ravel()*mask.ravel()))

    # use nan_to_num to turn nans to s/zeros/ones/ and ensure the profile continues
    # outward, but set the weights of these zero-points to be zero using the mask
    # (use 'ones' to avoid div-by-zero error; this is not obviously necessary)
    ch3oh_abundance_toprofile = ch3oh_abundance.copy()
    ch3oh_abundance_toprofile[~np.isfinite(ch3oh_abundance)] = 1.0
    # force the masked values to be positive.  Negative values cause bad things.
    ch3oh_abundance_toprofile[~mask] = 1.0

    radbins,radialprof = image_tools.radialprofile.azimuthalAverage(ch3oh_abundance_toprofile,
                                                                    weights=mask.astype('float'),
                                                                    returnradii=True,
                                                                    binsize=1,
                                                                    # to ensure self-consistency, use an identical center
                                                                    center=center,
                                                                    interpnan=True)

    pl.figure(8).clf()
    pl.scatter(rr_as.value[mask],
               ch3oh_abundance[mask],
               c=ch3ohT[mask], vmax=600, vmin=50, edgecolor='none', alpha=0.9)
    #pl.semilogy()
    pl.plot((radbins*pixscale).to(u.arcsec).value, radialprof, color='k', alpha=0.5, linewidth=2)
    pl.axis((0,3, #(rr.max()*pixscale).to(u.arcsec).value,
             5e-8, 5e-6))
    pl.xlabel("Separation from {0} (\")".format(source))
    pl.ylabel("$X$(CH$_3$OH)")
    pl.gca().ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
    cb = pl.colorbar()
    cb.set_label('CH$_3$OH-derived Temperature')
    pl.savefig(paths.fpath('chemistry/{0}_CH3OH_LTE_abundance_radial_profile.png'.format(source)),
               bbox_inches='tight')

    pl.figure(9).clf()
    pl.scatter(rr_as.value[mask],
               ch3ohT[mask],
               c=ch3oh_abundance[mask], vmax=5e-6, vmin=5e-8, edgecolor='none', alpha=0.8,
               norm=matplotlib.colors.LogNorm())
    #pl.scatter(rr_as.value[~mask],
    #           ch3ohT[~mask],
    #           c=ch3oh_abundance[~mask], vmax=5e-6, vmin=5e-8, edgecolor='k', alpha=0.5,
    #           norm=matplotlib.colors.LogNorm(), zorder=-1)
    #pl.semilogy()
    pl.axis((0,3,#(rr.max()*pixscale).to(u.arcsec).value,
             50,600))
    pl.xlabel("Separation from {0} (\")".format(source))
    cb = pl.colorbar()
    cb.set_label("$X$(CH$_3$OH)")
    pl.ylabel('CH$_3$OH-derived Temperature')
    pl.savefig(paths.fpath('chemistry/{0}_CH3OH_LTE_temperature_radial_profile.png'.format(source)),
               bbox_inches='tight')



    pl.figure(10).clf()
    pl.scatter(ch3ohT[mask],
               ch3ohN[mask],
               c=ch3oh_abundance[mask], vmax=5e-6, vmin=1e-7, edgecolor='none', alpha=0.8,
               norm=matplotlib.colors.LogNorm())
    pl.semilogy()
    #pl.axis((0,(rr.max()*pixscale).to(u.arcsec).value,
    #         50,600))
    pl.axis((50,600, 0.5e17,1e19))
    cb = pl.colorbar()
    cb.set_label("$X$(CH$_3$OH)")
    pl.xlabel('CH$_3$OH-derived Temperature')
    pl.ylabel('N(CH$_3$OH) [cm$^{-2}$]')
    pl.savefig(paths.fpath('chemistry/{0}_CH3OH_LTE_temperature_vs_column.png'.format(source)),
               bbox_inches='tight')


    pl.figure(11).clf()
    pl.scatter(ch3ohT[mask],
               ch3oh_abundance[mask],
               c=ch3ohN[mask], vmax=1e17, vmin=1e19, edgecolor='none', alpha=0.8,
               norm=matplotlib.colors.LogNorm())
    pl.semilogy()
    #pl.axis((0,(rr.max()*pixscale).to(u.arcsec).value,
    #         50,600))
    pl.axis((50,600, 5e-8, 5e-6,))
    cb = pl.colorbar()
    pl.ylabel("$X$(CH$_3$OH)")
    pl.xlabel('CH$_3$OH-derived Temperature')
    cb.set_label('N(CH$_3$OH) [cm$^{-2}$]')
    pl.savefig(paths.fpath('chemistry/{0}_CH3OH_LTE_temperature_vs_abundance.png'.format(source)),
               bbox_inches='tight')
