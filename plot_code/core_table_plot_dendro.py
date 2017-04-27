import numpy as np
import paths
from astropy.table import Table, Column
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy import coordinates
import powerlaw
import pylab as pl
import scipy
from matplotlib.patches import Circle

# can take time:
from volume_integrals import mass_scalings

pl.style.use('classic')
pl.matplotlib.rc_file('pubfiguresrc')
pl.mpl.rcParams['font.size'] = 14.0
pl.mpl.rcParams['axes.titlesize'] = 16.0
pl.mpl.rcParams['axes.labelsize'] = 15.0
pl.mpl.rcParams['axes.color_cycle'] = ('338ADD', '9A44B6', 'A60628', '467821',
                                       'CF4457', '188487', 'E24A33', 'b', 'r',
                                       'g', 'm', 'k')

#pruned_ppcat = Table.read(paths.tpath("dendrogram_continuum_catalog.ipac"), format='ascii.ipac')
dendro_merge = Table.read(paths.tpath('dendro_merge_continuum_and_line.ipac'), format='ascii.ipac')
corelike = dendro_merge['corelike'] == 'True'

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

ax1.hist(dendro_merge['peak_cont_flux'], log=True, bins=np.logspace(-3,-0.5,15))
ax1.set_xscale('log')
ax1.set_ylim(0.3, 51)


fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(211)

fit = powerlaw.Fit(dendro_merge['peak_cont_flux'])
fit.plot_ccdf(color='k')
fit.power_law.plot_ccdf(color='r', linestyle='--')
ax2.set_ylabel("Fraction of sources")

ax3 = fig2.add_subplot(212)

fit = powerlaw.Fit(dendro_merge['peak_cont_flux'])
# doesn't work at all fit.plot_pdf(color='k')
ax3.hist(dendro_merge['peak_cont_flux'], bins=np.logspace(-3,-0.5,15),
         color='k', facecolor='none', histtype='step')
ax3.set_xscale('log')
fit.power_law.plot_pdf(color='r', linestyle='--')
ax3.set_ylim(0.3, 51)
ax2.set_xlim(ax3.get_xlim())
ax3.set_xlabel("Peak flux density (Jy/beam)")
ax3.set_ylabel("Number of sources")
fig2.savefig(paths.fpath('coreplots/dendro_flux_powerlaw_histogram_fit.png'), bbox_inches='tight')

print("Fit parameters: alpha={0}".format(fit.power_law.alpha))

radii = (0.2,0.4,0.6,0.8,1.0,1.5)*u.arcsec
lines = np.array([[row['peak_cont_flux']] + [row['cont_flux{0}arcsec'.format(rad).replace(".","p")]
                                        for rad in radii.value]
                  for row in dendro_merge])
pradii = (0.0,0.2,0.4,0.6,0.8,1.0,1.5)

pl.clf()
pl.plot(pradii, lines.T)
pl.plot(pradii, lines[corelike].T)
pl.xlabel("Aperture Radius (\")")
pl.ylabel("Flux (Jy)")
pl.savefig(paths.fpath('coreplots/flux_vs_aperture_radius_alldendrosources.png'), bbox_inches='tight')

pl.clf()
# 0.05" pixels
ppbeam = 16.201645578251686
background = (lines[:,2]-lines[:,1])/(np.pi*(pradii[2]-pradii[1])**2/0.05**2) * ppbeam
pl.plot(dendro_merge['peak_cont_flux'], background, '.')
pl.plot(dendro_merge['peak_cont_flux'][corelike], background[corelike], '.')
pl.xlabel("Peak Flux")
pl.ylabel("Background Flux")
pl.savefig(paths.fpath('coreplots/dendro_peak_vs_background.png'), bbox_inches='tight')


pl.clf()
pl.plot(dendro_merge['peak_cont_flux'], (dendro_merge['peak_cont_flux']-background), '.')
pl.plot(dendro_merge['peak_cont_flux'][corelike], (dendro_merge['peak_cont_flux']-background)[corelike], '.')
pl.xlabel("Peak Flux")
pl.ylabel("Background-subtracted Flux")
pl.savefig(paths.fpath('coreplots/dendro_peak_vs_peak_minus_background.png'), bbox_inches='tight')


pl.clf()
pl.title("Sanity check: what fraction of the peak flux is recovered in a 0.2\" aperture?")

# integrate over a gaussian to determine aperture correction
def gaussian2d(x,y):
    return np.exp(-x**2/2.-y**2/2.)
def gauss2d_integral(radius):
    from scipy.integrate import dblquad
    return dblquad(gaussian2d, 0, radius, lambda x: 0, lambda x:
                   (radius**2-x**2)**0.5)[0] * 4

aperture_correction = gauss2d_integral(0.2 / (dendro_merge['beam_area'][0]**0.5 * 206265.)) / gauss2d_integral(100)
pl.plot(dendro_merge['cont_flux0p2arcsec'], (dendro_merge['peak_cont_flux']), '.')
pl.plot(dendro_merge['cont_flux0p2arcsec'][corelike], (dendro_merge['peak_cont_flux'])[corelike], '.')
pl.plot([0.,0.5], [0,0.5], 'k--', zorder=-5)
pl.plot([0.,0.5], [0,0.5*aperture_correction], 'k:', zorder=-5)
pl.xlabel("Peak Flux")
pl.ylabel("0.2\" aperture flux")
pl.savefig(paths.fpath('coreplots/dendro_peak_vs_0p2arcsec.png'), bbox_inches='tight')



fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(211)

fit = powerlaw.Fit(dendro_merge['peak_cont_mass'])
fit.plot_ccdf(color='k')
fit.power_law.plot_ccdf(color='r', linestyle='--')
ax2.set_ylabel("Fraction of sources")

ax3 = fig2.add_subplot(212)

fit = powerlaw.Fit(dendro_merge['peak_cont_mass'])
# doesn't work at all fit.plot_pdf(color='k')
bmin, bmax = 0.1, 350
bins = np.logspace(np.log10(bmin),np.log10(bmax),20)
#bins = np.linspace((bmin),(bmax),15)
H,L,P = ax3.hist(dendro_merge['peak_cont_mass'], bins=bins, color='k',
                 facecolor='none', histtype='step')
pdf = fit.power_law.pdf(bins)/fit.power_law.pdf(bins).max()*np.max(H)
ax3.plot(bins[bins>fit.power_law.xmin], pdf, 'r--')
#fit.power_law.plot_pdf(color='r', linestyle='--')
#ax3.set_ylim(0.03, 0.5)
ax3.set_xscale('log')
ax2.set_xlim(ax3.get_xlim())
#ax3.set_yscale('log')
ax3.set_xlabel("Temperature-corrected mass")
#ax3.set_ylabel("Fraction of sources")
ax3.set_ylabel("Number of sources")
fig2.savefig(paths.fpath('coreplots/dendro_tcorr_mass_powerlaw_histogram_fit.png'), bbox_inches='tight')

print("Fit parameters: alpha={0}".format(fit.power_law.alpha))



powerlaw_parameters = {}

fig2 = pl.figure(2)
fig2.clf()
#ax = fig2.gca()

bins = np.logspace(-3,1.05, 20)

ax = pl.subplot(7,1,1)

H,L,P = ax.hist(dendro_merge['peak_cont_flux'], bins=bins, facecolor='none',
                histtype='step')
ax.set_xscale('log')
ax.set_ylim(0, H.max()+1)
ax.set_xlim(0.001, 15)
ax.set_ylabel("Peak")

ax.xaxis.set_ticklabels([])
ax.set_xlabel("")

max_yticks=3
yloc = pl.MaxNLocator(max_yticks)
ax.yaxis.set_major_locator(yloc)

apertures = ('0p2', '0p4', '0p6', '0p8', '1p0', '1p5')
for ii,aperture in enumerate(apertures):
    flux = (dendro_merge['cont_flux{0}arcsec'.format(aperture)] -
            dendro_merge['KUbandcont_flux{0}arcsec'.format(aperture)])
    freefreedominated = (dendro_merge['KUbandcont_flux{0}arcsec'.format(aperture)] /
          dendro_merge['cont_flux{0}arcsec'.format(aperture)]) > 0.5


    fit = powerlaw.Fit(flux[~freefreedominated])
    print("Powerlaw fit for apertures {0}: {1}+/-{4}     xmin: {2}"
          "    n: {3}".format(aperture, fit.power_law.alpha,
                              fit.power_law.xmin, fit.power_law.n,
                              fit.power_law.sigma,
                             ))
    powerlaw_parameters[aperture] = {'alpha':fit.power_law.alpha,
                                     'e_alpha':fit.power_law.sigma,
                                     'xmin':fit.power_law.xmin,
                                     'number':fit.power_law.n,
                                    }


    ax = pl.subplot(7,1,ii+1)
    ax.set_ylabel("${0}''$".format(aperture.replace("p",".")))
    H,L,P = ax.hist(flux[~freefreedominated], bins=bins,
                    #facecolor='none',
                    #alpha=0.5,
                    histtype='step',
                    label=aperture)
    if ii < 5:
        ax.set_xticklabels([])
        ax.get_xaxis().set_visible(False)

    max_yticks=3
    yloc = pl.MaxNLocator(max_yticks)
    ax.yaxis.set_major_locator(yloc)

    ax.set_xscale('log')
    ax.set_xlim(0.002, 15)
    ax.set_ylim(0, H.max()+1)
ax.set_xlabel("Flux (Jy)")
pl.subplots_adjust(hspace=0)
pl.savefig(paths.fpath("coreplots/dendro_core_flux_histogram_apertureradius.png"), bbox_inches='tight')
#pl.legend(loc='best')

powerlaw_par_table = Table()
for colname in ('alpha','e_alpha','xmin','number'):
    powerlaw_par_table.add_column(Column([powerlaw_parameters[ap][colname]
                                          for ap in apertures],
                                         name=colname))
powerlaw_par_table.write(paths.tpath("dendro_powerlaw_fit_parameters.ipac"),
                         format='ascii.ipac')





fig2 = pl.figure(2)
fig2.clf()
ax3 = fig2.add_subplot(111)
bmin, bmax = 0.2, 6.0
bins = np.linspace((bmin),(bmax),15)
H,L,P = ax3.hist(dendro_merge['peak_cont_mass'], bins=bins*0.99, color='k',
                 facecolor='none', histtype='step', label='M($20$ K)',
                 linewidth=2, alpha=0.5)
#H,L,P = ax3.hist(cores_merge['peak_cont_mass'], bins=np.linspace(bmin, 130, 50), color='b',
#                 facecolor='none', histtype='step', label='M($20$K)',
#                 linewidth=2, alpha=0.5)
#peak_plot = P
#starless = Table.read('/Users/adam/work/catalogs/enoch_perseus/table1.dat',
#                      format='ascii.cds',
#                      readme='/Users/adam/work/catalogs/enoch_perseus/ReadMe')
#protostellar = Table.read('/Users/adam/work/catalogs/enoch_perseus/table2.dat',
#                          format='ascii.cds',
#                          readme='/Users/adam/work/catalogs/enoch_perseus/ReadMe')
#H,L,P = ax3.hist(starless['TMass'], bins=bins*0.98, color='r', linestyle='dashed',
#                 facecolor='none', histtype='step', label='Perseus Starless')
#H,L,P = ax3.hist(protostellar['TMass'], bins=bins*1.01, color='g', linestyle='dashed',
#                 facecolor='none', histtype='step', label='Perseus Protostellar')
#ax3.set_xlabel("Mass")
#ax3.set_ylabel("Number of sources")
#pl.legend(loc='best')
#fig2.savefig(paths.fpath('coreplots/mass_histograms.png'), bbox_inches='tight')
#peak_plot[0].set_visible(False)
#H,L,P = ax3.hist(cores_merge['peak_cont_mass'], bins=bins, color='b',
#                 facecolor='none', histtype='step', label='M($20$K)',
#                 linewidth=2, alpha=0.5)
#ax3.set_xlim(0,7)
#fig2.savefig(paths.fpath('coreplots/mass_histograms_low.png'), bbox_inches='tight')
#
#
#
#
# beam area is for the continuum data
beam_area = u.Quantity(np.array(dendro_merge['beam_area']), u.sr)
jy_to_k = (1*u.Jy).to(u.K, u.brightness_temperature(beam_area,
                                                    220*u.GHz)).mean()

m20kconv = float(dendro_merge.meta['keywords']['mass_conversion_factor']['value'].split()[1])
def m20ktickfunc(x):
    return ["{0:0.0f}".format(y*m20kconv) for y in x]

fig4 = pl.figure(4)
fig4.clf()
ax5 = fig4.gca()
#for species in np.unique(dendro_merge['PeakLineSpecies']):
#    if species != 'NONE':
#        mask = species == dendro_merge['PeakLineSpecies']
#        ax5.plot(dendro_merge['cont_flux0p4arcsec'][mask]-dendro_merge['cont_flux0p2arcsec'][mask],
#                 dendro_merge['peak_cont_flux'][mask], 's', label=species)
ax5.plot(dendro_merge['cont_flux0p4arcsec']-dendro_merge['cont_flux0p2arcsec'],
         dendro_merge['peak_cont_flux'], 'ks', alpha=0.75, label='')
# R1=2 R0 -> V1/V0 = 7
ax5.plot(np.array([0, 1.5])*mass_scalings['2-1to1-0'][0][0], np.array([0, 1.5]), 'b--', label='Constant density', zorder=-10)
ax5.plot(np.array([0, 1.5])*mass_scalings['2-1to1-0'][1][0], np.array([0, 1.5]), 'b:', label='$\\rho\\propto R^{-1}$', zorder=-10)
ax5.plot(np.array([0, 1.5])*mass_scalings['2-1to1-0'][2][0], np.array([0, 1.5]), 'b-.', label='$\\rho\\propto R^{-2}$', zorder=-10)
ax5.plot(np.array([0, 1.5])*mass_scalings['2-1to1-0'][3][0], np.array([0, 1.5]), 'b-', alpha=0.5, label='$\\rho\\propto R^{-3}$', zorder=-10)
ax5.set_ylabel("Peak continuum flux density (Jy/beam)")
ax5.set_xlabel("Background $1000 \\rm{AU} < r < 2000 \\rm{AU}$ continuum flux density (Jy)")
pl.legend(loc='best', prop={'size':16})
ax5.set_xlim(0, 1.5)
ax5.set_ylim(0, 0.5)
ax5x = ax5.twiny()
ax5x.set_xticklabels(m20ktickfunc(ax5.get_xticks()))
ax5x.set_xlabel("Background Mass $M(20\\rm{K})$")
ax5y = ax5.twinx()
ax5y.set_yticklabels(m20ktickfunc(ax5.get_yticks()))
ax5y.set_ylabel("Peak Mass $M(20\\rm{K})$")
fig4.savefig(paths.fpath('coreplots/dendro_continuum_background_vs_peak.png'), bbox_inches='tight')
ax5.set_xlim(0, 0.4)
ax5.set_ylim(0, 0.2)
ax5x.set_xticklabels(m20ktickfunc(ax5.get_xticks()))
ax5y.set_yticklabels(m20ktickfunc(ax5.get_yticks()))
fig4.savefig(paths.fpath('coreplots/dendro_continuum_background_vs_peak_zoom.png'), bbox_inches='tight')
ax5.set_xlim(0, 0.25)
ax5.set_ylim(0, 0.06)
ax5x.set_xticklabels(m20ktickfunc(ax5.get_xticks()))
ax5y.set_yticklabels(m20ktickfunc(ax5.get_yticks()))
fig4.savefig(paths.fpath('coreplots/dendro_continuum_background_vs_peak_zoom_more.png'), bbox_inches='tight')

fig4 = pl.figure(4)
fig4.clf()
ax5 = fig4.gca()
ax5.plot(dendro_merge['cont_flux0p6arcsec']-dendro_merge['cont_flux0p4arcsec'],
         dendro_merge['cont_flux0p4arcsec']-dendro_merge['cont_flux0p2arcsec'], 'ks', alpha=0.75, label='')
# R1=2 R0 -> V1/V0 = 7
ax5.plot(np.array([0, 1.5])*mass_scalings['3-2to2-1'][0][0], np.array([0, 1.5]), 'b--', label='Constant density', zorder=-10)
ax5.plot(np.array([0, 1.5])*mass_scalings['3-2to2-1'][1][0], np.array([0, 1.5]), 'b:', label='$\\rho\\propto R^{-1}$', zorder=-10)
ax5.plot(np.array([0, 1.5])*mass_scalings['3-2to2-1'][2][0], np.array([0, 1.5]), 'b-.', label='$\\rho\\propto R^{-2}$', zorder=-10)
ax5.plot(np.array([0, 1.5])*mass_scalings['3-2to2-1'][3][0], np.array([0, 1.5]), 'b-', alpha=0.5, label='$\\rho\\propto R^{-3}$', zorder=-10)
ax5.set_ylabel("Background $1000 \\rm{AU} < r < 2000 \\rm{AU}$ continuum flux density (Jy)")
ax5.set_xlabel("Background $2000 \\rm{AU} < r < 3000 \\rm{AU}$ continuum flux density (Jy)")
pl.legend(loc='best', prop={'size':16})
ax5x = ax5.twiny()
ax5.set_ylim(0, 1.5)
ax5.set_xlim(0, 1.5)
ax5x.set_xticklabels(m20ktickfunc(ax5.get_xticks()))
ax5x.set_xlabel("Big Annulus mass $M(20\\rm{K})$")
ax5y = ax5.twinx()
ax5y.set_yticklabels(m20ktickfunc(ax5.get_yticks()))
ax5y.set_ylabel("Small Annulus Mass $M(20\\rm{K})$")
fig4.savefig(paths.fpath('coreplots/dendro_continuum_background_vs_peak_shell1to2.png'), bbox_inches='tight')
ax5.set_xlim(0, 0.4)
ax5.set_ylim(0, 0.4)
ax5x.set_xticklabels(m20ktickfunc(ax5.get_xticks()))
ax5y.set_yticklabels(m20ktickfunc(ax5.get_yticks()))
fig4.savefig(paths.fpath('coreplots/dendro_continuum_background_vs_peak_shell1to2_zoom.png'), bbox_inches='tight')


########## REQUIRES LINE INFO ###########
fig3 = pl.figure(3)
fig3.clf()
ax4 = fig3.gca()
for species in np.unique(dendro_merge['PeakLineSpecies']):
    if species != 'NONE':
        mask = species == dendro_merge['PeakLineSpecies']
        ax4.plot(dendro_merge['peak_cont_flux'][mask], dendro_merge['PeakLineBrightness'][mask], 's', label=species)
ax4.plot([0,0.4], [0, 0.4*jy_to_k.value], 'k--')
ax4.set_xlabel("Continuum flux density (Jy/beam)")
ax4.set_ylabel("Peak line brightness (K)")
ax4.set_xlim([0, 0.4])
pl.legend(loc='best', fontsize=14)
fig3.savefig(paths.fpath('coreplots/dendro_peakTB_vs_continuum.png'), bbox_inches='tight')

fig3 = pl.figure(3)
fig3.clf()
ax4 = fig3.gca()
for species in np.unique(dendro_merge['PeakLineSpecies']):
    if species != 'NONE':
        mask = species == dendro_merge['PeakLineSpecies']
        ax4.plot(dendro_merge['PeakLineContinuumBG'][mask], dendro_merge['PeakLineBrightness'][mask], 's', label=species)
ax4.plot([0,6], [0, 6], 'k--')
ax4.set_xlabel("Continuum Brightness (K)")
ax4.set_ylabel("Peak line brightness (K)")
ax4.set_xlim([0, 6])
pl.legend(loc='best', fontsize=14)
fig3.savefig(paths.fpath('coreplots/dendro_peakTB_vs_selfconsistentcontinuum.png'), bbox_inches='tight')


fig6 = pl.figure(6)
fig6.clf()

flux02 = (dendro_merge['cont_flux0p2arcsec'] -
          dendro_merge['KUbandcont_flux0p2arcsec'])
flux04 = (dendro_merge['cont_flux0p4arcsec'] -
          dendro_merge['KUbandcont_flux0p4arcsec'])
flux06 = (dendro_merge['cont_flux0p6arcsec'] -
          dendro_merge['KUbandcont_flux0p6arcsec'])

beam_fwhm = 0.2/(8*np.log(2))**0.5
# apparently the integral of a 2D gaussian is approximately
# int_-x^+x e^(-x^2/ (2*sigma**2)) = 2 pi sigma^2 * erf(x * pi / 2.)**2
# or, in terms of FWHM,
# int_-x^+x e^(-x^2/ (2*fwhm**2/2.35**2)) = 2 pi sigma^2 * erf(x * pi / 2.)**2
# if x is 0.2" and the fwhm is 0.2", x = 2.35 sigma
gaussian_0p2 = scipy.special.erf(0.2/beam_fwhm * (2./np.pi))**2
gaussian_0p4 = scipy.special.erf(0.4/beam_fwhm * (2./np.pi))**2
gaussian_0p6 = scipy.special.erf(0.6/beam_fwhm * (2./np.pi))**2

pl.plot(flux02/flux04, flux04/flux06, '.')
pl.vlines(gaussian_0p2/gaussian_0p4, 0, 1, color='k', linestyle='--')
pl.hlines(gaussian_0p4/gaussian_0p6, 0, 1, color='k', linestyle='--')
#pl.xlabel("0.2\" aperture flux")
pl.xlabel("0.2\" over 0.4\" concentration parameter")
pl.ylabel("0.4\" over 0.6\" concentration parameter")
pl.axis((0,1.1,0,1.1))
#pl.ylim(0,1.5)


fig4 = pl.figure(4)
fig4.clf()
ax5 = fig4.gca()
ax5.plot(dendro_merge['peak_cont_mass'], dendro_merge['T_corrected_mass'], 's')
ylims = ax5.get_ylim()
ax5.plot([0,20], [0,20], 'k--')
ax5.set_ylim(ylims)
ax5.set_xlabel("Mass at 20K [M$_\\odot$]")
ax5.set_ylabel("Mass at peak $T_B$ [M$_\\odot$]")
fig4.savefig(paths.fpath('coreplots/dendro_mass20K_vs_massTB.png'), bbox_inches='tight')





ffile = fits.open(paths.dpath('w51_te_continuum_best.fits'))
mywcs = wcs.WCS(ffile[0].header)
img = np.isfinite(ffile[0].data)
x,y = img.max(axis=0), img.max(axis=1)
xlim = np.where(x)[0].min(),np.where(x)[0].max()
ylim = np.where(y)[0].min(),np.where(y)[0].max()

fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(1,1,1, projection=mywcs)

coords = coordinates.SkyCoord(dendro_merge['x_cen'].data, dendro_merge['y_cen'].data,
                              unit=(u.deg, u.deg),
                              frame='fk5')
ax2.contour(img, levels=[0.5], colors='k')
sc = ax2.scatter(coords.ra.deg, coords.dec.deg, marker='.',
                 s=120,
                 transform=ax2.get_transform('fk5'), c=dendro_merge['mean_velo'],
                 edgecolors='none',
                 cmap=pl.cm.jet)
cb = pl.colorbar(mappable=sc, ax=ax2)
cb.set_label('Velocity (km s$^{-1}$)')
ax2.set_xlim(*xlim)
ax2.set_ylim(*ylim)
#ax2.set_xlim(*ax2.get_xlim()[::-1])
#ax2.set_ylabel('Galactic Latitude')
#ax2.set_xlabel('Galactic Longitude')
ax2.set_ylabel('Right Ascension')
ax2.set_xlabel('Declination')
ax2.set_aspect(1)

# bigger circle; doesn't fit south sources as well
#circlecen = coordinates.SkyCoord('19:23:40.747 +14:30:22.75', frame='fk5',
#                                 unit=(u.hour, u.deg))
#ax2.add_patch(Circle([circlecen.ra.deg, circlecen.dec.deg], radius=0.01287,
#                     facecolor='none', edgecolor='b', zorder=-10, alpha=0.5,
#                     linewidth=4,
#                     transform=ax2.get_transform('fk5')))

circlecen = coordinates.SkyCoord('19:23:41.253 +14:30:32.11', frame='fk5',
                                 unit=(u.hour, u.deg))
ax2.add_patch(Circle([circlecen.ra.deg, circlecen.dec.deg], radius=0.0105479,
                     facecolor='none', edgecolor='k', zorder=-15, alpha=0.2,
                     linewidth=8, linestyle=':',
                     transform=ax2.get_transform('fk5')))
fig2.savefig(paths.fpath('coreplots/dendro_core_spatial_distribution.png'), bbox_inches='tight')


pl.draw()
pl.show()
