import numpy as np
import paths
from astropy.table import Table, Column
from astropy import units as u
from astropy import coordinates
import powerlaw
import pylab as pl
from astropy.io import fits
from astropy import wcs
from matplotlib.patches import Circle
import matplotlib
import os
from astropy import log
import contextlib
devnull = open(os.devnull,'w')

log.setLevel('CRITICAL')

pl.style.use('classic')
pl.matplotlib.rc_file('pubfiguresrc')
pl.mpl.rcParams['font.size'] = 14.0
pl.mpl.rcParams['axes.titlesize'] = 16.0
pl.mpl.rcParams['axes.labelsize'] = 15.0
pl.mpl.rcParams['axes.color_cycle'] = ('338ADD', '9A44B6', 'A60628', '467821',
                                       'CF4457', '188487', 'E24A33', 'b', 'r',
                                       'g', 'm', 'k')

core_velo_tbl = Table.read(paths.tpath("core_velocities.ipac"), format="ascii.ipac")
core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
cores_merge = Table.read(paths.tpath('core_continuum_and_line.ipac'), format='ascii.ipac')

beam_area = cores_merge['beam_area']
jy_to_k = (1*u.Jy).to(u.K, u.brightness_temperature(u.Quantity(beam_area, u.sr),
                                                    226*u.GHz)).mean()

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

ax1.hist(core_phot_tbl['peak'], log=True, bins=np.logspace(-3,-0.5,15))
ax1.set_xscale('log')
ax1.set_xlabel("Peak flux density (Jy)")
ax1.set_ylabel("Source Counts")
ax1.set_ylim(0.3, 11)


fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(211)

print()
print("Core distribution: peak")
with contextlib.redirect_stderr(devnull):
    fit = powerlaw.Fit(core_phot_tbl['peak'])
fit.plot_ccdf(color='k')
fit.power_law.plot_ccdf(color='r', linestyle='--')
ax2.set_ylabel("Fraction of sources")
ax2.xaxis.set_label_position('top')
ax2.xaxis.set_ticks_position('top')
def my_formatter_fun(x, p):
    brightness = x*jy_to_k
    return "$%0.1f$" % (brightness.value)
ax2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(my_formatter_fun))
ax2.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(
    np.array([1,3,5,10,20,40,100,200,300])/jy_to_k.value))
ax2.set_xlabel("Peak $T_B$ [K]")


ax3 = fig2.add_subplot(212)

with contextlib.redirect_stderr(devnull):
    fit = powerlaw.Fit(core_phot_tbl['peak'])
# doesn't work at all fit.plot_pdf(color='k')
ax3.hist(core_phot_tbl['peak'], bins=np.logspace(-3,-0.5,15),
         color='k', facecolor='none', histtype='step')
ax3.set_xscale('log')
fit.power_law.plot_pdf(color='r', linestyle='--')
ax3.set_ylim(0.3, 15)
ax3.set_xlabel("Peak flux density (Jy/beam)")
ax3.set_ylabel("Number of sources")
fig2.savefig(paths.fpath('coreplots/flux_powerlaw_histogram_fit.png'), bbox_inches='tight')

print("Flux Fit parameters: alpha={0}".format(fit.power_law.alpha))

fig = pl.figure(1)
fig.clf()
ax = fig.gca()
ax.hist(core_phot_tbl['peak'] * jy_to_k, histtype='step', facecolor='none', color='k',
        bins=np.linspace(0,100,26))
ax.set_xlabel("Peak $T_B$ [K]")
fig.savefig(paths.fpath('coreplots/peaktb_histogram.png'), bbox_inches='tight')


fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(211)

print()
print("Core distribution: tem-corrected peak mass")
with contextlib.redirect_stderr(devnull):
    fit = powerlaw.Fit(cores_merge['T_corrected_peakmass'])
print("xmin = {0}, alpha = {1}".format(fit.xmin, fit.alpha))
fit.plot_ccdf(color='k')
fit.power_law.plot_ccdf(color='r', linestyle='--')
ax2.set_ylabel("Fraction of sources")

ax3 = fig2.add_subplot(212)

with contextlib.redirect_stderr(devnull):
    fit = powerlaw.Fit(cores_merge['T_corrected_peakmass'])
# doesn't work at all fit.plot_pdf(color='k')
bmin, bmax = 0.2, 6.0
bmin, bmax = 1.0, 70
bins = np.logspace(np.log10(bmin),np.log10(bmax),21)
#bins = np.linspace((bmin),(bmax),15)
H,L,P = ax3.hist(cores_merge['T_corrected_peakmass'], bins=bins, color='k',
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
fig2.savefig(paths.fpath('coreplots/tcorr_mass_powerlaw_histogram_fit.png'), bbox_inches='tight')

print("Mass Fit parameters: alpha={0}".format(fit.power_law.alpha))

fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(211)

print()
print("Core distribution: 'masses to fit', which replaces <20K w/20K masses")
masses_to_fit = np.where(cores_merge['T_corrected_aperturemass'] < cores_merge['ApertureMass20K'],
                         cores_merge['T_corrected_aperturemass'], cores_merge['ApertureMass20K'])

with contextlib.redirect_stderr(devnull):
    fit = powerlaw.Fit(masses_to_fit)
print("xmin = {0}, alpha = {1}".format(fit.xmin, fit.alpha))
fit.plot_ccdf(color='k')
fit.power_law.plot_ccdf(color='r', linestyle='--')
ax2.set_ylabel("Fraction of sources")

ax3 = fig2.add_subplot(212)

with contextlib.redirect_stderr(devnull):
    fit = powerlaw.Fit(masses_to_fit)
# doesn't work at all fit.plot_pdf(color='k')
bmin, bmax = 0.1, 300.0
bins = np.logspace(np.log10(bmin),np.log10(bmax),15)
#bins = np.linspace((bmin),(bmax),15)
H,L,P = ax3.hist(masses_to_fit, bins=bins, color='k',
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
print("Aperture mass fit parameters: alpha={0}+/-{1}".format(fit.power_law.alpha, fit.power_law.sigma))
fig2.savefig(paths.fpath('coreplots/tcorr_aperture_mass_powerlaw_histogram_fit.png'), bbox_inches='tight')







fig2 = pl.figure(2)
fig2.clf()
ax3 = fig2.add_subplot(111)
bmin, bmax = 0.2, 6.0
bins = np.linspace((bmin),(bmax),15)
H,L,P = ax3.hist(cores_merge['T_corrected_peakmass'], bins=bins*0.99, color='k',
                 facecolor='none', histtype='step', label='M($T_B$)',
                 linewidth=2, alpha=0.5)
H,L,P = ax3.hist(cores_merge['peak_mass'], bins=np.linspace(bmin, 130, 50), color='b',
                 facecolor='none', histtype='step', label='M($20$K)',
                 linewidth=2, alpha=0.5)
peak_plot = P
starless = Table.read('/Users/adam/work/catalogs/enoch_perseus/table1.dat',
                      format='ascii.cds',
                      readme='/Users/adam/work/catalogs/enoch_perseus/ReadMe')
protostellar = Table.read('/Users/adam/work/catalogs/enoch_perseus/table2.dat',
                          format='ascii.cds',
                          readme='/Users/adam/work/catalogs/enoch_perseus/ReadMe')
H,L,P = ax3.hist(starless['TMass'], bins=bins*0.98, color='r', linestyle='dashed',
                 facecolor='none', histtype='step', label='Perseus Starless')
H,L,P = ax3.hist(protostellar['TMass'], bins=bins*1.01, color='g', linestyle='dashed',
                 facecolor='none', histtype='step', label='Perseus Protostellar')
ax3.set_xlabel("Mass")
ax3.set_ylabel("Number of sources")
pl.legend(loc='best')
fig2.savefig(paths.fpath('coreplots/mass_histograms.png'), bbox_inches='tight')
peak_plot[0].set_visible(False)
H,L,P = ax3.hist(cores_merge['peak_mass'], bins=bins, color='b',
                 facecolor='none', histtype='step', label='M($20$K)',
                 linewidth=2, alpha=0.5)
ax3.set_xlim(0,7)
fig2.savefig(paths.fpath('coreplots/mass_histograms_low.png'), bbox_inches='tight')





fig3 = pl.figure(3)
fig3.clf()
ax4 = fig3.gca()
ax4.plot(cores_merge['peak'], cores_merge['PeakLineBrightness'], 's')
ax4.plot([0,0.4], [0, 0.4*jy_to_k.value], 'k--')
ax4.set_xlabel("Continuum flux density (Jy/beam)")
ax4.set_ylabel("Peak line brightness (K)")
ax4.set_xlim([0, 0.4])
fig3.savefig(paths.fpath('coreplots/peakTB_vs_continuum.png'), bbox_inches='tight')

fig4 = pl.figure(4)
fig4.clf()
ax5 = fig4.gca()
ax5.plot(cores_merge['peak_mass'], cores_merge['T_corrected_peakmass'], 's')
ylims = ax5.get_ylim()
ax5.plot([0,250], [0,250], 'k--')
ax5.set_ylim(ylims)
ax5.set_xlabel("Mass at 20K [M$_\\odot$]")
ax5.set_ylabel("Mass at peak $T_B$ [M$_\\odot$]")
fig4.savefig(paths.fpath('coreplots/mass20K_vs_massTB.png'), bbox_inches='tight')



ffile = fits.open(paths.dpath('w51_te_continuum_best.fits'))
mywcs = wcs.WCS(ffile[0].header)
img = np.isfinite(ffile[0].data)
x,y = img.max(axis=0), img.max(axis=1)
xlim = np.where(x)[0].min(),np.where(x)[0].max()
ylim = np.where(y)[0].min(),np.where(y)[0].max()

ppbeam = cores_merge['beam_area'].to(u.deg**2) / np.abs(np.product(mywcs.pixel_scale_matrix.diagonal())*u.deg**2)
cores_merge['integrated'] = cores_merge['sum']/ppbeam



fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.gca()

print()
print("Core distribution: integrated flux")
with contextlib.redirect_stderr(devnull):
    fit = powerlaw.Fit(cores_merge['integrated'])
fit.plot_ccdf(color='k')
fit.power_law.plot_ccdf(color='r', linestyle='--')
ax2.set_ylabel("Fraction of sources")
ax2.xaxis.set_label_position('top')
ax2.xaxis.set_ticks_position('top')
ax2.set_xlabel("Integrated Flux (Jy)")

fig2.savefig(paths.fpath('coreplots/integrated_flux_powerlaw_histogram_fit.png'), bbox_inches='tight')

print("Integrated Flux Fit parameters: alpha={0}".format(fit.power_law.alpha))





fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(1,1,1, projection=mywcs)
ax2.imshow(ffile[0].data, cmap='gray_r', vmax=0.025, vmin=0.0001)


not_freefree = ~(cores_merge['is_freefree'].astype('bool'))

coords = coordinates.SkyCoord(cores_merge['RA'], cores_merge['Dec'],
                              frame='fk5')[not_freefree]
ax2.contour(img, levels=[0.5], colors='k')
sc = ax2.scatter(coords.ra.deg, coords.dec.deg, marker='.',
                 s=30*cores_merge['T_corrected_peakmass'][not_freefree],
                 transform=ax2.get_transform('fk5'), c=cores_merge['mean_velo'][not_freefree],
                 edgecolors='k',
                 alpha=0.5,
                 cmap=pl.cm.jet)
cb = pl.colorbar(mappable=sc, ax=ax2)
cb.set_label('Velocity (km s$^{-1}$)')
ax2.set_xlim(*xlim)
ax2.set_ylim(*ylim)
#ax2.set_xlim(*ax2.get_xlim()[::-1])
ax2.set_xlabel('Right Ascension')
ax2.set_ylabel('Declination')
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
fig2.savefig(paths.fpath('coreplots/core_spatial_distribution.png'), bbox_inches='tight')



powerlaw_parameters = {}

fig5 = pl.figure(5)
fig5.clf()
#ax = fig5.gca()

starless = cores_merge['Classification'] == 'StarlessCore'

bins = np.logspace(-3,1.05, 20)

ax = pl.subplot(7,1,1)

H,L,P = ax.hist(cores_merge['peak'][starless], bins=bins, facecolor='none',
                histtype='stepfilled', hatch='/')
H,L,P = ax.hist(cores_merge['peak'], bins=bins, facecolor='none',
                histtype='step')
ax.set_xscale('log')
ax.set_ylim(0, H.max()+1)
ax.set_xlim(0.001, 15)
ax.set_ylabel("Peak")

ax.xaxis.set_ticklabels([])
ax.set_xlabel("Peak $T_B$ [K]")
ax.set_xlabel("")

max_yticks=3
yloc = pl.MaxNLocator(max_yticks)
ax.yaxis.set_major_locator(yloc)


class_colors = {
    'DustyHII': 'w',
    'ExtendedColdCore': 'b',
    'ExtendedHotCore': 'k',
    #'HII': 'w',
    'HotCore': 'g',
    'StarlessCore': 'r',
    'UncertainCompact': 'c',
    'UncertainExtended': 'y',
}

print()
print("Fit each aperture excluding free-free fluxes")
apertures = ('0p2', '0p4', '0p6', '0p8', '1p0', '1p5',)
for ii,aperture in enumerate(apertures):
    flux = (cores_merge['cont_flux{0}arcsec'.format(aperture)] -
            cores_merge['KUbandcont_flux{0}arcsec'.format(aperture)])
    ff = (cores_merge['KUbandcont_flux{0}arcsec'.format(aperture)] /
          cores_merge['cont_flux{0}arcsec'.format(aperture)]) > 0.5


    with contextlib.redirect_stderr(devnull):
        fit = powerlaw.Fit(flux[~ff])
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


    ax = pl.subplot(7,1,ii+2)
    ax.set_ylabel("${0}''$".format(aperture.replace("p",".")))

    #nice idea, but it looks awful
    # histdata = [(cls, flux[cores_merge['Classification'] == cls])
    #             for cls in class_colors]

    # H,L,P = ax.hist([x[1] for x in histdata],
    #                 label=[x[0] for x in histdata],
    #                 color=[class_colors[x[0]] for x in histdata],
    #                 stacked=True)

    # do starless first b/c H.max() is used to set ylim
    H,L,P = ax.hist(flux[starless & ~ff], bins=bins,
                    facecolor='none',
                    #alpha=0.5,
                    histtype='stepfilled',
                    hatch='/',
                    label=aperture)


    H,L,P = ax.hist(flux[~ff], bins=bins,
                    #facecolor='none',
                    #alpha=0.5,
                    histtype='step',
                    label=aperture)


    if ii < 5:
        ax.set_xticklabels([])
        ax.get_xaxis().set_visible(False)

    yloc = pl.MaxNLocator(max_yticks)
    ax.yaxis.set_major_locator(yloc)

    ax.set_xscale('log')
    ax.set_xlim(0.001, 15)
    ax.set_ylim(0, H.max()+1)
ax.set_xlabel("Integrated or Peak Flux Density (Jy)")
pl.subplots_adjust(hspace=0)
pl.savefig(paths.fpath("coreplots/core_flux_histogram_apertureradius.png"), bbox_inches='tight')


pl.draw()
pl.show()


fig3 = pl.figure(3)
fig3.clf()

flux02 = (cores_merge['cont_flux0p2arcsec'.format(aperture)] -
          cores_merge['KUbandcont_flux0p2arcsec'.format(aperture)])
flux04 = (cores_merge['cont_flux0p4arcsec'.format(aperture)] -
          cores_merge['KUbandcont_flux0p4arcsec'.format(aperture)])

pl.plot(flux02, flux02/flux04, '.')
pl.xlabel("0.2\" aperture flux")
pl.ylabel("0.2\" over 0.4\" concentration parameter")
pl.ylim(0,1.5)


fig4 = pl.figure(4)
fig4.clf()
ax=fig4.gca()

isOK = np.isfinite(cores_merge['BrightestFittedPeakPixBrightness'])
histdata = []

# need to add back in continuum because we're concerned with the *absolute*
# brightness.  MOVED TO merge_spectral...py
#jtok_eq = u.brightness_temperature(cores_merge['beam_area'], 225*u.GHz)
#cont_brightness = (u.beam * cores_merge['sum']/cores_merge['npix']).to(u.K, jtok_eq)
#contincluded_line_brightness = cores_merge['BrightestFittedPeakPixBrightness'] + cont_brightness
cont_brightness = cores_merge['MeanContinuumBrightness']
contincluded_line_brightness = cores_merge['BrightestFittedPeakPixBrightnessWithcont']

linename_mapping = {'13CSv=0': '$^{13}$CS',
                    'CH3CH2CNv=0':'CH$_3$CH$_2$CN',
                    'CH3OCH3':'CH$_3$OCH$_3$',
                    'CH3OCHOv=0':'CH$_3$OCHO',
                    'CH3OHvt=0':'CH$_3$OH',
                    'HC3Nv=0':'HC$_3$N',
                    'NH2CHO':'NH$_2$CHO',
                    'SO3&Sigma;v=0':'SO',
                   }

for linename in np.unique(cores_merge['BrightestFittedLine']):
    if linename == '-':
        continue
    match = linename==cores_merge['BrightestFittedLine']
    histdata.append((linename_mapping[linename],
                     contincluded_line_brightness[match & isOK]))

histd, bins, patches = ax.hist([x[1] for x in histdata], label=[x[0] for x in
                                                                histdata],
                               stacked=True)
for pp, hh in zip(patches, ('//', '\\\\', '--', '||')*4):
    for art in pp:
        fc = art.get_facecolor()
        art.set_edgecolor(fc)
        fc = fc[:3] + (0.5,)
        art.set_facecolor(fc)
        art.set_hatch(hh)
        #art.set_linewidth(2)
#stack_hist(ax=ax, stacked_data=[x[1] for x in histdata],
#           sty_cycle=color_cycle+hatch_cycle,
#           labels=[x[0] for x in histdata],)

ax.legend(loc='best')
ax.set_xlabel("$T_B$ [K]")
fig4.savefig(paths.fpath('coreplots/brightest_line_histogram.png'), bbox_inches='tight')

fig6 = pl.figure(6)
fig6.clf()
ax6 = fig6.gca()
ax6.plot(cont_brightness, contincluded_line_brightness, 's')
ax6.plot([0,200], [0, 200], 'k--', alpha=0.5)
ax6.plot([0,200], [0, 200*10], 'k:', zorder=-1, alpha=0.5)
ax6.plot([0,200], [0+70, 200+70], 'k-', zorder=-1, alpha=0.2)
ax6.set_xlabel("Continuum brightness (K)")
ax6.set_ylabel("Peak line brightness (K)")
ax6.set_xlim([0, 105])
ax6.set_ylim([0, 175])
fig6.savefig(paths.fpath('coreplots/fittedpeakTB_vs_aperturecontinuum.png'), bbox_inches='tight')

fig7 = pl.figure(7)
fig7.clf()
ax7 = fig7.gca()
L, = ax7.plot(cores_merge['ApertureMass20K'][cores_merge['ApertureMass20K'] > cores_merge['T_corrected_aperturemass']],
              cores_merge['T_corrected_aperturemass'][cores_merge['ApertureMass20K'] > cores_merge['T_corrected_aperturemass']], 's')
ax7.plot(cores_merge['ApertureMass20K'][cores_merge['ApertureMass20K'] < cores_merge['T_corrected_aperturemass']],
         cores_merge['T_corrected_aperturemass'][cores_merge['ApertureMass20K'] < cores_merge['T_corrected_aperturemass']], 's',
         color=L.get_color(), alpha=0.25)
ax7.plot(cores_merge['ApertureMass20K'][cores_merge['ApertureMass20K'] < cores_merge['T_corrected_aperturemass']],
         cores_merge['ApertureMass20K'][cores_merge['ApertureMass20K'] < cores_merge['T_corrected_aperturemass']], 'o')
ylims = ax7.get_ylim()
ax7.set_xscale('log')
ax7.set_yscale('log')
ax7.plot([0,2500], [0,2500], 'k--', alpha=0.5, zorder=-5)
ax7.plot([0,2500], [0,2500/10.], 'k:', alpha=0.5, zorder=-5)
ax7.set_ylim(ylims)
ax7.set_xlim((4,2000))
ax7.set_xlabel("Mass at 20K [M$_\\odot$]")
ax7.set_ylabel("Mass at peak $T_B$ [M$_\\odot$]")
fig7.savefig(paths.fpath('coreplots/aperture_mass20K_vs_massTB.png'), bbox_inches='tight')



fig8 = pl.figure(8)
fig8.clf()

bins = np.logspace(-0.5,3.5, 30)

ax = pl.subplot(7,1,1)

H,L,P = ax.hist(cores_merge['T_corrected_peakmass'][starless], bins=bins,
                facecolor='none', histtype='stepfilled', hatch='/')
H,L,P = ax.hist(cores_merge['T_corrected_peakmass'], bins=bins,
                facecolor='none', histtype='step')
ax.set_xscale('log')
ax.set_ylim(0, H.max()+1)
ax.set_xlim(0.3, 3000)
ax.set_ylabel("Peak")

ax.xaxis.set_ticklabels([])
ax.set_xlabel("Mass at line-peak $T_B$ [$M_\\odot$]")
ax.set_xlabel("")

max_yticks=3
yloc = pl.MaxNLocator(max_yticks)
ax.yaxis.set_major_locator(yloc)


print()
print("Fit each aperture excluding free-free MASSES, THEN starless-only also ")
apertures = ('0p2', '0p4', '0p6', '0p8', '1p0', '1p5')
for ii,aperture in enumerate(apertures):
    mass = (cores_merge['T_corrected_{0}aperturemass'.format(aperture)])
    ff = (cores_merge['KUbandcont_flux{0}arcsec'.format(aperture)] /
          cores_merge['cont_flux{0}arcsec'.format(aperture)]) > 0.5


    with contextlib.redirect_stderr(devnull):
        fit = powerlaw.Fit(mass[~ff])
    print("Powerlaw mass fit for apertures {0}: {1}+/-{4}     xmin: {2}"
          "    n: {3}".format(aperture, fit.power_law.alpha,
                              fit.power_law.xmin, fit.power_law.n,
                              fit.power_law.sigma,
                             ))
    powerlaw_parameters[aperture] = {'alpha':fit.power_law.alpha,
                                     'e_alpha':fit.power_law.sigma,
                                     'xmin':fit.power_law.xmin,
                                     'number':fit.power_law.n,
                                    }

    with contextlib.redirect_stderr(devnull):
        fit = powerlaw.Fit(mass[starless & ~ff])
    print("Powerlaw mass fit for starless apertures {0}: {1}+/-{4}     xmin: {2}"
          "    n: {3}".format(aperture, fit.power_law.alpha,
                              fit.power_law.xmin, fit.power_law.n,
                              fit.power_law.sigma,
                             ))
    powerlaw_parameters['starless'+aperture] = {'alpha':fit.power_law.alpha,
                                                'e_alpha':fit.power_law.sigma,
                                                'xmin':fit.power_law.xmin,
                                                'number':fit.power_law.n, }


    ax = pl.subplot(7,1,ii+2)
    ax.set_ylabel("${0}''$".format(aperture.replace("p",".")))

    # do starless first b/c H.max() is used to set ylim
    H,L,P = ax.hist(mass[starless & ~ff], bins=bins,
                    facecolor='none',
                    #alpha=0.5,
                    histtype='stepfilled',
                    hatch='/',
                    label=aperture)


    H,L,P = ax.hist(mass[~ff], bins=bins,
                    #facecolor='none',
                    #alpha=0.5,
                    histtype='step',
                    label=aperture)


    if ii < 5:
        ax.set_xticklabels([])
        ax.get_xaxis().set_visible(False)

    yloc = pl.MaxNLocator(max_yticks)
    ax.yaxis.set_major_locator(yloc)

    ax.set_xscale('log')
    ax.set_xlim(0.3, 3000)
    ax.set_ylim(0, H.max()+1)
ax.set_xlabel("Mass at peak $T_B$ [$M_\odot$]")
pl.subplots_adjust(hspace=0)
pl.savefig(paths.fpath("coreplots/core_mass_histogram_apertureradius.png"), bbox_inches='tight')
