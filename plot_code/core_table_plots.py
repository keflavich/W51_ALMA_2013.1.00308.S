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

pl.matplotlib.rc_file('pubfiguresrc')

core_velo_tbl = Table.read(paths.tpath("core_velocities.ipac"), format="ascii.ipac")
core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
cores_merge = Table.read(paths.tpath('core_continuum_and_line.ipac'), format='ascii.ipac')

beam_area = cores_merge['beam_area']
jy_to_k = (1*u.Jy).to(u.K, u.brightness_temperature(beam_area,
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

fit = powerlaw.Fit(core_phot_tbl['peak'])
# doesn't work at all fit.plot_pdf(color='k')
ax3.hist(core_phot_tbl['peak'], bins=np.logspace(-3,-0.5,15),
         color='k', facecolor='none', histtype='step')
ax3.set_xscale('log')
fit.power_law.plot_pdf(color='r', linestyle='--')
ax3.set_ylim(0.3, 15)
ax3.set_xlabel("Peak flux density (Jy/beam)")
ax3.set_ylabel("Number of sources")
fig2.savefig(paths.fpath('coreplots/flux_powerlaw_histogram_fit.png'))

print("Fit parameters: alpha={0}".format(fit.power_law.alpha))

fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(211)

fit = powerlaw.Fit(cores_merge['T_corrected_mass'])
fit.plot_ccdf(color='k')
fit.power_law.plot_ccdf(color='r', linestyle='--')
ax2.set_ylabel("Fraction of sources")

ax3 = fig2.add_subplot(212)

fit = powerlaw.Fit(cores_merge['T_corrected_mass'])
# doesn't work at all fit.plot_pdf(color='k')
bmin, bmax = 0.2, 6.0
bins = np.logspace(np.log10(bmin),np.log10(bmax),15)
bins = np.linspace((bmin),(bmax),15)
H,L,P = ax3.hist(cores_merge['T_corrected_mass'], bins=bins, color='k',
                 facecolor='none', histtype='step')
pdf = fit.power_law.pdf(bins)*np.max(H)
ax3.plot(bins[bins>fit.power_law.xmin], pdf, 'r--')
#fit.power_law.plot_pdf(color='r', linestyle='--')
#ax3.set_ylim(0.03, 0.5)
#ax3.set_xscale('log')
#ax3.set_yscale('log')
ax3.set_xlabel("Temperature-corrected mass")
#ax3.set_ylabel("Fraction of sources")
ax3.set_ylabel("Number of sources")
fig2.savefig(paths.fpath('coreplots/tcorr_mass_powerlaw_histogram_fit.png'))

print("Fit parameters: alpha={0}".format(fit.power_law.alpha))

fig2 = pl.figure(2)
fig2.clf()
ax3 = fig2.add_subplot(111)
bmin, bmax = 0.2, 6.0
bins = np.linspace((bmin),(bmax),15)
H,L,P = ax3.hist(cores_merge['T_corrected_mass'], bins=bins*0.99, color='k',
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
fig2.savefig(paths.fpath('coreplots/mass_histograms.png'))
peak_plot[0].set_visible(False)
H,L,P = ax3.hist(cores_merge['peak_mass'], bins=bins, color='b',
                 facecolor='none', histtype='step', label='M($20$K)',
                 linewidth=2, alpha=0.5)
ax3.set_xlim(0,7)
fig2.savefig(paths.fpath('coreplots/mass_histograms_low.png'))





fig3 = pl.figure(3)
fig3.clf()
ax4 = fig3.gca()
ax4.plot(cores_merge['peak'], cores_merge['PeakLineBrightness'], 's')
ax4.plot([0,0.4], [0, 0.4*jy_to_k.value], 'k--')
ax4.set_xlabel("Continuum flux density (Jy/beam)")
ax4.set_ylabel("Peak line brightness (K)")
ax4.set_xlim([0, 0.4])
fig3.savefig(paths.fpath('coreplots/peakTB_vs_continuum.png'))

fig4 = pl.figure(4)
fig4.clf()
ax5 = fig4.gca()
ax5.plot(cores_merge['peak_mass'], cores_merge['T_corrected_mass'], 's')
ylims = ax5.get_ylim()
ax5.plot([0,250], [0,250], 'k--')
ax5.set_ylim(ylims)
ax5.set_xlabel("Mass at 20K [M$_\\odot$]")
ax5.set_ylabel("Mass at peak $T_B$ [M$_\\odot$]")
fig4.savefig(paths.fpath('coreplots/mass20K_vs_massTB.png'))



ffile = fits.open(paths.dpath('w51_te_continuum_best.fits'))
mywcs = wcs.WCS(ffile[0].header)
img = np.isfinite(ffile[0].data)
x,y = img.max(axis=0), img.max(axis=1)
xlim = np.where(x)[0].min(),np.where(x)[0].max()
ylim = np.where(y)[0].min(),np.where(y)[0].max()

fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(1,1,1, projection=mywcs)

coords = coordinates.SkyCoord(cores_merge['RA'], cores_merge['Dec'],
                              frame='fk5')
ax2.contour(img, levels=[0.5], colors='k')
sc = ax2.scatter(coords.ra.deg, coords.dec.deg, marker='.',
                 s=120,
                 transform=ax2.get_transform('fk5'), c=cores_merge['mean_velo'],
                 edgecolors='none',
                 cmap=pl.cm.jet)
cb = pl.colorbar(mappable=sc, ax=ax2)
cb.set_label('Velocity (km s$^{-1}$)')
ax2.set_xlim(*xlim)
ax2.set_ylim(*ylim)
#ax2.set_xlim(*ax2.get_xlim()[::-1])
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
fig2.savefig(paths.fpath('coreplots/core_spatial_distribution.png'))



powerlaw_parameters = {}

fig2 = pl.figure(2)
fig2.clf()
#ax = fig2.gca()

bins = np.logspace(-3,1.05, 20)

ax = pl.subplot(7,1,1)

H,L,P = ax.hist(core_phot_tbl['peak'], bins=bins, facecolor='none',
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

apertures = ('0p2', '0p4', '0p6', '0p8', '1p0', '1p5')
for ii,aperture in enumerate(apertures):
    flux = (cores_merge['cont_flux{0}arcsec'.format(aperture)] -
            cores_merge['KUbandcont_flux{0}arcsec'.format(aperture)])
    ff = (cores_merge['KUbandcont_flux{0}arcsec'.format(aperture)] /
          cores_merge['cont_flux{0}arcsec'.format(aperture)]) > 0.5


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
pl.savefig(paths.fpath("coreplots/core_flux_histogram_apertureradius.png"))


pl.draw()
pl.show()
