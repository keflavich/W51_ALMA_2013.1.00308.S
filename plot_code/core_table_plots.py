import numpy as np
import paths
from astropy.table import Table, Column
from astropy import units as u
from astropy import coordinates
import powerlaw
import pylab as pl

core_velo_tbl = Table.read(paths.tpath("core_velocities.ipac"), format="ascii.ipac")
core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
cores_merge = Table.read(paths.tpath('core_continuum_and_line.ipac'), format='ascii.ipac')

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

ax1.hist(core_phot_tbl['peak'], log=True, bins=np.logspace(-3,-0.5,15))
ax1.set_xscale('log')
ax1.set_ylim(0.3, 11)


fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(211)

fit = powerlaw.Fit(core_phot_tbl['peak'])
fit.plot_ccdf(color='k')
fit.power_law.plot_ccdf(color='r', linestyle='--')
ax2.set_ylabel("Fraction of sources")

ax3 = fig2.add_subplot(212)

fit = powerlaw.Fit(core_phot_tbl['peak'])
# doesn't work at all fit.plot_pdf(color='k')
ax3.hist(core_phot_tbl['peak'], bins=np.logspace(-3,-0.5,15),
         color='k', facecolor='none', histtype='step')
ax3.set_xscale('log')
fit.power_law.plot_pdf(color='r', linestyle='--')
ax3.set_ylim(0.3, 11)
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




beam_area = cores_merge['beam_area']
jy_to_k = (1*u.Jy).to(u.K, u.brightness_temperature(beam_area,
                                                    220*u.GHz)).mean()

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
ax5.plot([0,20], [0,20], 'k--')
ax5.set_ylim(ylims)
ax5.set_xlabel("Mass at 20K [M$_\\odot$]")
ax5.set_ylabel("Mass at peak $T_B$ [M$_\\odot$]")
fig4.savefig(paths.fpath('coreplots/mass20K_vs_massTB.png'))

fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.gca()

coords = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                              frame='fk5').galactic
ax2.plot(coords.l, coords.b, '.')
ax2.set_xlim(*ax2.get_xlim()[::-1])
ax2.set_ylabel('Galactic Latitude')
ax2.set_xlabel('Galactic Longitude')
ax2.set_aspect(1)
fig2.savefig(paths.fpath('coreplots/core_spatial_distribution.png'))


pl.draw()
pl.show()
