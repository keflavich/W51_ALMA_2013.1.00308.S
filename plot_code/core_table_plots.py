import numpy as np
import paths
from astropy.table import Table, Column
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



fig3 = pl.figure(3)
fig3.clf()
ax4 = fig3.gca()
ax4.plot(cores_merge['peak'], cores_merge['PeakLineBrightness'], 's')
ax4.set_xlabel("Continuum flux density (Jy/beam)")
ax4.set_ylabel("Peak line brightness (K)")

fig4 = pl.figure(4)
fig4.clf()
ax5 = fig4.gca()
ax5.plot(cores_merge['peak_mass'], cores_merge['T_corrected_mass'], 's')
ylims = ax5.get_ylim()
ax5.plot([0,20], [0,20], 'k--')
ax5.set_ylim(ylims)
ax5.set_xlabel("Mass at 20K [M$_\\odot$]")
ax5.set_ylabel("Mass at peak $T_B$ [M$_\\odot$]")

pl.draw()
pl.show()
