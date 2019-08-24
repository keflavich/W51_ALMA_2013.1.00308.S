import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
import astropy.visualization as vis
from astropy.wcs import utils as wcsutils
import pylab as pl
import pyspeckit
import paths
from astropy import modeling
from astropy import stats


cube = SpectralCube.read('/Users/adam/work/w51/alma/FITS/longbaseline/velo_cutouts/w51e2e_csv0_j2-1_r0.5_medsub.fits')
cs21cube = subcube = cube.spectral_slab(16*u.km/u.s, 87*u.km/u.s)[::-1]

norm = vis.ImageNormalize(subcube,
                          interval=vis.ManualInterval(-0.002, 0.010),
                          stretch=vis.AsinhStretch(),
                         )

pl.rcParams['font.size'] = 12

szinch = 18
fig = pl.figure(1, figsize=(szinch,szinch))
pl.pause(0.1)
for ii in range(5):
    fig.set_size_inches(szinch, szinch)
    pl.pause(0.1)
    try:
        assert np.all(fig.get_size_inches() == np.array([szinch,szinch]))
        break
    except AssertionError:
        continue

chanmaps = subcube.plot_channel_maps(4, 4, np.arange(4, 20), norm=norm,
                                     cmap='gray_r', fig=fig,
                                     contourkwargs=dict(
                                         levels=[0.02,0.04,0.06,0.08,0.1,0.125,0.15],
                                         colors=['w']*10,
                                         linewidths=[0.5]*10,), textyloc=0.8,
                                     decimals=1,)
#chanmaps[12].coords['ra'].tick_params(rotation='vertical')
chanmaps[12].coords['ra'].set_ticklabel(exclude_overlapping=True)

for cmn in range(9, 15):
    chanmaps[cmn].title.set_color('k')

#pl.setp(chanmaps[12].xaxis.get_majorticklabels(), rotation=90)
#chanmaps[12].coords['ra'].ticklabels.set_rotation(20)



fig.savefig(paths.fpath("CS_maser_channel_maps.pdf"), bbox_inches='tight')


cs21peak = np.unravel_index(subcube.argmax(), subcube.shape)
cs21pkpos = subcube.wcs.celestial.wcs_pix2world(cs21peak[1], cs21peak[2], 0)
print(f"Center of CS 2-1: {cs21pkpos}")
peakspec = subcube[:,cs21peak[1],cs21peak[2]]


pl.rcParams['font.size'] = 16



fig2 = pl.figure(2)
fig2.set_size_inches(8,5)
fig2.clf()
ax = pl.gca()
peakspec.quicklook()
ax.set_xlabel("$V_{LSR,Radio}$ [km s$^{-1}$]")
ax.set_ylabel('Intensity [Jy beam$^{-1}$]')

ax.yaxis.tick_right()
ax2 = ax.twinx()
ax2.set_ylabel('Brightness Temperature [K]')

avgbm = subcube.average_beams(0.1)
jtok = avgbm.jtok_equiv(subcube.with_spectral_unit(u.GHz).spectral_axis.mean())
peak_K = peakspec.to(u.K, jtok).with_spectral_unit(u.km/u.s) #.quicklook()
lin, = ax2.plot(peak_K.spectral_axis.value, peak_K.value)
lin.set_visible(False)
ylim = ax2.get_ylim()
ax2.vlines(55.5, -2000, 8000, color='k', linewidth=4, zorder=-5, alpha=0.25)
ax2.set_ylim(*ylim)

pl.savefig(paths.fpath('spectra/CS2-1_maser_JyandK.pdf'), bbox_inches='tight')

sp21 = pyspeckit.Spectrum.from_hdu(peak_K.hdu)
sp21.plotter(figure=fig2)
sp21.specfit(guesses=[6000, 35, 5])
sp21.specfit(guesses=[6000, 35, 5])





cs10cube = SpectralCube.read('/Users/adam/work/w51/vla_q/FITS/cutouts/W51e2e_CS_cutout.fits').spectral_slab(16*u.km/u.s, 87*u.km/u.s)

cs10peak = np.unravel_index(cs10cube.argmax(), cs10cube.shape)
cs10pkpos = cs10cube.wcs.celestial.wcs_pix2world(cs10peak[1], cs10peak[2], 0)
print(f"Center of CS 1-0: {cs10pkpos}")
peakspec = cs10cube[:,cs10peak[1],cs10peak[2]]

fig3 = pl.figure(3)
fig3.set_size_inches(8,5)
fig3.clf()
ax = pl.gca()
peakspec.quicklook()
ax.set_xlabel("$V_{LSR,Radio}$ [km s$^{-1}$]")
ax.set_ylabel('Intensity [Jy beam$^{-1}$]')

ax.yaxis.tick_right()
ax2 = ax.twinx()
ax2.set_ylabel('Brightness Temperature [K]')

avgbm = cs10cube.average_beams(0.1)
jtok = avgbm.jtok_equiv(cs10cube.with_spectral_unit(u.GHz).spectral_axis.mean())
peak_K = peakspec.to(u.K, jtok).with_spectral_unit(u.km/u.s) #.quicklook()
lin, = ax2.plot(peak_K.spectral_axis.value, peak_K.value)
lin.set_visible(False)
ylim = ax2.get_ylim()
ax2.vlines(55.5, -2000, 8000, color='k', linewidth=4, zorder=-5, alpha=0.25)
ax2.set_ylim(*ylim)

pl.savefig(paths.fpath('spectra/CS1-0_maser_JyandK.pdf'), bbox_inches='tight')




sp10 = pyspeckit.Spectrum.from_hdu(peak_K.hdu)
sp10.plotter(figure=fig3)
sp10.error[:] = 850
sp10.specfit(guesses=[sp10.data.max(), 64, 1.5],
             limits=[(0,10000), (50, 75), (1, 20)], limited=[(True,True)]*3,
             fixed=[True, False, False],
            )

cs10pk = cs10cube.max(axis=0)
gf10 = modeling.models.Gaussian2D(amplitude=cs10pk.max(), x_mean=cs10peak[2], y_mean=cs10peak[1])
yinds,xinds = np.indices(cs10pk.shape)
fitter10 = modeling.fitting.LevMarLSQFitter()
gf10fit = fitter10(model=gf10, x=xinds, y=yinds, z=cs10pk, weights=1/stats.mad_std(cs10pk.value))

# correct cs10 wcs
ww_cs10_corr = cs10pk.wcs
ww_cs10_corr.wcs.radesys = 'ICRS'
center10 = ww_cs10_corr.wcs_pix2world(gf10fit.x_mean, gf10fit.y_mean, 0)
pixscale10 = wcsutils.proj_plane_pixel_area(cs10pk.wcs)**0.5*u.deg
centererr10 = np.array([fitter10.fit_info['param_cov'][(1,1)],fitter10.fit_info['param_cov'][(2,2)]])**0.5*pixscale10

cs21pk = cs21cube.max(axis=0)
gf21 = modeling.models.Gaussian2D(amplitude=cs21pk.max(), x_mean=cs21peak[2], y_mean=cs21peak[1])
yinds,xinds = np.indices(cs21pk.shape)
fitter21 = modeling.fitting.LevMarLSQFitter()
fitter21(model=gf21, x=xinds, y=yinds, z=cs21pk)
gf21fit = fitter21(model=gf21, x=xinds, y=yinds, z=cs21pk, weights=1/stats.mad_std(cs21pk.value))
center21 = cs21pk.wcs.wcs_pix2world(gf21fit.x_mean, gf21fit.y_mean, 0)
pixscale21 = wcsutils.proj_plane_pixel_area(cs21pk.wcs)**0.5*u.deg
centererr21 = np.array([fitter21.fit_info['param_cov'][(1,1)],fitter21.fit_info['param_cov'][(2,2)]])**0.5*pixscale21


fwhm_scale = np.sqrt(8*np.log(2))
with open('../paper_cs_maser/linefit_table.tex', 'w') as fh:
    fh.write("CS J=1-0 &")
    fh.write(f"{int(np.round(sp10.specfit.parinfo[0].value, -2)):10d} &")
    fh.write(f"${sp10.specfit.parinfo[1].value:10.1f} \pm {sp10.specfit.parinfo[1].error:10.1f}$ &")
    fh.write(f"${sp10.specfit.parinfo[2].value*fwhm_scale:10.1f} \pm {sp10.specfit.parinfo[2].error*fwhm_scale:10.1f}$ &")
    fh.write(f"${center10[0]:15.{-int(np.floor(np.log10(centererr21[0].value)))}f} \pm {centererr10[0].value:15.{-int(np.floor(np.log10(centererr21[0].value)))}f}$ &")
    fh.write(f"${center10[1]:15.{-int(np.floor(np.log10(centererr21[1].value)))}f} \pm {centererr10[1].value:15.{-int(np.floor(np.log10(centererr21[1].value)))}f}$ \\\\\n")

    fh.write("CS J=2-1 &")
    fh.write(f"${int(np.round(sp21.specfit.parinfo[0].value, -2)):10d} \pm {int(np.round(sp21.specfit.parinfo[0].error, -2)):10d}$ &")
    fh.write(f"${sp21.specfit.parinfo[1].value:10.2f} \pm {sp21.specfit.parinfo[1].error:10.2f}$ &")
    fh.write(f"${sp21.specfit.parinfo[2].value*fwhm_scale:10.2f} \pm {sp21.specfit.parinfo[2].error*fwhm_scale:10.2f}$ &")
    fh.write(f"${center21[0]:15.{-int(np.floor(np.log10(centererr21[0].value)))}f} \pm {centererr10[0].value:15.{-int(np.floor(np.log10(centererr21[0].value)))}f}$ &")
    fh.write(f"${center21[1]:15.{-int(np.floor(np.log10(centererr21[1].value)))}f} \pm {centererr10[1].value:15.{-int(np.floor(np.log10(centererr21[1].value)))}f}$ \\\\\n")

separation = (((center10[0] - center21[0])**2 + (center10[1] - center21[1])**2)**0.5 * u.deg).to(u.arcsec)
seperr = (((centererr10**2).sum() + (centererr21**2).sum())**0.5).to(u.arcsec)

print(f"Separation is {separation.value:0.3f} +/- {seperr.value:0.3f} arcsec")
print(f"Separation is {separation.value*5400:0.1f} +/- {seperr.value*5400:0.1f} AU")

