import numpy as np
from astropy import units as u
import os
import pylab as pl
import aplpy
import paths
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs
import paths
from outflow_meta import e2e, e2e_reference_vector, e8fil, north
from spectral_cube import SpectralCube
from astropy.visualization import SqrtStretch,AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

co21fn = paths.dpath('12m/w51_12CO_21_contsub_hires.image.pbcor.fits')
if not os.path.exists(co21fn):
    co21fn = paths.dpath('12m/w51_12CO_21_contsub_hires.image.pbcor.fits.gz')
so65fn = paths.dpath('12m/w51_SO_65-54_contsub.fits')


# I manually did this, which was pointless, but 7 km/s slices look best
slabs = [#(-55,-48),
         (-48,-41),
         (-41,-34),
         (-34,-27),
         (-27,-24),
         (-20,-13),
         (-13,-6),
         (-6,1),
         (1,8),
         (8,15),
         (15,22),
         (22,29),
         (29,36),
         (36,43),
         (43,50),
         (50,57),
         (57,64),
         (64,71),
         (71,78),
         (78,85),
         (85,92),
         (92,99),
         (99,106),
         (106,113),
         (113,120),
         (120,127),
         (127,134),
         (134,141),
         (141,148),
         (148,155),
         (155,162),
        ]*u.km/u.s


for source, sourcename, radius, refvec in zip((e2e, e8fil, north),
                                              ('e2e', 'e8', 'north'),
                                              (10*u.arcsec, 20*u.arcsec, 30*u.arcsec),
                                              (e2e_reference_vector, None, None)):
    for species,cubefn,mn in zip(('CO2-1','SO6-5'), (co21fn, so65fn), (-100, -2)):
        cube = SpectralCube.read(cubefn)
        cutout = Cutout2D(cube[0,:,:], source, radius, wcs=cube.wcs.celestial)

        scube = cube_cutout = cube[(slice(None),)+cutout.slices_original]
        ghzaxis = scube.with_spectral_unit(u.GHz).spectral_axis
        scube = scube.to(u.K,
                         scube.beam.jtok_equiv(ghzaxis[:,None,None]))

        if refvec is None:
            refvec_pix = [], []
        else:
            refvec_pix = scube.wcs.celestial.wcs_world2pix(refvec.spherical.lon.deg,
                                                           refvec.spherical.lat.deg,
                                                           0,
                                                          )

        # Begin channel map code here
        Nrows = 5
        Ncols = 6
        figsize = ((Ncols/Nrows)*12,12)
        fig3 = pl.figure(3)
        if any(fig3.get_size_inches() != figsize):
            pl.close(3)
        pl.figure(3, figsize=figsize).clf()
        fig, axes = pl.subplots(Nrows, Ncols,
                                sharex=True,
                                sharey=True, num=3)
        
        layers = [scube.spectral_slab(v0,v1).moment0()
                  if len(scube.spectral_slab(v0,v1)) > 1
                  else np.zeros(scube.shape[1:])*u.K
                  for v0,v1 in slabs
                 ]
        # Determine the maximum value to display
        mx = np.max([np.nanmax(x).value for x in layers])
        mn = np.max([mn, np.min([np.nanmin(x).value for x in layers])])


        for ii,(v1,v2) in enumerate(slabs):
            #pl.subplot(4, 4, ii+1)
            layer = layers[ii]
            ax = axes[ii / Ncols, ii % Ncols]
            im = ax.imshow(layer.value, norm=ImageNormalize(vmin=mn, vmax=mx,
                                                            stretch=AsinhStretch(),),
                           cmap=pl.cm.gray_r)
            axlims = ax.axis()
            ax.plot(refvec_pix[0], refvec_pix[1], 'r--', alpha=0.5)
            ax.annotate("${0:d} < v < {1:d}$".format(int(v1.value),
                                                     int(v2.value)), (0.1,
                                                                      0.8),
                        xycoords='axes fraction', color='k', fontsize='large')
            ax.axis(axlims)

        #fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.92, 0.05, 0.04, 0.9])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label("K km s$^{-1}$")

        pl.subplots_adjust(hspace=0,
                           wspace=0)

        for i in range(Nrows):
            for j in range(Ncols):
                if i == 0:
                    axes[i,j].xaxis.set_ticks_position('top')
                    pl.setp(axes[i,j].get_xticklabels(), visible=False)
                    axes[i,j].xaxis.set_ticklabels([])
                elif i == Nrows-1:
                    axes[i,j].xaxis.set_ticks_position('bottom')
                    pl.setp(axes[i,j].get_xticklabels(), visible=True)
                else:
                    axes[i,j].xaxis.set_ticks_position('none')
                    pl.setp(axes[i,j].get_xticklabels(), visible=False)
                    axes[i,j].xaxis.set_ticklabels([])

                if j == 0:
                    axes[i,j].yaxis.set_ticks_position('left')
                elif j == Ncols-1:
                    axes[i,j].yaxis.set_ticks_position('right')
                    pl.setp(axes[i,j].get_yticklabels(), visible=False)
                    axes[i,j].yaxis.set_ticklabels([])
                else:
                    axes[i,j].yaxis.set_ticks_position('none')
                    pl.setp(axes[i,j].get_yticklabels(), visible=False)
                    axes[i,j].yaxis.set_ticklabels([])

        pl.subplots_adjust(hspace=0,
                           wspace=0)
        pl.savefig(paths.fpath('outflows/{0}_{1}_channelmaps.png'.format(sourcename,species)),
                   bbox_inches='tight')
        pl.savefig(paths.fpath('outflows/{0}_{1}_channelmaps.pdf'.format(sourcename,species)),
                   bbox_inches='tight')
