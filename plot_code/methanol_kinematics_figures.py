import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from astropy import units as u
from astropy.io import fits
import pylab as pl
import paths
import reproject

blue_fits_fn = paths.dpath('moments/w51_12co2-1_blue0to45_masked.fits')
red_fits_fn = paths.dpath('moments/w51_12co2-1_red73to130_masked.fits')
red_fits = fits.open(red_fits_fn)
blue_fits = fits.open(blue_fits_fn)

vrange = {'e2wide': (53,61),
          'e8wide': (54,64),
          'northwide': (56,64),
         }
threshold= {'e2wide': 40,
            'e8wide': 30,
            'northwide': 30,
           }
wrange = {'e2wide': (6,20),
          'e8wide': (6,15),
          'northwide': (6,14),
         }
cutout = {'e2wide': [slice(20,-20),slice(15,-25)],
          'e8wide': [slice(10,-20),slice(20,-20)],
          'northwide': [slice(58,120),slice(28,88)],
         }
figsize = {'e2wide': (8,7.75),
           'e8wide': (8,8.725),
           'northwide': (8,8),
          }
labelline = {'e2wide': ([8,18],[65,65]),
             'e8wide': ([8,18],[80,80]),
             'northwide': ([8,18],[50,50]),
            }
outflow_levels = {'e2wide': ([1],[1]),
                  'e8wide': ([0.5],[1]),
                  'northwide': ([0.5,1,2],[0.5,1,2]),
                 }

for sourcename in ('e2wide', 'e8wide', 'northwide'):
    linename = 'CH3OH808-716'
    suffix = ''

    cslice = cutout[sourcename]
    wmin,wmax = wrange[sourcename]
    vmin,vmax = vrange[sourcename]
    linex, liney = labelline[sourcename]
    levels_r,levels_b = outflow_levels[sourcename]

    m1fitsfn = paths.dpath("chemslices/chemical_m1_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
    m2fitsfn = paths.dpath("chemslices/chemical_m2_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
    maxfitsfn = paths.dpath("chemslices/chemical_max_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))

    ch3oh87_m1 = fits.open(m1fitsfn)[0]
    ch3oh87_m2 = fits.open(m2fitsfn)[0]
    ch3oh87_max = fits.open(maxfitsfn)[0]
    ch3oh87_mask = ch3oh87_max.data > threshold[sourcename]
    ch3oh87_m1.data[~ch3oh87_mask] = np.nan
    ch3oh87_m2.data[~ch3oh87_mask] = np.nan
    ch3oh87_m2_fwhm = ch3oh87_m2.data**0.5 * np.sqrt(8*np.log(2))

    linename = 'CH3OH1029-936'
    m1fitsfn = paths.dpath("chemslices/chemical_m1_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
    m2fitsfn = paths.dpath("chemslices/chemical_m2_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))
    maxfitsfn = paths.dpath("chemslices/chemical_max_slabs_{0}_{1}{2}.fits".format(sourcename, linename, suffix))

    ch3oh109_m1 = fits.open(m1fitsfn)[0]
    ch3oh109_m2 = fits.open(m2fitsfn)[0]
    ch3oh109_max = fits.open(maxfitsfn)[0]
    ch3oh109_mask = ch3oh109_max.data > threshold[sourcename]
    ch3oh109_m1.data[~ch3oh109_mask] = np.nan
    ch3oh109_m2.data[~ch3oh109_mask] = np.nan
    ch3oh109_m2_fwhm = ch3oh109_m2.data**0.5 * np.sqrt(8*np.log(2))

    co_red,_ = reproject.reproject_interp(red_fits, ch3oh87_m1.header)
    co_blue,_ = reproject.reproject_interp(blue_fits, ch3oh87_m1.header)

    pl.close(1)
    fig1 = pl.figure(1, figsize=figsize[sourcename])
    fig1.clf()

    gs1 = gridspec.GridSpec(2,2)
    gs1.update(wspace=0.0, hspace=0.0)

    ax1 = fig1.add_subplot(gs1[0])
    ax1.imshow(ch3oh87_m1.data[cslice], vmin=vmin, vmax=vmax, origin='lower', cmap='seismic', interpolation='nearest')

    ax2 = fig1.add_subplot(gs1[1])
    mappable2 = ax2.imshow(ch3oh109_m1.data[cslice], vmin=vmin, vmax=vmax, origin='lower', cmap='seismic', interpolation='nearest')

    ax3 = fig1.add_subplot(gs1[2])
    ax3.imshow(ch3oh87_m2_fwhm[cslice], vmin=wmin, vmax=wmax, origin='lower', cmap='inferno', interpolation='nearest')

    ax4 = fig1.add_subplot(gs1[3])
    mappable4 = ax4.imshow(ch3oh109_m2_fwhm[cslice], vmin=wmin, vmax=wmax, origin='lower', cmap='inferno', interpolation='nearest')

    ax_lims = ax4.axis()

    ax4.plot(linex, liney, 'k-')
    length = (10 * ch3oh87_m1.header['CDELT2'] * u.deg).to(u.arcsec)
    ax4.annotate('0.5 arcsec',(np.mean(linex),np.mean(liney)+3), fontsize=10, horizontalalignment='center')
    ax4.annotate('2700 AU',(np.mean(linex),np.mean(liney)-5), fontsize=10, horizontalalignment='center')
    ax4.axis(ax_lims)

    axlims = ax4.axis()
    ax4.contour(co_blue[cslice], levels=levels_b, colors=['b']*len(levels_b))
    ax4.contour(co_red[cslice], levels=levels_r, colors=['r']*len(levels_r))
    ax4.axis(axlims)

    for ax in (ax1,ax2,ax3,ax4):
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

    box = ax2.get_position()
    axColor = pl.axes([box.x0*1.02 + box.width, box.y0+0.02, 0.02, box.height*0.9])
    cb1 = pl.colorbar(mappable2, cax=axColor, orientation="vertical")
    cb1.set_label("$V_{LSR}$ [km s$^{-1}$]")

    box = ax4.get_position()
    axColor = pl.axes([box.x0*1.02 + box.width, box.y0+0.02, 0.02, box.height*0.9])
    cb2 = pl.colorbar(mappable4, cax=axColor, orientation="vertical")
    cb2.set_label("$\sigma_{\\mathrm{FWHM}}$ [km s$^{-1}$]")

    #pl.tight_layout()
    pl.savefig(paths.fpath("methanol_velocity_moments_{sourcename}.png".format(sourcename=sourcename)), bbox_inches='tight', dpi=200)
    pl.draw()
    pl.show()
