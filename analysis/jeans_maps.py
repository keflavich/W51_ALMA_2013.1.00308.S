import numpy as np
import radio_beam
import pyregion
import paths
#import image_tools
from astropy import wcs
from astropy import constants
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import units as u
from astropy import coordinates
from astropy.convolution import convolve, Gaussian2DKernel
import pylab as pl
import itertools
import masscalc
import dust_emissivity
import reproject

def jeans_maps(regions, size=u.Quantity([1.25,1.25], u.arcsec), smooth=0):
    names = [r.attr[1]['text'] for r in regions]
    center_positions = coordinates.SkyCoord([r.coord_list
                                             for r in regions],
                                            unit=(u.deg, u.deg),
                                            frame='fk5')

    fn = "W51_te_continuum_best.fits"
    fh = fits.open(paths.dpath(fn))
    mywcs = wcs.WCS(fh[0].header)
    bm = radio_beam.Beam.from_fits_header(fh[0].header)
    pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg
    pixscale_cm = (pixscale * masscalc.distance).to(u.cm, u.dimensionless_angles())
    #ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam

    mass_maps = {}

    for ii,(name,position) in enumerate(zip(names, center_positions)):
        cutout = Cutout2D(fh[0].data, position, size, wcs=mywcs)

        source = 'e2' if name == 'e2e' else name

        temperature_map_fn = paths.dpath('12m/moments/CH3OH_{0}_cutout_temperaturemap.fits'.format(source))
        temperature_map_fh = fits.open(temperature_map_fn)

        # this whole section is copied from overlay_contours_on_ch3oh
        ch3ohN_hdul = fits.open(paths.dpath('12m/moments/CH3OH_{0}_cutout_columnmap.fits'.format(source)))
        #ch3ohT_hdul = fits.open(paths.dpath('12m/moments/CH3OH_{0}_cutout_temperaturemap.fits'.format(source)))
        #bigwcs = wcs.WCS(ch3ohT_hdul[0].header)
        #bigpixscale = (bigwcs.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg
        #ch3ohN = ch3ohN_hdul[0].data
        #ch3ohT = ch3ohT_hdul[0].data
        dust_brightness,wts = reproject.reproject_interp(fits.open(paths.dpath('W51_te_continuum_best.fits')),
                                                         ch3ohN_hdul[0].header)

        temwcs = wcs.WCS(temperature_map_fh[0].header)
        temcutout = Cutout2D(temperature_map_fh[0].data, position, size, wcs=temwcs)
        #maskcutout = Cutout2D(mask.astype('float'), position, size, wcs=bigwcs)

        bm_cm = ((bm.major**2 + bm.minor**2)**0.5 * masscalc.distance).to(u.cm, u.dimensionless_angles())

        if smooth != 0:
            stddev_pix = ((smooth**2 - bm_cm**2)**0.5 / pixscale_cm).decompose()
            print('stddev_pix: {0}'.format(stddev_pix.decompose()))
            kernel = Gaussian2DKernel(stddev_pix)

            kernel.normalize('integral')
            smdata = convolve(cutout.data, kernel)

            # temperature should be the average temperature
            kernel.normalize('integral')
            temmap = convolve(temcutout.data, kernel)
            mass_map = (smdata * masscalc.mass_conversion_factor(TK=temcutout.data)).to(u.M_sun)

            bm_cm = smooth
        else:
            mass_map = (cutout.data * masscalc.mass_conversion_factor(TK=temcutout.data)).to(u.M_sun)
            temmap = temcutout.data

        # volume of a gaussian: sqrt(2 pi)^N r^N
        volume = (2*np.pi)**(1.5) * bm_cm**3
        density_map = (mass_map / volume / (2.8*u.Da)).to(u.cm**-3)
        print("Scale: {0}".format(bm_cm.to(u.au)))

        c_s_map = ((constants.k_B * temmap*u.K / (2.4*u.Da))**0.5).to(u.km/u.s)
        MJ_map = (np.pi/6. * c_s_map**3 / (constants.G**1.5 *
                                           (2.8*u.Da*density_map)**0.5)).to(u.M_sun)

        fig = pl.figure(ii)
        fig.clf()
        ax1=pl.subplot(2,3,1)
        im1 = ax1.imshow(mass_map.value, cmap='gray', vmin=0, vmax=10)
        ax1.set_title("Measured mass")
        fig.colorbar(im1)

        ax2=pl.subplot(2,3,2)
        im2 = ax2.imshow(MJ_map.value, cmap='gray', vmin=0, vmax=10)
        ax2.set_title("Jeans mass")
        fig.colorbar(im2)
        ax3=pl.subplot(2,3,3)
        ax3.set_title("Black: unstable")
        # 1 -> white
        # 0 -> black
        # mass > M_J = unstable
        # M_J > mass = stable = 1 = white
        ax3.imshow(MJ_map > mass_map, cmap='gray', vmin=0, vmax=1)

        ax4=pl.subplot(2,3,4)
        ax4.set_title("log Density")
        im4 = ax4.imshow(np.log10(density_map.value), cmap='gray', vmin=5, vmax=7.5)
        fig.colorbar(im4)

        ax5=pl.subplot(2,3,5)
        ax5.set_title("Temperature")
        im5 = ax5.imshow(temmap, cmap='gray', vmin=0, vmax=700)
        fig.colorbar(im5)

        for jj in range(1,6):
            pl.subplot(2,3,jj).xaxis.set_major_formatter(pl.NullFormatter())
            pl.subplot(2,3,jj).yaxis.set_major_formatter(pl.NullFormatter())


        if smooth == 0:
            fig.savefig(paths.fpath("jeans_maps_{0}.png".format(name)))
        else:
            fig.savefig(paths.fpath("jeans_maps_{0}_smooth{1}.png".format(name, smooth)))

                
        mass_maps[name] = (mass_map, MJ_map)

    return mass_maps

regions = pyregion.open(paths.rpath("hmcore_centroids.reg"))
jmps = jeans_maps(regions)
jmps2 = jeans_maps(regions, smooth=1700*u.au)
jmps3 = jeans_maps(regions, smooth=2000*u.au)
jmps4 = jeans_maps(regions, smooth=3000*u.au)
