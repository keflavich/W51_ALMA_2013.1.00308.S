"""
Make moment maps of cutouts
"""
import re
import numpy as np
from astropy import units as u
from astropy import coordinates
from spectral_cube import SpectralCube
from astropy.io import fits
from astropy import wcs
import line_to_image_list
import pyregion
import paths
import pylab as pl

regions = pyregion.open(paths.rpath("e2e8northcutouts.reg"))

corners = {reg.attr[1]['text']: {'lowerleft': coordinates.SkyCoord([reg.coord_list[:2]], frame='fk5', unit=(u.deg, u.deg)),
                                 'upperright': coordinates.SkyCoord([reg.coord_list[2:4]], frame='fk5', unit=(u.deg, u.deg)),}
           for reg in regions
          }

files = [
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW1_ALL_medsub_cutout.fits",
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW2_ALL_medsub_cutout.fits",
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW3_ALL_medsub_cutout.fits",
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW4_ALL_medsub_cutout.fits",
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW5_ALL_medsub_cutout.fits",
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW6_ALL_medsub_cutout.fits",
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW7_ALL_medsub_cutout.fits",
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW8_ALL_medsub_cutout.fits",
"/Volumes/INTENSO/w51-longbaseline/W51e2cax.SPW9_ALL_medsub_cutout.fits",
]

vrange = 56-8, 56+8

cont = fits.open(paths.dpath('longbaseline/W51e2cax.cont.image.pbcor.fits'))
contwcs = wcs.WCS(cont[0].header)
source = 'e2'
lowerleft, upperright = corners[source]['lowerleft'],corners[source]['upperright'],
bl_x, bl_y = contwcs.celestial.wcs_world2pix(lowerleft.ra, lowerleft.dec, 0)
tr_x, tr_y = contwcs.celestial.wcs_world2pix(upperright.ra, upperright.dec, 0)
assert tr_y > bl_y+2
assert tr_x > bl_x+2
cont_cut = fits.PrimaryHDU(data=cont[0].data.squeeze()[bl_y:tr_y, bl_x:tr_x],
                           header=contwcs.celestial[bl_y:tr_y, bl_x:tr_x].to_header())

for fn in files:
    cube = SpectralCube.read(fn).with_spectral_unit(u.GHz)

    for linename, freq, _, _ in line_to_image_list.line_to_image_list:
        frq = float(freq.strip("GHz"))*u.GHz

        vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
                                        rest_value=frq)

        if ((vcube.spectral_axis.min() < vrange[0]*u.km/u.s) and
            (vcube.spectral_axis.max() > vrange[1]*u.km/u.s)):

            slab = vcube.spectral_slab(vrange[0]*u.km/u.s, vrange[1]*u.km/u.s)

            print(linename, frq, slab)

            emi = cont_cut.data < 0.0015
            absorb = ~emi

            m0 = slab.moment0()
            m0std = m0.std()
            mx = slab.max(axis=0)
            mn = slab.min(axis=0)
            slabstd = slab.std()

            m1emi = slab.with_mask((slab>slabstd)).with_mask(emi).moment1()
            m1abs = slab.with_mask((slab<-slabstd)).with_mask(absorb).moment1()
            m1emi[absorb] = m1abs[absorb]
            #m2 = slab.moment2()

            for ii in pl.get_fignums():
                pl.close(ii)

            m0.quicklook()
            m0.FITSFigure.show_contour(cont_cut, levels=[0.0015, 0.006, 0.012],
                                       colors=['r']*3)
            m0.FITSFigure.save(paths.fpath('longbaseline/moments/e2e_{0}_mom0.png'.format(linename)))

            mx.quicklook()
            mx.FITSFigure.show_contour(cont_cut, levels=[0.0015, 0.006, 0.012],
                                       colors=['r']*3)
            mx.FITSFigure.save(paths.fpath('longbaseline/moments/e2e_{0}_max.png'.format(linename)))
            
            # mask out low significance pixels
            m1emi[(mx < slabstd*3) & (mn > -slabstd*3)] = np.nan

            m1emi.quicklook()
            m1emi.FITSFigure.show_contour(cont_cut, levels=[0.0015, 0.006, 0.012],
                                          colors=['k']*3, alpha=0.5)
            m1emi.FITSFigure.show_colorscale(vmin=vrange[0], vmax=vrange[1])
            m1emi.FITSFigure.save(paths.fpath('longbaseline/moments/e2e_{0}_mom1.png'.format(linename)))
