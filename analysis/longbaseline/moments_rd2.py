"""
Make moment maps of cutouts

Meant to operate on full-window cubes
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

files = {'e2':[
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW1_ALL_medsub_cutout.fits",
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW2_ALL_medsub_cutout.fits",
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW3_ALL_medsub_cutout.fits",
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW4_ALL_medsub_cutout.fits",
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW5_ALL_medsub_cutout.fits",
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW6_ALL_medsub_cutout.fits",
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW7_ALL_medsub_cutout.fits",
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW8_ALL_medsub_cutout.fits",
"/Volumes/passport/alma/w51/longbaseline/W51e2cax.SPW9_ALL_medsub_cutout.fits",
],}
files['e8'] = [x.replace('e2cax','e8cax') for x in files['e2']]
files['north'] = [x.replace('e2cax','northcax') for x in files['e2']]
files['northwest'] = [x.replace('e2cax','northwestcax') for x in files['e2']]

vrange_ = {'e2': (56-8, 56+8),
           'north': (61-8, 61+8),
           'northwest': (61-8, 61+8),
           'e8': (56-8, 56+8),
          }

for source,stretch in zip(('northwest','north','e2','e8',),
                          ((57,63), (58,62), (53,58), (54,58))):
    cont = fits.open(paths.dpath('longbaseline/W51{0}cax.cont.image.pbcor.fits').format('n' if 'north' in source else 'e2'))
    contwcs = wcs.WCS(cont[0].header)
    lowerleft, upperright = corners[source]['lowerleft'],corners[source]['upperright'],
    bl_x, bl_y = contwcs.celestial.wcs_world2pix(lowerleft.ra, lowerleft.dec, 0)
    tr_x, tr_y = contwcs.celestial.wcs_world2pix(upperright.ra, upperright.dec, 0)
    assert tr_y > bl_y+2
    assert tr_x > bl_x+2
    cont_cut = fits.PrimaryHDU(data=cont[0].data.squeeze()[int(bl_y):int(tr_y), int(bl_x):int(tr_x)],
                               header=contwcs.celestial[int(bl_y):int(tr_y), int(bl_x):int(tr_x)].to_header())
    assert cont_cut.data.size > 0

    vrange = vrange_[source]

    for fn in files[source]:
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

                bluvcube = vcube.spectral_slab((vrange[0]-20)*u.km/u.s, vrange[0]*u.km/u.s)
                redvcube = vcube.spectral_slab((vrange[1])*u.km/u.s, (vrange[1]+20)*u.km/u.s)
                blumax = bluvcube.max(axis=0)
                redmax = redvcube.max(axis=0)

                m1emi = slab.with_mask((slab>slabstd)).with_mask(emi).moment1()
                m1abs = slab.with_mask((slab<-slabstd)).with_mask(absorb).moment1()
                m1emi[absorb] = m1abs[absorb]
                m2emi = slab.with_mask((slab>slabstd)).with_mask(emi).moment2()
                m2abs = slab.with_mask((slab<-slabstd)).with_mask(absorb).moment2()
                m2emi[absorb] = m2abs[absorb]

                for ii in pl.get_fignums():
                    pl.close(ii)

                m0.quicklook()
                m0.FITSFigure.show_contour(cont_cut, levels=[0.0015, 0.006, 0.012],
                                           colors=['r']*3)
                m0.FITSFigure.save(paths.fpath('longbaseline/moments/{1}_{0}_mom0.png'.format(linename, source)))
                m0.write(paths.dpath('longbaseline/moments/{1}_{0}_mom0.fits'.format(linename, source)), overwrite=True)

                mx.quicklook()
                mx.FITSFigure.show_contour(cont_cut, levels=[0.0015, 0.006, 0.012],
                                           colors=['c']*3)
                mx.FITSFigure.save(paths.fpath('longbaseline/moments/{1}_{0}_max.png'.format(linename, source)))
                mx.write(paths.dpath('longbaseline/moments/{1}_{0}_max.fits'.format(linename, source)), overwrite=True)

                mx.FITSFigure.show_contour(blumax.hdu, levels=[0.005, 0.0075, 0.010],
                                           colors=['b']*3)
                mx.FITSFigure.show_contour(redmax.hdu, levels=[0.005, 0.0075, 0.010],
                                           colors=['r']*3)
                mx.FITSFigure.save(paths.fpath('longbaseline/moments/{1}_{0}_max_outflow.png'.format(linename, source)))
                mx.write(paths.dpath('longbaseline/moments/{1}_{0}_max_outflow.fits'.format(linename, source)), overwrite=True)
                
                # mask out low significance pixels
                m1emi[(mx < slabstd*2) & (mn > -slabstd*2)] = np.nan

                m1emi.quicklook()
                m1emi.FITSFigure.show_contour(cont_cut, levels=[0.0015, 0.006, 0.012],
                                              colors=['k']*3, alpha=0.5)
                m1emi.FITSFigure.show_colorscale(vmin=vrange[0], vmax=vrange[1])
                m1emi.FITSFigure.save(paths.fpath('longbaseline/moments/{1}_{0}_mom1.png'.format(linename, source)))
                m1emi.FITSFigure.show_colorscale(vmin=stretch[0], vmax=stretch[1])
                m1emi.FITSFigure.save(paths.fpath('longbaseline/moments/{1}_{0}_mom1_rescale.png'.format(linename, source)))

                m2emi[(mx < slabstd*2) & (mn > -slabstd*2)] = np.nan

                m2emi.quicklook()
                m2emi.FITSFigure.show_contour(cont_cut, levels=[0.0015, 0.006, 0.012],
                                              colors=['k']*3, alpha=0.5)
                m2emi.FITSFigure.show_colorscale(vmin=0, vmax=10)
                m2emi.FITSFigure.save(paths.fpath('longbaseline/moments/{1}_{0}_mom2.png'.format(linename, source)))
