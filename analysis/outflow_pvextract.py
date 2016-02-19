import numpy as np
import paths
import aplpy
import pylab as pl
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import coordinates
from astropy import wcs
from astropy.io import fits

import pvextractor
from pvextractor import extract_pv_slice
from pvextractor.geometry import Path

e8mm = coordinates.SkyCoord(290.93289, 14.507833, frame='fk5', unit=(u.deg,
                                                                     u.deg))
e8blue_endpoint = coordinates.SkyCoord(290.93538, 14.506403, frame='fk5',
                                       unit=(u.deg, u.deg))
e8red_endpoint = coordinates.SkyCoord(290.92946, 14.509999, frame='fk5',
                                      unit=(u.deg, u.deg))
endpoints = coordinates.SkyCoord([e8blue_endpoint, e8red_endpoint])
e8flowxy = Path(endpoints, width=1.5*u.arcsec)

#lacyreg = pyregion.open(paths.rpath('lacyjet_segment_trace.reg'))
lacypath = pvextractor.pvregions.paths_from_regfile(paths.rpath('lacyjet_segment_trace.reg'))[0]
lacypath.width = 1.5*u.arcsec
lacyorigin = coordinates.SkyCoord(290.91564, 14.518122, frame='fk5',
                                  unit=(u.deg, u.deg))

h2velomap = fits.open('/Users/adam/work/w51/sinfoni/h2_velocity_map.fits')
h2velowidthmap = fits.open('/Users/adam/work/w51/sinfoni/h2_velocity_width_map.fits')
h2wcs = wcs.WCS(h2velomap[0].header)

def offset_to_point(xx, yy, path):
    """
    Determine the offset along the path to the nearest point on a path to the
    specified point
    """
    import shapely.geometry as geom

    if hasattr(path,'_coords'):
        line = geom.LineString(zip(path._coords.ra.deg, path._coords.dec.deg))
    else:
        raise ValueError("Path doesn't have coords - maybe need _xy?")
    point = geom.Point(xx, yy)
    return line.project(point)

parameters = {'e8':
              {'path': e8flowxy,
               'cx': 0.0035,
               'cv': 65e3,
               'wx': 0.007,
               'wv': 140e3,
               'origin': e8mm,
              },
              'lacy':
              {'path': lacypath,
               'cx': 0.0025,
               'cv': 65e3,
               'wx': 0.005,
               'wv': 300e3,
               'origin': lacyorigin,
              }
             }

for ii, (fn, stretch, vmin, vmax, source) in enumerate(
    (
     ('H77a_BDarray_speccube_briggs0_contsub_cvel_big.fits',
      'arcsinh', -0.0005, 0.003, 'lacy'),
     ('w51_SO_65-54_contsub.fits', 'arcsinh', -0.005, 0.15, 'lacy'),
     ('w51_12CO_21_contsub_hires.image.pbcor.fits', 'arcsinh', -0.02, 0.2,
      'lacy'),
     ('w51_SO_65-54_contsub.fits', 'arcsinh', -0.01, 0.3, 'e8'),
     ('w51_12CO_21_contsub_hires.image.pbcor.fits', 'arcsinh', -0.01, 0.2,
      'e8'),
    )):

    pars = parameters[source]

    cube = SpectralCube.read(paths.dpath(fn))
    cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    cube = cube.spectral_slab(-200*u.km/u.s, 200*u.km/u.s)

    outname = '{0}_{1}.{{extension}}'.format(fn[:-5], source)

    pv = extract_pv_slice(cube.hdu, pars['path'])
    #pv.data -= np.nanmin(pv.data) - 1e-3
    pv.writeto(paths.dpath('pv/'+outname.format(extension='fits')),
               clobber=True)

    origin = offset_to_point(pars['origin'].ra.deg,
                             pars['origin'].dec.deg,
                             pars['path'])

    fig = pl.figure(ii)
    fig.clf()
    FF = aplpy.FITSFigure(pv, figure=fig)
    FF.show_grayscale(aspect='auto', invert=True, stretch=stretch,
                      vmin=vmin, vmax=vmax)
    FF.show_lines([np.array([[origin, origin],
                             [pars['cv']-pars['wv']/2.,
                              pars['cv']+pars['wv']/2.,]])
                  ],
                  color='r')
    FF.recenter(pars['cx'], pars['cv'], width=pars['wx'], height=pars['wv'])

    # show() is unfortunately required before the text labels are set
    pl.draw()
    pl.show()
    FF._ax1.set_yticklabels([str(float(L.get_text())/1e3) for L in
                             FF._ax1.get_yticklabels()])
    FF._ax1.set_ylabel("Velocity (km/s)")
    FF._ax1.set_xticklabels([str(float(L.get_text())*3600) for L in
                             FF._ax1.get_xticklabels()])
    FF._ax1.set_xlabel("Offset (arcsec)")
    FF.save(paths.fpath('outflows_pv/'+outname.format(extension='png')))

    h2pixels = pars['path'].sample_points(1.0, h2wcs)
    h2velos = [h2velomap[0].data[x,y] for y,x in zip(*h2pixels)
               if x>0 and y>0 and x<h2velomap[0].data.shape[1] and y<h2velomap[0].data.shape[0]]
    h2radec = [h2wcs.wcs_pix2world(x,y,0) for y,x in zip(*h2pixels)
               if x>0 and y>0 and x<h2velomap[0].data.shape[1] and y<h2velomap[0].data.shape[0]]
    h2offset = [offset_to_point(ra,dec, pars['path'])
                for ra,dec in h2radec]
    #h2widths = h2velowidthmap[h2pixels]
    if h2offset:
        FF.show_markers(h2offset, (np.array(h2velos)*1000).tolist(), marker='x')

    outname = '{0}_{1}_with_h2.{{extension}}'.format(fn[:-5], source)
    FF.save(paths.fpath('outflows_pv/'+outname.format(extension='png')))
