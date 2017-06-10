import paths
import radio_beam
import regions
from astropy import coordinates
from astropy.io import fits
from astropy import units as u
from astropy import wcs

import warnings
warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)

objects = {'north': {'center': coordinates.SkyCoord('19:23:40.054', '14:31:05.498', unit=(u.hour, u.deg), frame='fk5'),
                     'radius': 0.15*u.arcsec,
                     'filename': 'W51n_cont_uniform.image.tt0.pbcor.fits'},
           'd2': {'center': coordinates.SkyCoord('19:23:39.820', '14:31:04.866', unit=(u.hour, u.deg), frame='fk5'),
                  'radius': 0.10*u.arcsec,
                  'filename': 'W51n_cont_uniform.image.tt0.pbcor.fits'},
           'e2e': {'center': coordinates.SkyCoord('19:23:43.969', '14:30:34.525', unit=(u.hour, u.deg), frame='fk5'),
                   'radius': 0.25*u.arcsec,
                   'filename': 'W51e2_cont_uniform.image.tt0.pbcor.fits'},
           'e8': {'center': coordinates.SkyCoord('19:23:43.906', '14:30:28.269', unit=(u.hour, u.deg), frame='fk5'),
                  'radius': 0.16*u.arcsec,
                  'filename': 'W51e2_cont_uniform.image.tt0.pbcor.fits'},
          }
for key in objects:
    obj = objects[key]
    objects[key]['aperture'] = regions.CircleSkyRegion(center=obj['center'],
                                                       radius=obj['radius'])

for objname in objects:
    obj = objects[objname]
    fh = fits.open(paths.lbpath(obj['filename']))[0]
    datawcs = wcs.WCS(fh.header)
    data = fh.data
    beam = radio_beam.Beam.from_fits_header(fh.header)
    pixreg = obj['aperture'].to_pixel(datawcs)
    mask = pixreg.to_mask()
    cutout = mask.cutout(data) * mask.data

    pkflux = cutout.max() * u.Unit(fh.header['BUNIT'])
    #pkbright = pkflux.to(u.K, beam.jtok_equiv(fh.header['CRVAL3']*u.Unit(fh.header['CUNIT3'])))
    pkbright = pkflux.to(u.K, beam.jtok_equiv(225*u.GHz))

    print("%{2:10s} Peak flux: {0:8.2f}  Peak brightness: {1:8.1f} Beam major: {3:0.3f}"
          .format(pkflux.to(u.mJy/u.beam), pkbright, objname, beam.major.to(u.arcsec)))
