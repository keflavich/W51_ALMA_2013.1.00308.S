"""
Approximate execution time 1-3 minutes
"""
import os
import numpy as np
from astropy import units as u
from astropy.utils.console import ProgressBar
from astropy import log
from astropy import coordinates
from astropy.nddata.utils import Cutout2D
import photutils
import paths
import pyregion
from astropy.io import fits
from astropy import wcs
from astropy.table import Table,Column
import masscalc
import radio_beam
import files

regions = pyregion.open(paths.rpath('cores.reg'))

contfile = fits.open(files.continuum_file)
data = contfile[0].data
beam = radio_beam.Beam.from_fits_header(files.continuum_file)

radiofilename = os.path.join(paths.vlapath,
                             'data/W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits')
radio_image = fits.open(radiofilename)

radii = (0.2,0.4,0.6,0.8,1.0,1.5)*u.arcsec

results = {}
units = {'peak':u.Jy/u.beam,
         'sum':u.Jy/u.beam,
         'npix':u.dimensionless_unscaled,
         'beam_area':u.sr,
         'peak_mass':u.M_sun,
         'peak_col':u.cm**-2,
         'RA': u.deg,
         'Dec': u.deg,
         'PeakRA': u.deg,
         'PeakDec': u.deg,
        }

for reg in ProgressBar(regions):
    if 'text' not in reg.attr[1]:
        continue
    if reg.name == 'point':
        continue

    shreg = pyregion.ShapeList([reg])
    name = reg.attr[1]['text']
    log.info(name)

    mask = shreg.get_mask(hdu=contfile[0])

    data = contfile[0].data
    results[name] = {'peak': data[mask].max(),
                     'sum': data[mask].sum(),
                     'npix': mask.sum(),
                     'beam_area': beam.sr,
                     'RA': reg.coord_list[0],
                     'Dec': reg.coord_list[1],
                    }
    results[name]['peak_mass'] = masscalc.mass_conversion_factor()*results[name]['peak']
    results[name]['peak_col'] = masscalc.col_conversion_factor(beamomega=beam.sr.value)*results[name]['peak']


    ppbeam_dict = {}
    for image, imname in ((contfile,''), (radio_image,'KUband')):
        data = image[0].data.squeeze()
        mywcs = wcs.WCS(image[0].header).celestial

        beam = radio_beam.Beam.from_fits_header(image[0].header)
        pixel_scale = np.abs(mywcs.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
        pixel_scale_as = pixel_scale.to(u.arcsec).value
        ppbeam = (beam.sr/(pixel_scale**2)).decompose().value
        if imname == '':
            ppbeam_dict['mm'] = ppbeam
        else:
            ppbeam_dict[imname] = ppbeam

        keys = ['{1}cont_flux{0}arcsec'.format(rr.value, imname) for rr in radii]
        for k in keys:
            results[name][k] = []
        log.info("Doing aperture photometry on {0}".format(imname))

        size = max(radii)*2.2
        xc,yc = reg.coord_list[:2]
        position = coordinates.SkyCoord(xc, yc, frame='fk5', unit=(u.deg,u.deg))
        cutout = Cutout2D(data, position, size, mywcs, mode='partial')

        yy,xx = np.indices(cutout.data.shape)
        xcp, ycp = cutout.wcs.wcs_world2pix(position.ra.deg, position.dec.deg, 0)
        rgrid = ((yy-ycp)**2+(xx-xcp)**2)**0.5

        for rr in radii:
            aperture = photutils.SkyCircularAperture(positions=position,
                                                     r=rr).to_pixel(cutout.wcs)
            aperture_data = photutils.aperture_photometry(data=cutout.data,
                                                          apertures=aperture,
                                                          method='exact')
            flux_jybeam = aperture_data[0]['aperture_sum']
            #flux_jybeam_average = flux_jybeam / aperture.area()
            #flux_jysr_average = flux_jybeam_average / beam.sr.value
            #flux_jysr = flux_jysr_average * aperture.area()
            #flux_jy = (flux_jysr * (pixel_scale**2).to(u.sr).value)
            colname = '{1}cont_flux{0}arcsec'.format(rr.value, imname)
            results[name][colname].append(flux_jybeam/ppbeam)

            # annoyingly, have to do this from scratch because photutils isn't
            # granular enough
            phot_mask = rgrid < (rr/pixel_scale).decompose()

            #phot_cutout = photutils.aperture.get_cutouts(cutout.data, aperture)[0][1]
            #assert phot_cutout.shape == cutouts.data.shape
            brightest_pixel = np.unravel_index(np.argmax(cutout.data*phot_mask), cutout.data.shape)
            brightest_position = cutout.wcs.wcs_pix2world(brightest_pixel[1],
                                                          brightest_pixel[0],
                                                          0)
            results[name]['PeakRA'] = brightest_position[0]
            results[name]['PeakDec'] = brightest_position[1]

# invert the table to make it parseable by astropy...
# (this shouldn't be necessary....)
results_inv = {'name':{}}
columns = {'name':[]}
for k,v in results.items():
    results_inv['name'][k] = k
    columns['name'].append(k)
    for kk,vv in v.items():
        kk = kk.replace(".","p")
        if kk in results_inv:
            results_inv[kk][k] = vv
            columns[kk].append(vv)
        else:
            results_inv[kk] = {k:vv}
            columns[kk] = [vv]

for c in columns:
    if c in units:
        columns[c] = columns[c] * units[c]
        
keys = ['{1}cont_flux{0}arcsec'.format(rr.value, imname).replace(".",'p')
        for rr in radii
        for imname in ("", "KUband")
       ]

tbl = Table([Column(data=columns[k],
                    name=k)
             for k in ['name', 'RA', 'Dec', 'PeakRA', 'PeakDec', 'peak', 'sum', 'npix', 'beam_area',
                       'peak_mass', 'peak_col']+keys])
tbl.meta = {'keywords': {'ppbeam_mm': {'value': ppbeam_dict['mm']},
                         'ppbeam_cm': {'value': ppbeam_dict['KUband']},
                         'beam_area_sr': {'value': beam.sr.value},
                         'pixel_scale_as': {'value': pixel_scale_as},
                         'mass_conversion_factor': {'value':
                                                    masscalc.mass_conversion_factor()},
                         'column_conversion_factor': {'value':
                                                      masscalc.col_conversion_factor(beam.sr)},
                        }
           }

tbl.sort('peak_mass')
tbl.write(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac',
          overwrite=True)
