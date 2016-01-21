import numpy as np
from astropy import units as u
from astropy import log
from astropy import wcs
import paths
import pyregion
from astropy.io import fits
from astropy.table import Table,Column
import masscalc
import constants
from spectral_cube import SpectralCube

regions = pyregion.open(paths.rpath('outflow_ellipses.reg'))

cofile = SpectralCube.read(paths.dpath('w51_12CO_21_contsub_hires.image.pbcor.fits'))
cube = cofile.with_spectral_unit(u.km/u.s)
cube = cube.to(u.K, equivalencies=u.brightness_temperature(cube.beam.sr,
                                                           cube.wcs.wcs.restfrq*u.Hz))
pixsize = wcs.utils.proj_plane_pixel_area(cube.wcs.celestial)**0.5*u.deg
pixsize_phys = (pixsize * constants.distance).to(u.pc, u.dimensionless_angles())
dv = cube.spectral_axis.diff()[0]

results = {}
units = {'peak':u.K,
         'mean':u.K,
         'integ':u.K*u.km/u.s,
         'npix':u.dimensionless_unscaled,
         'peak_mass':u.M_sun,
         'peak_col':u.cm**-2}

for reg in regions:
    if 'text' not in reg.attr[1]:
        continue

    shreg = pyregion.ShapeList([reg])
    name,v1,v2 = reg.attr[1]['text'].split()
    log.info(name)

    scube = cube.spectral_slab(v1*u.km/u.s,
                               v2*u.km/u.s)
    scube = scube.subcube_from_ds9region(shreg)

    results[name] = {'peak': scube.max()*u.K,
                     'mean': scube.mean()*u.K,
                     'integ': scube.moment0(axis=0).mean(),
                     'npix': scube[0,:,:].mask.include().sum(),
                     'pixsize': pixsize,
                     'pixsize_phys': pixsize_phys,
                    }
    results[name]['peak_col'] = masscalc.co21_conversion_factor()*results[name]['peak']*dv
    results[name]['peak_mass'] = (results[name]['peak_col'] * pixsize_phys**2 *
                                  constants.mh2).to(u.M_sun)
    results[name]['mean_col'] = results[name]['integ'] * masscalc.co21_conversion_factor()
    results[name]['total_mass'] = (results[name]['mean_col'] * pixsize_phys**2 *
                                   constants.mh2).to(u.M_sun)

# invert the table to make it parseable by astropy...
# (this shouldn't be necessary....)
results_inv = {'name':{}}
columns = {'name':[]}
for k,v in results.items():
    results_inv['name'][k] = k
    columns['name'].append(k)
    for kk,vv in v.items():
        if kk in results_inv:
            results_inv[kk][k] = vv
            columns[kk].append(vv)
        else:
            results_inv[kk] = {k:vv}
            columns[kk] = [vv]

for c in columns:
    if c in units:
        columns[c] = columns[c] * units[c]

tbl = Table([Column(data=columns[k],
                    name=k)
             for k in ['name','peak','mean','integ','npix','peak_mass','peak_col',
                       'mean_col', 'total_mass',]])

tbl.sort('total_mass')
tbl.write(paths.tpath("outflow_co_photometry.ipac"), format='ascii.ipac')

