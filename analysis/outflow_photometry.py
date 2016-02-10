import string
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
pixsize = wcs.utils.proj_plane_pixel_area(cube.wcs.celestial)**0.5*u.deg
pixsize_phys = (pixsize * constants.distance).to(u.pc, u.dimensionless_angles())
dv = cube.spectral_axis.diff()[0]
tb_equiv = u.brightness_temperature(cube.beam.sr, cube.wcs.wcs.restfrq*u.Hz)

results = {}
units = {'peak':u.K,
         'mean':u.K,
         'integ':u.K*u.km/u.s,
         'npix':u.dimensionless_unscaled,
         'peak_mass':u.M_sun,
         'peak_col':u.cm**-2,
         'mean_col':u.cm**-2,
         'total_mass':u.M_sun,
         'pixsize':u.deg,
         'pixsize_phys':u.pc,
         'v1': u.km/u.s,
         'v2': u.km/u.s,
         'momentum': u.M_sun*u.km/u.s,
        }

suffixes = string.ascii_letters

centralv_table = Table.read(paths.tpath('core_velocities.ipac'),
                            format='ascii.ipac')

for reg in regions:
    if 'text' not in reg.attr[1]:
        continue

    shreg = pyregion.ShapeList([reg])
    name_,v1,v2 = reg.attr[1]['text'].split()
    v1 = float(v1)*u.km/u.s
    v2 = float(v2)*u.km/u.s

    # get the name before a suffix is appended
    if name_ in centralv_table['SourceID']:
        central_velo = centralv_table['mean_velo'][centralv_table['SourceID']==name_][0] * u.km/u.s
    else:
        central_velo = None

    ii = 0
    name = name_+'_a'
    while name in results:
        ii = ii+1
        name = name_+"_"+suffixes[ii]

    log.info(name)

    scube = cube.spectral_slab(v1, v2)
    scube = scube.subcube_from_ds9region(shreg)
    scube = scube.to(u.K, equivalencies=tb_equiv)

    sumspec = scube.sum(axis=(1,2))

    results[name] = {'peak': scube.max(),
                     'mean': scube.mean(),
                     'integ': np.nanmean(scube.moment0(axis=0).value)*u.K*u.km/u.s,
                     'npix': scube.mask.include(view=(0, slice(None), slice(None))).sum(),
                     'pixsize': pixsize,
                     'pixsize_phys': pixsize_phys,
                     'v1': v1,
                     'v2': v2,
                     'SourceID': name_,
                    }
    results[name]['peak_col'] = masscalc.co21_conversion_factor(results[name]['peak'])*results[name]['peak']*dv/masscalc.co_abund
    results[name]['peak_mass'] = (results[name]['peak_col'] * pixsize_phys**2 *
                                  constants.mh2).to(u.M_sun)
    results[name]['mean_col'] = results[name]['integ'] * masscalc.co21_conversion_factor(results[name]['peak'])/masscalc.co_abund
    results[name]['total_mass'] = (results[name]['mean_col'] * pixsize_phys**2 *
                                   results[name]['npix'] *
                                   constants.mh2).to(u.M_sun)
    channel_to_mass = (masscalc.co21_conversion_factor(results[name]['peak']) /
                       masscalc.co_abund * pixsize_phys**2 *
                       results[name]['npix'] * constants.mh2)
    if central_velo is not None:
        results[name]['momentum'] = u.Quantity((np.abs(scube.spectral_axis - central_velo) * sumspec * dv * channel_to_mass).sum()).to(u.Msun * u.km/u.s)
    else:
        results[name]['momentum'] = np.nan * u.Msun * u.km/u.s

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
        if hasattr(columns[c], 'value'):
            pass
        elif hasattr(columns[c][0], 'value'):
            columns[c] = np.array([x.value for x in columns[c]]) * units[c]
        else:
            columns[c] = np.array(columns[c])

tbl = Table([Column(data=columns[k],
                    name=k)
             for k in ['name', 'SourceID', 'peak', 'mean', 'integ', 'npix',
                       'pixsize', 'pixsize_phys', 'peak_mass', 'peak_col',
                       'mean_col', 'total_mass', 'momentum']])

tbl.sort('total_mass')
tbl.write(paths.tpath("outflow_co_photometry.ipac"), format='ascii.ipac')
