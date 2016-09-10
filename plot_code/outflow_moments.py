"""
(late added) script to make outflow moment maps for various regions

these files:
#'/Users/adam/work/w51/alma/FITS/moments/w51_12co2-1_blue0to45_masked.fits'
#'/Users/adam/work/w51/alma/FITS/moments/w51_12co2-1_red73to130_masked.fits'

were made by hand
"""

import paths
from spectral_cube import SpectralCube
from astropy import units as u

# abort trap 6
from astropy import log
log.setLevel('DEBUG')

vrange_dict = {'LacyJet':
               {'SO65-54': {'blue': [30,50],
                            'red': [70,95],},
                '12CO2-1': {'blue': [30,50],
                            'red': [70,95],},
               },
              }
cutout_dict = {'LacyJet': 'north'}

for objname in vrange_dict:
    for species in vrange_dict[objname]:
        for shift in vrange_dict[objname][species]:

            fn = paths.dpath('12m/cutouts/W51_b6_12M.{0}.image.pbcor_{1}cutout.fits'
                             .format(species, cutout_dict[objname]))
            cube = SpectralCube.read(fn).with_spectral_unit(u.km/u.s,
                                                            velocity_convention='radio')
            cube.beam_threshold = 1.0

            med = cube.percentile(25, axis=0)
            medsub = cube-med

            vrange = vrange_dict[objname][species][shift]
            slab = medsub.spectral_slab(*(vrange*u.km/u.s))

            m0 = slab.moment0()

            m0.write(paths.dpath('12m/moments/{objname}_{species}_{shift}'
                                 '{v1}_{v2}.fits'
                                 .format(objname=objname,
                                         species=species,
                                         shift=shift,
                                         v1=vrange[0],
                                         v2=vrange[1],)
                                 .replace("-","m")),
                     overwrite=True,
                    )
