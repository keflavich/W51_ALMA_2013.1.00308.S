from astropy.table import Table
from astropy import units as u
from spectral_cube import SpectralCube
import paths

ch3cn_linelist = Table.read('ch3cn_line_list.txt', format='ascii')

for cubefn, prefix, velo, dv in ((paths.dpath('longbaseline/W51northcax.SPW2_ALL_cutout_medsub_K.fits'), 'north', 60.7*u.km/u.s, 10*u.km/u.s), ):
    cube = SpectralCube.read(cubefn)
    cutouts = {}
    for line in ch3cn_linelist:
        vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
                                        rest_value=line['freq']*u.GHz)
        cutout = vcube.spectral_slab(velo-dv, velo+dv)
        name = str(line['name']).replace("=","_").replace("(","_").replace(")","_")
        if cutout.shape[0] > 1:
            cutouts[name] = cutout
            cutout.write(paths.dpath('longbaseline/velo_cutouts/{0}_{1}.fits'.format(prefix, name)),
                         overwrite=True)
            m1 = cutout.moment1(axis=0)
            m1.write(paths.dpath('longbaseline/velo_cutouts/moments/{0}_{1}_m1.fits'.format(prefix,
                                                                                            name)),
                     overwrite=True)
