import numpy as np
import os
import paths
from spectral_cube import SpectralCube
import pyspeckit
from astropy import units as u

p303 = paths.dpath('merge/W51_b6_7M_12M_natural.H2CO303_202.regrid_medsub.southcluster_cutout.fits')
p321 = paths.dpath('merge/W51_b6_7M_12M_natural.H2CO321_220.regrid_medsub.southcluster_cutout.fits')
p322 = paths.dpath('merge/W51_b6_7M_12M_natural.H2CO322_221.regrid_medsub.southcluster_cutout.fits')

cube303 = SpectralCube.read(p303)
pcube303 = pyspeckit.Cube(cube=cube303)#[:,90:130,80:120])
cube321 = SpectralCube.read(p321)
pcube321 = pyspeckit.Cube(cube=cube321)#[:,90:130,80:120])
cube322 = SpectralCube.read(p322)
pcube322 = pyspeckit.Cube(cube=cube322)#[:,90:130,80:120])

std = cube303[-10:].std(axis=0)
mask = cube303.max(axis=0) > 3*std

modelfn303 = 'model_fits/h2co_303_south_modelfits.fits'
if os.path.exists(modelfn303):
    pcube303.load_model_fit(modelfn303, npars=3, npeaks=3)
else:
    pcube303.fiteach(guesses=[0.1772, 49.474, 1.4699, 0.0451, 55.494, 2.1582, 0.0653, 62.482, 1.9253,],
                     mask=mask,
                     errmap=std.value,
                     start_from_point=(99,114),
                     #start_from_point=(10,10),
                    )
    pcube303.write_fit(modelfn303, clobber=True)

modelfn321 = 'model_fits/h2co_321_south_modelfits.fits'
if os.path.exists(modelfn321):
    pcube321.load_model_fit(modelfn321, npars=3, npeaks=3)
else:
    pcube321.fiteach(guesses=pcube303.parcube,
                     mask=mask,
                     errmap=std.value,
                     start_from_point=(99,114),
                     #start_from_point=(10,10),
                    )
    pcube321.write_fit(modelfn321, clobber=True)

modelfn322 = 'model_fits/h2co_322_south_modelfits.fits'
if os.path.exists(modelfn322):
    pcube322.load_model_fit(modelfn322, npars=3, npeaks=3)
else:
    pcube322.fiteach(guesses=pcube303.parcube,
                     mask=mask,
                     errmap=std.value,
                     start_from_point=(99,114),
                     #start_from_point=(10,10),
                    )
    pcube322.write_fit(modelfn322, clobber=True)

bad_0 = np.abs((pcube321.parcube[2,:,:] - pcube303.parcube[2,:,:]) / pcube303.parcube[2,:,:]) > 0.5
bad_1 = np.abs((pcube321.parcube[5,:,:] - pcube303.parcube[5,:,:]) / pcube303.parcube[5,:,:]) > 0.5
bad_2 = np.abs((pcube321.parcube[8,:,:] - pcube303.parcube[8,:,:]) / pcube303.parcube[8,:,:]) > 0.5

r321303_0 = pcube321.parcube[0,:,:] / pcube303.parcube[0,:,:] * (pcube321.parcube[2,:,:] / pcube303.parcube[2,:,:])**2
r321303_1 = pcube321.parcube[3,:,:] / pcube303.parcube[3,:,:] * (pcube321.parcube[5,:,:] / pcube303.parcube[5,:,:])**2
r321303_2 = pcube321.parcube[6,:,:] / pcube303.parcube[6,:,:] * (pcube321.parcube[8,:,:] / pcube303.parcube[8,:,:])**2
r321303_0[bad_0] = np.nan
r321303_1[bad_1] = np.nan
r321303_2[bad_2] = np.nan
