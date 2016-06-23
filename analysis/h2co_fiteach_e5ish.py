import numpy as np
import os
import paths
from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar
import pyspeckit
from astropy import units as u
from astropy.io import fits

p303 = paths.dpath('merge/W51_b6_7M_12M.H2CO303_202.image.pbcor.fits')
p321 = paths.dpath('merge/W51_b6_7M_12M.H2CO321_220.image.pbcor.fits')
p322 = paths.dpath('merge/W51_b6_7M_12M.H2CO322_221.image.pbcor.fits')

### POSSIBLE OFFSET BY 1-2 PIXELS: no guarantee the 3 cubes are on identical grids

cube303 = SpectralCube.read(p303)[:,1411:1620,968:1202]
cube303.beam_threshold=0.1
beam303 = cube303._average_beams(1)
med303 = cube303.with_mask(((cube303.spectral_axis < 35*u.km/u.s) |
                            (cube303.spectral_axis >
                             85*u.km/u.s))[:,None,None]).median(axis=0)
pcube303 = pyspeckit.Cube(cube=SpectralCube(cube303._data-med303.value, cube303.wcs,
                                            cube303.mask, beam=beam303,
                                            read_beam=False).with_spectral_unit(u.km/u.s))
cube321 = SpectralCube.read(p321)[:,1411:1620,968:1202]
cube321.beam_threshold=0.1
beam321 = cube321._average_beams(1)
med321 = cube321.with_mask(((cube321.spectral_axis < 35*u.km/u.s) |
                            (cube321.spectral_axis >
                             85*u.km/u.s))[:,None,None]).median(axis=0)
pcube321 = pyspeckit.Cube(cube=SpectralCube(cube321._data-med321.value, cube321.wcs,
                                            cube321.mask, beam=beam321,
                                            read_beam=False).with_spectral_unit(u.km/u.s))
cube322 = SpectralCube.read(p322)[:,1411:1620,968:1202]
cube322.beam_threshold=0.1
beam322 = cube322._average_beams(1)
med322 = cube322.with_mask(((cube322.spectral_axis < 35*u.km/u.s) |
                            (cube322.spectral_axis >
                             85*u.km/u.s))[:,None,None]).median(axis=0)
pcube322 = pyspeckit.Cube(cube=SpectralCube(cube322._data-med322.value, cube322.wcs,
                                            cube322.mask, beam=beam322,
                                            read_beam=False).with_spectral_unit(u.km/u.s))

std = cube303[-10:].std(axis=0)
mask = (cube303.max(axis=0) > 3*std) & (cube303.max(axis=0) > 100*u.mJy)

modelfn303 = 'model_fits/h2co_303_e5ish_modelfits.fits'
if os.path.exists(modelfn303):
    pcube303.load_model_fit(modelfn303, npars=3, npeaks=2)
else:
    pcube303.fiteach(guesses=[0.05, 54, 1.7, 0.21, 60.5, 1.3253,],
                     maskmap=mask,
                     errmap=std.value,
                     limited=[(True,True)]*6,
                     limits=[(0,5), (47, 56), (0.2, 4),
                             (0,5), (57, 63), (0.2, 4),
                            ],
                     integral=False,
                     start_from_point=(99,114),
                     #start_from_point=(10,10),
                    )
    pcube303.write_fit(modelfn303, clobber=True)

guesses = pcube303.parcube
guesses[1,:,:][(guesses[1,:,:]<47)  | (guesses[1,:,:]>56)] = 55
guesses[2,:,:][(guesses[2,:,:]<0.2) | (guesses[2,:,:]>4)] = 1.0
guesses[4,:,:][(guesses[4,:,:]<57)  | (guesses[4,:,:]>63)] = 60
guesses[5,:,:][(guesses[5,:,:]<0.2) | (guesses[5,:,:]>4)] = 1.0

modelfn321 = 'model_fits/h2co_321_e5ish_modelfits.fits'
if os.path.exists(modelfn321):
    pcube321.load_model_fit(modelfn321, npars=3, npeaks=2)
else:
    pcube321.fiteach(guesses=guesses,
                     maskmap=mask,
                     errmap=std.value,
                     limited=[(True,True)]*6,
                     integral=False,
                     limits=[(0,5), (47, 56), (0.2, 4),
                             (0,5), (57, 63), (0.2, 4),
                            ],
                     start_from_point=(99,114),
                     #start_from_point=(10,10),
                    )
    pcube321.write_fit(modelfn321, clobber=True)

modelfn322 = 'model_fits/h2co_322_e5ish_modelfits.fits'
if os.path.exists(modelfn322):
    pcube322.load_model_fit(modelfn322, npars=3, npeaks=2)
else:
    pcube322.fiteach(guesses=guesses,
                     maskmap=mask,
                     errmap=std.value,
                     limited=[(True,True)]*6,
                     integral=False,
                     limits=[(0,5), (47, 56), (0.2, 4),
                             (0,5), (57, 63), (0.2, 4),
                            ],
                     start_from_point=(99,114),
                     #start_from_point=(10,10),
                    )
    pcube322.write_fit(modelfn322, clobber=True)

# mask out components where the fitted line width or line centroid disagree
# because of blended components, more liberal here than in South
bad_0 = ((np.abs((pcube321.parcube[2,:,:] - pcube303.parcube[2,:,:]) / pcube303.parcube[2,:,:]) > 1.0) |
         (np.abs(pcube321.parcube[1] - pcube303.parcube[1]) > 1.0) |
         (pcube321.parcube[0] < 0) | (pcube303.parcube[0] < 0))
bad_1 = ((np.abs((pcube321.parcube[5,:,:] - pcube303.parcube[5,:,:]) / pcube303.parcube[5,:,:]) > 1.0) |
         (np.abs(pcube321.parcube[4] - pcube303.parcube[4]) > 1.0) |
         (pcube321.parcube[3] < 0) | (pcube303.parcube[3] < 0))

# peak intensity times line width squared (2 pi divides out)
r321303_0 = pcube321.parcube[0,:,:] / pcube303.parcube[0,:,:] * (pcube321.parcube[2,:,:] / pcube303.parcube[2,:,:])**2
r321303_1 = pcube321.parcube[3,:,:] / pcube303.parcube[3,:,:] * (pcube321.parcube[5,:,:] / pcube303.parcube[5,:,:])**2
r321303_0[bad_0] = np.nan
r321303_1[bad_1] = np.nan

# repeat for 322
bad_0 = ((np.abs((pcube322.parcube[2,:,:] - pcube303.parcube[2,:,:]) / pcube303.parcube[2,:,:]) > 1.0) |
         (np.abs(pcube322.parcube[1] - pcube303.parcube[1]) > 1.0) |
         (pcube322.parcube[0] < 0) | (pcube303.parcube[0] < 0))
bad_1 = ((np.abs((pcube322.parcube[5,:,:] - pcube303.parcube[5,:,:]) / pcube303.parcube[5,:,:]) > 1.0) |
         (np.abs(pcube322.parcube[4] - pcube303.parcube[4]) > 1.0) |
         (pcube322.parcube[3] < 0) | (pcube303.parcube[3] < 0))

# peak intensity times line width squared (2 pi divides out)
r322303_0 = pcube322.parcube[0,:,:] / pcube303.parcube[0,:,:] * (pcube322.parcube[2,:,:] / pcube303.parcube[2,:,:])**2
r322303_1 = pcube322.parcube[3,:,:] / pcube303.parcube[3,:,:] * (pcube322.parcube[5,:,:] / pcube303.parcube[5,:,:])**2
r322303_0[bad_0] = np.nan
r322303_1[bad_1] = np.nan

from h2co_modeling import lte_model
modelratio321303 = lte_model.T_321/lte_model.T_303
modelratio322303 = lte_model.T_322/lte_model.T_303

#def ratio_to_temperature(data, modelratio=modelratio321303):
#    vals = data[np.isfinite(data)]
#    tems = np.interp(vals, modelratio[np.isfinite(modelratio)], np.array(lte_model.tem)[np.isfinite(modelratio)])
#    newr = data.copy()
#    newr[np.isfinite(data)] = tems
#    return newr
#
#t321303_0 = ratio_to_temperature(r321303_0, modelratio=modelratio321303)
#t321303_1 = ratio_to_temperature(r321303_1, modelratio=modelratio321303)
#t321303_2 = ratio_to_temperature(r321303_2, modelratio=modelratio321303)
#t322303_0 = ratio_to_temperature(r322303_0, modelratio=modelratio322303)
#t322303_1 = ratio_to_temperature(r322303_1, modelratio=modelratio322303)
#t322303_2 = ratio_to_temperature(r322303_2, modelratio=modelratio322303)

from h2co.constrain_parameters import paraH2COmodel
pmod = paraH2COmodel()

ygrid,xgrid = np.indices(r321303_0.shape)


jtok303 = beam303.jtok(cube303.with_spectral_unit(u.GHz).spectral_axis).mean().value
jtok321 = beam321.jtok(cube321.with_spectral_unit(u.GHz).spectral_axis).mean().value
jtok322 = beam322.jtok(cube322.with_spectral_unit(u.GHz).spectral_axis).mean().value

rcomps = {0: r321303_0,
          1: r321303_1,
         }

for component in (1,0):
    pb = ProgressBar(xgrid.size)

    rat = rcomps[component]
    denstemcol = np.zeros([3, rat.shape[0], rat.shape[1]])
    denstemcol[:] = np.nan

    for xx,yy in zip(xgrid.flat, ygrid.flat):
        if not np.isnan(rat[yy,xx]):
            pmod.set_constraints(taline303=pcube303.parcube[component*3,yy,xx]*jtok303,
                                 etaline303=pcube303.errcube[component*3,yy,xx]*jtok303,
                                 taline321=pcube321.parcube[component*3,yy,xx]*jtok321,
                                 etaline321=pcube321.errcube[component*3,yy,xx]*jtok321,
                                 taline322=pcube322.parcube[component*3,yy,xx]*jtok322,
                                 etaline322=pcube322.errcube[component*3,yy,xx]*jtok322,
                                 linewidth=pcube303.parcube[component*3+1,yy,xx],
                                 fit_intensity=True)
            constraints = pmod.get_parconstraints()
            denstemcol[:, yy, xx] = (constraints['density_chi2'], constraints['temperature_chi2'], constraints['column_chi2'])
        pb.update()

    hdu = fits.PrimaryHDU(data=denstemcol, header=cube303.header)
    hdu.writeto(paths.dpath('h2co_fitted_denstemcol_e5ish_component{0}.fits'.format(component)), clobber=True)
