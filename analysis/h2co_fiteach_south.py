import numpy as np
import os
import paths
from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar
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
                     limited=[(True,True)]*9,
                     limits=[(0,5), (40, 70), (0.5, 4)]*3,
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
                     limited=[(True,True)]*9,
                     limits=[(0,5), (40, 70), (0.5, 4)]*3,
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
                     limited=[(True,True)]*9,
                     limits=[(0,5), (40, 70), (0.5, 4)]*3,
                     start_from_point=(99,114),
                     #start_from_point=(10,10),
                    )
    pcube322.write_fit(modelfn322, clobber=True)

# mask out components where the fitted line width or line centroid disagree
bad_0 = ((np.abs((pcube321.parcube[2,:,:] - pcube303.parcube[2,:,:]) / pcube303.parcube[2,:,:]) > 0.5) |
         (np.abs(pcube321.parcube[1] - pcube303.parcube[1]) > 1.0) |
         (pcube321.parcube[0] < 0) | (pcube303.parcube[0] < 0))
bad_1 = ((np.abs((pcube321.parcube[5,:,:] - pcube303.parcube[5,:,:]) / pcube303.parcube[5,:,:]) > 0.5) |
         (np.abs(pcube321.parcube[4] - pcube303.parcube[4]) > 1.0) |
         (pcube321.parcube[3] < 0) | (pcube303.parcube[3] < 0))
bad_2 = ((np.abs((pcube321.parcube[5,:,:] - pcube303.parcube[5,:,:]) / pcube303.parcube[5,:,:]) > 0.5) |
         (np.abs(pcube321.parcube[7] - pcube303.parcube[7]) > 1.0) |
         (pcube321.parcube[6] < 0) | (pcube303.parcube[6] < 0))

# peak intensity times line width squared (2 pi divides out)
r321303_0 = pcube321.parcube[0,:,:] / pcube303.parcube[0,:,:] * (pcube321.parcube[2,:,:] / pcube303.parcube[2,:,:])**2
r321303_1 = pcube321.parcube[3,:,:] / pcube303.parcube[3,:,:] * (pcube321.parcube[5,:,:] / pcube303.parcube[5,:,:])**2
r321303_2 = pcube321.parcube[6,:,:] / pcube303.parcube[6,:,:] * (pcube321.parcube[8,:,:] / pcube303.parcube[8,:,:])**2
r321303_0[bad_0] = np.nan
r321303_1[bad_1] = np.nan
r321303_2[bad_2] = np.nan

# repeat for 322
bad_0 = ((np.abs((pcube322.parcube[2,:,:] - pcube303.parcube[2,:,:]) / pcube303.parcube[2,:,:]) > 0.5) |
         (np.abs(pcube322.parcube[1] - pcube303.parcube[1]) > 1.0) |
         (pcube322.parcube[0] < 0) | (pcube303.parcube[0] < 0))
bad_1 = ((np.abs((pcube322.parcube[5,:,:] - pcube303.parcube[5,:,:]) / pcube303.parcube[5,:,:]) > 0.5) |
         (np.abs(pcube322.parcube[4] - pcube303.parcube[4]) > 1.0) |
         (pcube322.parcube[3] < 0) | (pcube303.parcube[3] < 0))
bad_2 = ((np.abs((pcube322.parcube[5,:,:] - pcube303.parcube[5,:,:]) / pcube303.parcube[5,:,:]) > 0.5) |
         (np.abs(pcube322.parcube[7] - pcube303.parcube[7]) > 1.0) |
         (pcube322.parcube[6] < 0) | (pcube303.parcube[6] < 0))

# peak intensity times line width squared (2 pi divides out)
r322303_0 = pcube322.parcube[0,:,:] / pcube303.parcube[0,:,:] * (pcube322.parcube[2,:,:] / pcube303.parcube[2,:,:])**2
r322303_1 = pcube322.parcube[3,:,:] / pcube303.parcube[3,:,:] * (pcube322.parcube[5,:,:] / pcube303.parcube[5,:,:])**2
r322303_2 = pcube322.parcube[6,:,:] / pcube303.parcube[6,:,:] * (pcube322.parcube[8,:,:] / pcube303.parcube[8,:,:])**2
r322303_0[bad_0] = np.nan
r322303_1[bad_1] = np.nan
r322303_2[bad_2] = np.nan

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


jtok303 = cube303.beam.jtok(cube303.with_spectral_unit(u.GHz).spectral_axis).mean().value
jtok321 = cube321.beam.jtok(cube321.with_spectral_unit(u.GHz).spectral_axis).mean().value
jtok322 = cube322.beam.jtok(cube322.with_spectral_unit(u.GHz).spectral_axis).mean().value

rcomps = {0: r321303_0,
          1: r321303_1,
          2: r321303_2,
         }

for component in (1,2,0):
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

    hdu = cube303.hdu
    hdu.data = denstemcol
    hdu.writeto(paths.dpath('h2co_fitted_denstemcol_southeast_component{0}.fits'.format(component)), clobber=True)
