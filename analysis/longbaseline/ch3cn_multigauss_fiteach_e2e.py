from ch3cn_fits import SpectralCube, pyspeckit, fits, u, np
import scipy.stats
import os
import paths
T=True
F=False

cubefn = paths.dpath('longbaseline/W51e2e_CH3CN_cutout.fits')
cube = SpectralCube.read(cubefn).minimal_subcube()
contcubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                                  cube.wcs.wcs.restfrq*u.Hz))
cubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                              cube.wcs.wcs.restfrq*u.Hz))
med = cubeK.percentile(50, axis=0)
#cubeK.allow_huge_operations=True
#cubeK = cubeK - med

# determine where absorption...
skew = cubeK.apply_numpy_function(scipy.stats.skew, axis=0)

# BAD error estimate
err = cubeK.std(axis=0)
err[:] = 5*u.K
peak = (cubeK).max(axis=0)
nadir = (cubeK).min(axis=0)
#mask = (peak > 200*u.K) & (skew > 0.1)
absorption_mask = (skew < -0.1)
#mask = mask & (~absorption_mask)

pcube = pyspeckit.Cube(cube=cubeK)

if os.path.exists('e2e_multigauss_fits.fits'):
    pcube.load_model_fit('e2e_multigauss_fits.fits', npars=4)
else:
    #vguesses = 62*u.km/u.s
    #widths = np.ones_like(mask)*2.0
    #guesses = np.array([vguesses.value, widths, temguesses, colguesses])
    #guesses = [62, 2.0, 250., 1e16]
    guesses = np.zeros([17, cubeK.shape[1], cubeK.shape[2]])
    guesses[0,:,:] = 62
    guesses[1,:,:] = 2.0
    guesses[2:,absorption_mask] = -100
    guesses[2:,~absorption_mask] = 100
    guesses[16,:,:] = med

    # For laptop
    mask = ((peak>(200*u.K+med)) & (skew > 0.1)) | ((nadir < (med-200*u.K)) & (skew < -0.1))
    print("Fitting {0} points".format(mask.sum()))

    start_point = (43,43)#np.unravel_index(np.nanargmax(peak*mask), peak.shape)

    position_order = 1./peak.value
    position_order[np.isnan(peak)] = np.inf

    sp = pcube.get_spectrum(*start_point)
    sp.plotter()
    sp.specfit(fittype='ch3cn_spw', guesses=guesses[:,43,43])

    pcube.fiteach(fittype='ch3cn_spw', guesses=guesses, integral=False,
                  verbose_level=0, start_from_point=start_point,
                  use_neighbor_as_guess=True, position_order=position_order,
                  fixed=[F]*17,
                  signal_cut=0,
                  maskmap=mask,
                  errmap=err.value, multicore=4)
    pcube.write_fit('e2e_multigauss_fits.fits', clobber=True)
