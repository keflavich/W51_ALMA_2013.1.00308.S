from ch3cn_fits import (SpectralCube, pyspeckit, fits, u, np, line_name_dict,
                        line_aij, line_deg, frequencies, line_names, line_eu,
                        fit_tex, nupper_of_kkms,
                        fit_all_tex)
from astropy import constants
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
    guesses[16,:,:] = med * (med > 0)

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
                  limited=[(T,T)]*17,
                  maxpars=[70,6]+[600]*15,
                  minpars=[50,0.1]+[-600]*14+[0],
                  fixed=[F]*17,
                  signal_cut=0,
                  maskmap=mask,
                  errmap=err.value, multicore=4)
    pcube.write_fit('e2e_multigauss_fits.fits', clobber=True)


# do some rotational diagram things...
ch3cn_inds = [ii for ii, name in enumerate(line_names)
              if 'CH3CN' in name]
ch3cn_freqs = np.array(frequencies)[ch3cn_inds]
ch3cn_energy = ([line_eu[name] for name in line_names if 'CH3CN' in name]*u.erg).to(u.K, u.temperature_energy())
ch3cn_degs = [line_deg[name] for name in line_names if 'CH3CN' in name]
ch3cn_aij = [10**line_aij[name] for name in line_names if 'CH3CN' in name]

# CH3CN cube is in K km/s
ch3cn_cube = pcube.parcube[np.array(ch3cn_inds)+2,:,:] * np.sqrt(2*np.pi) * pcube.parcube[1,:,:]
# but for absorption lines, we want the positive component, so we subtract the
# (negative) flux from the background
background_integral = pcube.parcube[-1,:,:] * np.sqrt(2*np.pi) * pcube.parcube[1,:,:]
absorption_mask = ch3cn_cube[0,:,:] < 0
ch3cn_cube[:,absorption_mask] = background_integral[absorption_mask] - ch3cn_cube[:,absorption_mask]

amp_errs = pcube.errcube[np.array(ch3cn_inds)+2,:,:]
sigma_err = pcube.errcube[1,:,:]
amps = pcube.parcube[np.array(ch3cn_inds)+2,:,:]
sigma = pcube.parcube[1,:,:]
ch3cn_errs = (2*np.pi*(sigma*np.abs(amps))**2*((amp_errs/amps)**2 + (sigma_err/sigma)**2))**0.5
ch3cn_pcube = pyspeckit.Cube(cube=ch3cn_cube, error=ch3cn_errs,
                             xarr=ch3cn_energy)

# single spec example
ii = 75
jj = 75
fit_tex(ch3cn_energy, nupper_of_kkms(ch3cn_cube[:,ii,jj], ch3cn_freqs,
                                     ch3cn_aij, ch3cn_degs).value,
        errors=nupper_of_kkms(ch3cn_errs[:,ii,jj], ch3cn_freqs, ch3cn_aij,
                              ch3cn_degs).value,
        plot=True
       )

tmap,Nmap = fit_all_tex(xaxis=ch3cn_energy, cube=ch3cn_cube,
                        cubefrequencies=ch3cn_freqs,
                        degeneracies=ch3cn_degs,
                        einsteinAij=ch3cn_aij,
                        errorcube=ch3cn_errs,)
