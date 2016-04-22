from ch3cn_fits import SpectralCube, pyspeckit, fits, u, np
import os
T=True
F=False

cubefn = '../FITS/longbaseline/W51e2e_CH3CN_cutout.fits'
cube = SpectralCube.read(cubefn).minimal_subcube()
med = cube.percentile(25, axis=0)
cube.allow_huge_operations=True
cube = cube - med
cubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                              cube.wcs.wcs.restfrq*u.Hz))
# BAD error estimate
err = cubeK.std(axis=0)
err[:] = 5*u.K
peak = cubeK.max(axis=0)
mask = (peak > 200*u.K)# & (peak > 6*err)

pcube = pyspeckit.Cube(cube=cubeK[:400,:,:]) # crop out k=0,1

vguesses = 62*u.km/u.s
colguesses = np.ones_like(mask)*1e16
temguesses = np.ones_like(mask)*250.
widths = np.ones_like(mask)*2.0
#guesses = np.array([vguesses.value, widths, temguesses, colguesses])
guesses = [62, 2.0, 250., 1e16]

# For laptop
#mask &= (peak>10*u.K)

start_point = (43,43)#np.unravel_index(np.nanargmax(peak*mask), peak.shape)

position_order = 1./peak.value
position_order[np.isnan(peak)] = np.inf

sp = pcube.get_spectrum(*start_point)
sp.plotter()
sp.specfit(fittype='ch3cn', guesses=guesses)

pcube.fiteach(fittype='ch3cn', guesses=guesses, integral=False,
              verbose_level=3, start_from_point=start_point,
              use_neighbor_as_guess=True, position_order=position_order,
              limitedmax=[T,T,T,T],
              maxpars=[100,5,1500,1e18],
              minpars=[0,0.1,50,1e13],
              signal_cut=0,
              maskmap=mask,
              errmap=err.value, multicore=4)
