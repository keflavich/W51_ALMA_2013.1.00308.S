from ch3cn_fits import SpectralCube, pyspeckit, fits, u, np
import os
import paths
T=True
F=False

cubefn = paths.dpath('longbaseline/W51northcax.SPW2_ALL_cutout.fits')
cubeKfn = paths.dpath('longbaseline/W51northcax.SPW2_ALL_cutout_medsub_K.fits')
medKfn = paths.dpath('longbaseline/W51northcax.SPW2_ALL_cutout_med_K.fits')
if not os.path.exists(cubeKfn):
    cube = SpectralCube.read(cubefn).minimal_subcube()
    contcubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                                      cube.wcs.wcs.restfrq*u.Hz))
    cube.allow_huge_operations = True
    cubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                                  cube.wcs.wcs.restfrq*u.Hz))
    med = cubeK.percentile(25, axis=0)
    cubeK.allow_huge_operations=True
    cubeK = cubeK - med
    med.write(medKfn)
    cubeK.write(cubeKfn)
else:
    cubeK = SpectralCube.read(cubeKfn)
    cubeK.allow_huge_operations = True
    med = fits.getdata(medKfn) * u.K
    contcubeK = cubeK + med
    contcubeK.allow_huge_operations = True

# BAD error estimate
err = cubeK.std(axis=0)
err[:] = 5*u.K
peak = (cubeK).max(axis=0)
mask = (peak > 200*u.K)# & (peak > 6*err)
absorption_mask = cubeK.min(axis=0) < -150*u.K
mask = mask & (~absorption_mask)


min_background = 100
background_guess = med.value
background_guess[background_guess < min_background] = min_background
guesses = np.empty((4,)+cubeK.shape[1:], dtype='float')
guesses[0,:] = background_guess
guesses[1,:] = -1
guesses[2,:] = 61
guesses[3,:] = 1.5

vcontcube_K7 = contcubeK.with_spectral_unit(u.km/u.s,
                                            rest_value=220.53932*u.GHz,
                                            velocity_convention='radio').spectral_slab(50*u.km/u.s,
                                                                                       72*u.km/u.s)
pcube_cont_K7 = pyspeckit.Cube(cube=vcontcube_K7)
start_point = (302,341) # np.unravel_index(np.nanargmax(peak*mask), peak.shape)
sp = pcube_cont_K7.get_spectrum(start_point[0], start_point[1])
sp.plotter()
sp.specfit(fittype='vheightgaussian', guesses=guesses[:,302,341],
           limitedmax=[T,T,T,T,T], limitedmin=[T,T,T,T,T],
           maxpars=[500,0,66,5],
           minpars=[0,-500,55,0],)

k7fitfn = 'north_CH3CN_K7_Gaussian_Absorption_fits'
if os.path.exists(k7fitfn):
    pcube_cont_K7.load_model_fit(k7fitfn, npars=4)
else:
    pcube_cont_K7.fiteach(fittype='vheightgaussian', guesses=guesses, integral=False,
                       verbose_level=3, start_from_point=start_point,
                       use_neighbor_as_guess=True,
                       limitedmax=[T,T,T,T,T],
                       limitedmin=[T,T,T,T,T],
                       maxpars=[500,0,66,5],
                       minpars=[0,-500,55,0],
                       signal_cut=0,
                       maskmap=absorption_mask,
                       errmap=err.value, multicore=4)
    pcube_cont_K7.write_fit(k7fitfn, clobber=True)

from astropy import coordinates
north = coordinates.SkyCoord("19:23:40.052", "14:31:05.50", frame='fk5', unit=(u.hour, u.deg))
pcube_cont_K7.show_fit_param(2,vmin=56,vmax=63)
pcube_cont_K7.mapplot.FITSFigure.recenter(north.ra.deg, north.dec.deg, 0.3/3600.)





vcontcube_K3 = contcubeK.with_spectral_unit(u.km/u.s,
                                            rest_value=220.70902*u.GHz,
                                            velocity_convention='radio').spectral_slab(50*u.km/u.s,
                                                                                       72*u.km/u.s)
pcube_cont_K3 = pyspeckit.Cube(cube=vcontcube_K3)
start_point = (302,341) # np.unravel_index(np.nanargmax(peak*mask), peak.shape)
sp = pcube_cont_K3.get_spectrum(start_point[0], start_point[1])
sp.plotter()
sp.specfit(fittype='vheightgaussian', guesses=guesses[:,302,341],
           limitedmax=[T,T,T,T,T], limitedmin=[T,T,T,T,T],
           maxpars=[500,0,66,5],
           minpars=[0,-500,55,0],)

k3fitfn = 'north_CH3CN_K3_Gaussian_Absorption_fits'
if os.path.exists(k3fitfn):
    pcube_cont_K3.load_model_fit(k3fitfn, npars=4)
else:
    pcube_cont_K3.fiteach(fittype='vheightgaussian', guesses=guesses, integral=False,
                       verbose_level=3, start_from_point=start_point,
                       use_neighbor_as_guess=True,
                       limitedmax=[T,T,T,T,T],
                       limitedmin=[T,T,T,T,T],
                       maxpars=[500,0,66,5],
                       minpars=[0,-500,55,0],
                       signal_cut=0,
                       maskmap=absorption_mask,
                       errmap=err.value, multicore=4)
    pcube_cont_K3.write_fit(k3fitfn, clobber=True)

from astropy import coordinates
north = coordinates.SkyCoord("19:23:40.052", "14:31:05.50", frame='fk5', unit=(u.hour, u.deg))
pcube_cont_K3.show_fit_param(2,vmin=56,vmax=63)
pcube_cont_K3.mapplot.FITSFigure.recenter(north.ra.deg, north.dec.deg, 0.3/3600.)
