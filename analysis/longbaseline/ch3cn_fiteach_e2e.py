from ch3cn_fits import SpectralCube, pyspeckit, fits, u, np
import os
T=True
F=False

cubefn = '../FITS/longbaseline/W51e2e_CH3CN_cutout.fits'
cube = SpectralCube.read(cubefn).minimal_subcube()
contcubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                                  cube.wcs.wcs.restfrq*u.Hz))
cubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                              cube.wcs.wcs.restfrq*u.Hz))
med = cubeK.percentile(25, axis=0)
cubeK.allow_huge_operations=True
cubeK = cubeK - med

# BAD error estimate
err = cubeK.std(axis=0)
err[:] = 5*u.K
peak = (cubeK).max(axis=0)
mask = (peak > 200*u.K)# & (peak > 6*err)
absorption_mask = cubeK.min(axis=0) < -150*u.K
mask = mask & (~absorption_mask)

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
              limitedmin=[T,T,T,T],
              maxpars=[100,5,1500,1e18],
              minpars=[0,0.1,50,1e13],
              signal_cut=0,
              maskmap=mask,
              errmap=err.value, multicore=4)
pcube.write_fit('e2e_CH3CN_Emission_fits.fits', clobber=True)

min_background = 100
background_guess = med.value
background_guess[background_guess < min_background] = min_background
guesses = np.empty((5,)+cube.shape[1:], dtype='float')
guesses[0,:] = 62
guesses[1,:] = 2
guesses[2,:] = 250.
guesses[3,:] = 1e16
guesses[4,:] = background_guess

# again, try cropping out the k=0,1 lines under the assumption that they do not
# trace the disk
pcube_cont = pyspeckit.Cube(cube=contcubeK[:400,:,:])
start_point = (71,66)#np.unravel_index(np.nanargmax(peak*mask), peak.shape)
sp = pcube_cont.get_spectrum(71,66)
sp.specfit(fittype='ch3cn_absorption', guesses=guesses[:,66,71],
           limitedmax=[T,T,T,T,T], limitedmin=[T,T,T,T,T],
           maxpars=[100,5,1500,1e18,10000],
           minpars=[0,0.1,50,1e13,100],)

pcube_cont.fiteach(fittype='ch3cn_absorption', guesses=guesses, integral=False,
                   verbose_level=3, start_from_point=start_point,
                   use_neighbor_as_guess=True, position_order=position_order,
                   limitedmax=[T,T,T,T,T],
                   limitedmin=[T,T,T,T,T],
                   maxpars=[70,5,1500,1e18,10000],
                   minpars=[40,0.1,50,1e13,min_background],
                   signal_cut=0,
                   maskmap=absorption_mask,
                   errmap=err.value, multicore=4)
pcube_cont.write_fit('e2e_CH3CN_Absorption_fits.fits', clobber=True)

from kinematic_analysis_pv_LB import diskycoords, outflowpath
import pvextractor
import pylab as pl

pl.figure(5).clf()
for width in (None, 0.05*u.arcsec, 0.1*u.arcsec, 0.15*u.arcsec):
    diskypath = pvextractor.Path(diskycoords, width)
    extracted_disky = pvextractor.extract_pv_slice(pcube_cont.parcube[0:1,:,:], diskypath, wcs=pcube_cont.wcs)

    pl.plot(extracted_disky.data.squeeze(), label=str(width))

pl.xlabel("Offset (pixels)")
pl.ylabel("Velocity (km/s)")
pl.ylim(55,60)
pl.xlim(855,890)

pl.legend(loc='best')

pl.figure(6).clf()
diskypath = pvextractor.Path(diskycoords, width=None)
extracted_disky = pvextractor.extract_pv_slice(pcube_cont.parcube[0:1,:,:], diskypath, wcs=pcube_cont.wcs)
extracted_disky_width = pvextractor.extract_pv_slice(np.abs(pcube_cont.parcube[1:2,:,:]), diskypath, wcs=pcube_cont.wcs)

pl.fill_between(np.arange(len(extracted_disky.data.squeeze())),
                extracted_disky.data.squeeze()-extracted_disky_width.data.squeeze(),
                extracted_disky.data.squeeze()+extracted_disky_width.data.squeeze(),
                alpha=0.5)
pl.plot(extracted_disky.data.squeeze())


pl.xlabel("Offset (pixels)")
pl.ylabel("Velocity (km/s)")
pl.ylim(52,62)
pl.xlim(855,890)

