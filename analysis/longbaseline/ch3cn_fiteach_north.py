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

# SKIP EMISSION pcube = pyspeckit.Cube(cube=cubeK[:400,:,:]) # crop out k=0,1
# SKIP EMISSION 
# SKIP EMISSION vguesses = 62*u.km/u.s
# SKIP EMISSION colguesses = np.ones_like(mask)*1e16
# SKIP EMISSION temguesses = np.ones_like(mask)*250.
# SKIP EMISSION widths = np.ones_like(mask)*2.0
# SKIP EMISSION #guesses = np.array([vguesses.value, widths, temguesses, colguesses])
# SKIP EMISSION guesses = [62, 1.0, 250., 1e16, 55, 1.0, 250., 1e16]
# SKIP EMISSION 
# SKIP EMISSION # For laptop
# SKIP EMISSION #mask &= (peak>10*u.K)
# SKIP EMISSION 
# SKIP EMISSION start_point = np.unravel_index(np.nanargmax(peak*mask), peak.shape)
# SKIP EMISSION start_point = (281,306)
# SKIP EMISSION 
# SKIP EMISSION position_order = 1./peak.value
# SKIP EMISSION position_order[np.isnan(peak)] = np.inf
# SKIP EMISSION 
# SKIP EMISSION sp = pcube.get_spectrum(*start_point)
# SKIP EMISSION sp.plotter()
# SKIP EMISSION sp.specfit(fittype='ch3cn', guesses=guesses)
# SKIP EMISSION 
# SKIP EMISSION if os.path.exists('north_CH3CN_Emission_fits.fits'):
# SKIP EMISSION     pcube.load_model_fit('north_CH3CN_Emission_fits', npars=4)
# SKIP EMISSION else:
# SKIP EMISSION     pcube.fiteach(fittype='ch3cn', guesses=guesses, integral=False,
# SKIP EMISSION                   verbose_level=3, start_from_point=start_point,
# SKIP EMISSION                   use_neighbor_as_guess=True, position_order=position_order,
# SKIP EMISSION                   limitedmax=[T,T,T,T],
# SKIP EMISSION                   limitedmin=[T,T,T,T],
# SKIP EMISSION                   maxpars=[100,5,1500,1e18],
# SKIP EMISSION                   minpars=[0,0.1,50,1e13],
# SKIP EMISSION                   signal_cut=0,
# SKIP EMISSION                   maskmap=mask,
# SKIP EMISSION                   errmap=err.value, multicore=4)
# SKIP EMISSION     pcube.write_fit('north_CH3CN_Emission_fits.fits', clobber=True)

min_background = 100
background_guess = med.value
background_guess[background_guess < min_background] = min_background
guesses = np.empty((5,)+cubeK.shape[1:], dtype='float')
guesses[0,:] = 61
guesses[1,:] = 2
guesses[2,:] = 250.
guesses[3,:] = 1e16
guesses[4,:] = background_guess

# again, try cropping out the k=0,1 lines under the assumption that they do not
# trace the disk
pcube_cont = pyspeckit.Cube(cube=contcubeK[:400,:,:])
start_point = (302,341) # np.unravel_index(np.nanargmax(peak*mask), peak.shape)
sp = pcube_cont.get_spectrum(start_point[0], start_point[1])
sp.plotter()
sp.specfit(fittype='ch3cn_absorption', guesses=guesses[:,66,71],
           limitedmax=[T,T,T,T,T], limitedmin=[T,T,T,T,T],
           maxpars=[100,5,1500,1e18,10000],
           minpars=[0,0.1,50,1e13,100],)

if os.path.exsits('north_CH3CN_Absorption_fits'):
    pcube_cont.load_model_fit('north_CH3CN_Absorption_fits', npars=5)
else:
    pcube_cont.fiteach(fittype='ch3cn_absorption', guesses=guesses, integral=False,
                       verbose_level=3, start_from_point=start_point,
                       use_neighbor_as_guess=True,
                       limitedmax=[T,T,T,T,T],
                       limitedmin=[T,T,T,T,T],
                       maxpars=[70,5,1500,1e18,10000],
                       minpars=[40,0.1,50,1e13,min_background],
                       signal_cut=0,
                       maskmap=absorption_mask,
                       errmap=err.value, multicore=4)
    pcube_cont.write_fit('north_CH3CN_Absorption_fits.fits', clobber=True)

from astropy import coordinates
north = coordinates.SkyCoord("19:23:40.052", "14:31:05.50", frame='fk5')
pcube_cont.show_fit_param(0,vmin=58,vmax=63)
pcube_cont.mapplot.FITSFigure.recenter(north.ra.deg, north.dec.deg, 0.3/3600.)
pcube.show_fit_param(0, vmin=60, vmax=70)
pcube.mapplot.FITSFigure.recenter(north.ra.deg, north.dec.deg, 0.3/3600.)

# from kinematic_analysis_pv_LB import diskycoords, outflowpath
# import pvextractor
# import pylab as pl
# 
# pl.figure(5).clf()
# for width in (None, 0.05*u.arcsec, 0.1*u.arcsec, 0.15*u.arcsec):
#     diskypath = pvextractor.Path(diskycoords, width)
#     extracted_disky = pvextractor.extract_pv_slice(pcube_cont.parcube[0:1,:,:], diskypath, wcs=pcube_cont.wcs)
# 
#     pl.plot(extracted_disky.data.squeeze(), label=str(width))
# 
# pl.xlabel("Offset (pixels)")
# pl.ylabel("Velocity (km/s)")
# pl.ylim(55,60)
# pl.xlim(855,890)
# 
# pl.legend(loc='best')
# 
# pl.figure(6).clf()
# diskypath = pvextractor.Path(diskycoords, width=None)
# extracted_disky = pvextractor.extract_pv_slice(pcube_cont.parcube[0:1,:,:], diskypath, wcs=pcube_cont.wcs)
# extracted_disky_width = pvextractor.extract_pv_slice(np.abs(pcube_cont.parcube[1:2,:,:]), diskypath, wcs=pcube_cont.wcs)
# 
# pl.fill_between(np.arange(len(extracted_disky.data.squeeze())),
#                 extracted_disky.data.squeeze()-extracted_disky_width.data.squeeze(),
#                 extracted_disky.data.squeeze()+extracted_disky_width.data.squeeze(),
#                 alpha=0.5)
# pl.plot(extracted_disky.data.squeeze())
# 
# 
# pl.xlabel("Offset (pixels)")
# pl.ylabel("Velocity (km/s)")
# pl.ylim(52,62)
# pl.xlim(855,890)
# 
