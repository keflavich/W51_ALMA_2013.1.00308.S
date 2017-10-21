from ch3cn_fits import (SpectralCube, pyspeckit, fits, u, np, line_name_dict,
                        line_aij, line_deg, frequencies, line_names, line_eu,
                        fit_tex, nupper_of_kkms,
                        fit_all_tex)
import pyregion
from astropy import constants
import scipy.stats
import os
import paths
from astropy import coordinates
from astropy import units as u
T=True
F=False

e2e_peak = coordinates.SkyCoord(290.9331916067147*u.deg, 14.50958361099483*u.deg, frame='icrs')

spectra = []
for fn in ('W51e2cax.SPW{0}_ALL_cutout.fits'.format(ii) for ii in range(10)):

    cubefn = paths.dpath('longbaseline/linked/{0}'.format(fn))
    cube = SpectralCube.read(cubefn)
    #contcubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
    #                                                  cube.wcs.wcs.restfrq*u.Hz))
    #cubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
    #                                              cube.wcs.wcs.restfrq*u.Hz))


    xx,yy = cube.wcs.celestial.wcs_world2pix(e2e_peak.ra.deg, e2e_peak.dec.deg, 0)
    spectrum_Jy = cube[:, int(np.round(yy)), int(np.round(xx))]
    spectrum_K = spectrum_Jy.to(u.K, cube.beam.jtok_equiv(cube.spectral_axis))

    spectra.append(spectrum_K)


# already done at import stage
#pyspeckit.spectrum.fitters.default_Registry.add_fitter('ch3cn_absorption',ch3cn_absorption_fitter(),5)

guesses = [58.205, 2.25, 40, 5e16, 500,
           56.205, 2.25, 250, 5e16, 500,
          ]

sp = pyspeckit.Spectra([pyspeckit.Spectrum.from_hdu(s.hdu) for s in spectra
                        if ((s.spectral_axis[0] > 220*u.GHz) and
                            (s.spectral_axis[0] < 222*u.GHz))])
# we have to zero the data before fitting for this model...
sp.data -= 500
sp.plotter()
sp.specfit(fittype='ch3cn_absorption', guesses=guesses,
           limitedmax=[T,T,T,T,T]*2, limitedmin=[T,T,T,T,T]*2,
           fixed=[F,F,F,F,T]*2,
           maxpars=[100,5,1500,1e18,10000]*2,
           minpars=[0,0.1,20,1e13,100]*2,)

#    #med = cubeK.percentile(50, axis=0)
#
#    ## determine where absorption...
#    #skew = cubeK.apply_numpy_function(scipy.stats.skew, axis=0)
#
#    # BAD error estimate
#    #err = cubeK.std(axis=0)
#    #err[:] = 5*u.K
#    #peak = (cubeK).max(axis=0)
#    #nadir = (cubeK).min(axis=0)
#    #mask = (peak > 200*u.K) & (skew > 0.1)
#    #absorption_mask = (skew < -0.1)
#    #mask = mask & (~absorption_mask)
#
#    min_background = 100
#    background_guess = med.value
#    background_guess[background_guess < min_background] = min_background
#    guesses = np.empty((4,)+cubeK.shape[1:], dtype='float')
#    guesses[0,:] = background_guess
#    guesses[1,:] = -100
#    guesses[2,:] = 56
#    guesses[3,:] = 1.5
#
#    vcontcube_K7 = contcubeK.with_spectral_unit(u.km/u.s,
#                                                rest_value=220.53932*u.GHz,
#                                                velocity_convention='radio').spectral_slab(42*u.km/u.s,
#                                                                                           72*u.km/u.s)
#    pcube_cont_K7 = pyspeckit.Cube(cube=vcontcube_K7)
#    start_point = (43,43) # np.unravel_index(np.nanargmax(peak*mask), peak.shape)
#    sp = pcube_cont_K7.get_spectrum(start_point[0], start_point[1])
#    sp.plotter()
#    sp.specfit(fittype='vheightgaussian', guesses=guesses[:,43,43],
#               limitedmax=[T,T,T,T,T], limitedmin=[T,T,T,T,T],
#               maxpars=[500,0,61,6],
#               minpars=[0,-500,52,0],)
#
#    k7fitfn = 'e2e_CH3CN_K7_Gaussian_Absorption_fits'
#    if os.path.exists(k7fitfn):
#        pcube_cont_K7.load_model_fit(k7fitfn, npars=4)
#    else:
#        pcube_cont_K7.fiteach(fittype='vheightgaussian', guesses=guesses, integral=False,
#                           verbose_level=3, start_from_point=start_point,
#                           use_neighbor_as_guess=True,
#                           limitedmax=[T,T,T,T,T],
#                           limitedmin=[T,T,T,T,T],
#                           maxpars=[500,0,61,6],
#                           minpars=[0,-500,52,0],
#                           signal_cut=0,
#                           maskmap=absorption_mask,
#                           errmap=err.value, multicore=4)
#        pcube_cont_K7.write_fit(k7fitfn, clobber=True)
#
#    from astropy import coordinates
#    #e2e = coordinates.SkyCoord("19:23:43.939", "14:30:34.57", frame='fk5', unit=(u.hour, u.deg))
#    pcube_cont_K7.show_fit_param(2,vmin=50,vmax=63)
#    #pcube_cont_K7.mapplot.FITSFigure.recenter(e2e.ra.deg, e2e.dec.deg, 0.3/3600.)
#
#
#
#
#    guesses[0,:] = 1
#    guesses[1,:] = -1
#    guesses[2,:] = 56
#    guesses[3,:] = 1.5
#
#    cube_k8 = SpectralCube.read(paths.dpath('longbaseline/velo_cutouts/e2e_CH3CNv_0_12_8_-11_8_.fits'))
#    cube_k8 = cube_k8.subcube_from_ds9region(pyregion.open(paths.rpath('w51e2box_ch3cn.reg')))
#
#    pcube_cont_K8 = pyspeckit.Cube(cube=cube_k8)
#    start_point = (43,43) # np.unravel_index(np.nanargmax(peak*mask), peak.shape)
#    sp = pcube_cont_K8.get_spectrum(start_point[0], start_point[1])
#    sp.plotter()
#    sp.specfit(fittype='vheightgaussian', guesses=guesses[:,43,43],
#               limitedmax=[T,T,T,T,T], limitedmin=[T,T,T,T,T],
#               maxpars=[5,0,61,6],
#               minpars=[0,-5,52,0],)
#
#    k8fitfn = 'e2e_CH3CN_K8_Gaussian_Absorption_fits'
#    if os.path.exists(k8fitfn):
#        pcube_cont_K8.load_model_fit(k8fitfn, npars=4)
#    else:
#        pcube_cont_K8.fiteach(fittype='vheightgaussian', guesses=guesses, integral=False,
#                           verbose_level=3, start_from_point=start_point,
#                           use_neighbor_as_guess=True,
#                           limitedmax=[T,T,T,T,T],
#                           limitedmin=[T,T,T,T,T],
#                           maxpars=[5,0,61,6],
#                           minpars=[0,-5,52,0],
#                           signal_cut=0,
#                           maskmap=absorption_mask,
#                           errmap=err.value, multicore=4)
#        pcube_cont_K8.write_fit(k8fitfn, clobber=True)
#
#    #e2e = coordinates.SkyCoord("19:23:43.939", "14:30:34.57", frame='fk5', unit=(u.hour, u.deg))
#    pcube_cont_K8.show_fit_param(2,vmin=50,vmax=63)
#    #pcube_cont_K8.mapplot.FITSFigure.recenter(e2e.ra.deg, e2e.dec.deg, 0.3/3600.)
#
#
#
#
#
#
#    guesses[0,:] = background_guess
#    guesses[1,:] = -100
#    guesses[2,:] = 56
#    guesses[3,:] = 1.5
#
#    vcontcube_K3 = contcubeK.with_spectral_unit(u.km/u.s,
#                                                rest_value=220.70902*u.GHz,
#                                                velocity_convention='radio').spectral_slab(42*u.km/u.s,
#                                                                                           72*u.km/u.s)
#    pcube_cont_K3 = pyspeckit.Cube(cube=vcontcube_K3)
#    start_point = (43,43) # np.unravel_index(np.nanargmax(peak*mask), peak.shape)
#    sp = pcube_cont_K3.get_spectrum(start_point[0], start_point[1])
#    sp.plotter()
#    sp.specfit(fittype='vheightgaussian', guesses=guesses[:,43,43],
#               limitedmax=[T,T,T,T,T], limitedmin=[T,T,T,T,T],
#               maxpars=[500,0,66,5],
#               minpars=[0,-500,55,0],)
#
#    k3fitfn = 'e2e_CH3CN_K3_Gaussian_Absorption_fits'
#    if os.path.exists(k3fitfn):
#        pcube_cont_K3.load_model_fit(k3fitfn, npars=4)
#    else:
#        pcube_cont_K3.fiteach(fittype='vheightgaussian', guesses=guesses, integral=False,
#                           verbose_level=3, start_from_point=start_point,
#                           use_neighbor_as_guess=True,
#                           limitedmax=[T,T,T,T,T],
#                           limitedmin=[T,T,T,T,T],
#                           maxpars=[500,0,66,5],
#                           minpars=[0,-500,55,0],
#                           signal_cut=0,
#                           maskmap=absorption_mask,
#                           errmap=err.value, multicore=4)
#        pcube_cont_K3.write_fit(k3fitfn, clobber=True)
#
#    from astropy import coordinates
#    pcube_cont_K3.show_fit_param(2,vmin=51,vmax=63)
#    #pcube_cont_K3.mapplot.FITSFigure.recenter(e2e.ra.deg, e2e.dec.deg, 0.3/3600.)
#
