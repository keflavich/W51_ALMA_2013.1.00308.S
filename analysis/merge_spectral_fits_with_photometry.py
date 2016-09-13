"""
Generally meant to be run from merge_outflow_core_tables.py
"""
import numpy as np
import paths
from astropy.table import Table, Column
from astropy import table
from astropy import units as u
import masscalc
import scipy.special

spectral_line_fit_tbl = Table.read(paths.tpath('spectral_lines_and_fits.csv'))

cores_merge = Table.read(paths.tpath('core_continuum_and_line.ipac'), format='ascii.ipac')
ppbeam = cores_merge.meta['keywords']['ppbeam_mm']['value']

molcld_exclude_names = ['13COv=0', 'C18O', 'H2CO', 'COv=0']
molcld_exclude = np.array([any(nm in row['Species'] for nm in molcld_exclude_names)
                           for row in spectral_line_fit_tbl])

brightest_noncld_lines = []
brightest_noncld_qns = []
brightest_noncld_fluxes = []
brightest_fitted_brightness = [] 

for row in cores_merge:
    src = row['SourceID']

    amplitudes = spectral_line_fit_tbl['{0}FittedAmplitude'.format(src)]
    noncld_amplitudes = amplitudes * ~molcld_exclude

    brightest_ind = np.argmax(amplitudes)
    brightest_noncloud_ind = np.argmax(noncld_amplitudes)
    if noncld_amplitudes.max() > 0:
        brightest_noncld_lines.append(spectral_line_fit_tbl['Species'][brightest_noncloud_ind])
        brightest_noncld_qns.append(spectral_line_fit_tbl['Resolved QNs'][brightest_noncloud_ind])
        brightest_noncld_fluxes.append(amplitudes[brightest_noncloud_ind])
        brightest_fitted_brightness.append(amplitudes[brightest_noncloud_ind]*spectral_line_fit_tbl['{0}JtoK'.format(src)][brightest_noncloud_ind])
    else:
        brightest_noncld_lines.append('-')
        brightest_noncld_qns.append('-')
        brightest_noncld_fluxes.append(np.nan)
        brightest_fitted_brightness.append(np.nan)


cores_merge.add_column(Column(brightest_noncld_lines, 'BrightestFittedLine'))
cores_merge.add_column(Column(brightest_noncld_qns, 'BrightestFittedQNs'))
cores_merge.add_column(Column(brightest_noncld_fluxes, 'BrightestFittedPeakPixFlux', unit=u.Jy))
cores_merge.add_column(Column(brightest_fitted_brightness, 'BrightestFittedPeakPixBrightness', unit=u.K))

# need to add back in continuum because we're concerned with the *absolute*
# brightness
# moved to other merge jtok_eq = u.brightness_temperature(cores_merge['beam_area'], 225*u.GHz)
# moved to other merge cont_brightness = (u.beam * cores_merge['sum']/cores_merge['npix']).to(u.K, jtok_eq)
cont_brightness = cores_merge['MeanContinuumBrightness']
contincluded_line_brightness = cores_merge['BrightestFittedPeakPixBrightness'] + cont_brightness
cores_merge.add_column(Column(contincluded_line_brightness, 'BrightestFittedPeakPixBrightnessWithcont', unit=u.K))

temperature_corrected_aperturemass = Column([(masscalc.mass_conversion_factor(20).value
                                          if np.isnan(row['BrightestFittedPeakPixBrightnessWithcont'])
                                          or (row['BrightestFittedPeakPixBrightnessWithcont'] < 20)
                                          else masscalc.mass_conversion_factor(row['BrightestFittedPeakPixBrightnessWithcont']).value)
                                         * row['sum']/ppbeam for row in cores_merge],
                                        name='T_corrected_aperturemass',
                                        unit=u.M_sun)
cores_merge.add_column(temperature_corrected_aperturemass)

apertures = ('0p2', '0p4', '0p6', '0p8', '1p0', '1p5')
for ap in apertures:
    tcm = Column([(masscalc.mass_conversion_factor(20).value
                   if np.isnan(row['BrightestFittedPeakPixBrightnessWithcont'])
                   or row['BrightestFittedPeakPixBrightnessWithcont'] < 20 
                   else
                   masscalc.mass_conversion_factor(row['BrightestFittedPeakPixBrightnessWithcont']).value)
                  # already ppbeam corrected
                  * row['cont_flux{0}arcsec'.format(ap)] for row in
                  cores_merge],
                 name='T_corrected_{0}aperturemass'.format(ap), unit=u.M_sun)
    cores_merge.add_column(tcm)


# classify each object based on some fixed criteria
aperture = '0p2'
freefreedominated = (cores_merge['KUbandcont_flux{0}arcsec'.format(aperture)] /
                     cores_merge['cont_flux{0}arcsec'.format(aperture)]) > 0.5
freefreecontaminated = (cores_merge['KUbandcont_flux{0}arcsec'.format(aperture)] /
                        cores_merge['cont_flux{0}arcsec'.format(aperture)]) > 0.1

hot = (cores_merge['BrightestFittedPeakPixBrightnessWithcont'] > 50)
cold = (cores_merge['BrightestFittedPeakPixBrightnessWithcont'] < 20)

# "concentration parameter" is the ratio of the smallest 1-fwhm aperture to the
# annulus around it, divided by 3 because the area of that annulus is 3x the
# area of the inner aperture
concentration_parameter = (cores_merge['cont_flux0p2arcsec'] /
                           ((cores_merge['cont_flux0p4arcsec'] -
                             cores_merge['cont_flux0p2arcsec'])/3.))
cores_merge.add_column(Column(data=concentration_parameter, name='ConcentrationParameter'))
# basing the "threshold" on ALMAmm4 (concentrated, compact) vs ALMAmm6 (no
# central source at all)
compact = concentration_parameter > 2

# approximate beam FHWM
beam_fwhm = 0.2/(8*np.log(2))**0.5
# apparently the integral of a 2D gaussian is approximately
# int_-x^+x e^(-x^2/ (2*sigma**2)) = 2 pi sigma^2 * erf(x * pi / 2.)**2
# or, in terms of FWHM,
# int_-x^+x e^(-x^2/ (2*fwhm**2/2.35**2)) = 2 pi sigma^2 * erf(x * pi / 2.)**2
# if x is 0.2" and the fwhm is 0.2", x = 2.35 sigma
gaussian_0p2 = scipy.special.erf(0.2/beam_fwhm * (2./np.pi))**2
gaussian_0p4 = scipy.special.erf(0.4/beam_fwhm * (2./np.pi))**2
gaussian_0p6 = scipy.special.erf(0.6/beam_fwhm * (2./np.pi))**2

print("Concentration of an unresolved source: {0}"
      .format(gaussian_0p2/(gaussian_0p4-gaussian_0p2)))

cores_merge.add_column(Column(data=compact.astype('int8'),
                              name='is_compact'))
cores_merge.add_column(Column(data=hot.astype('int8'),
                              name='is_hot'))
cores_merge.add_column(Column(data=cold.astype('int8'),
                              name='is_cold'))
cores_merge.add_column(Column(data=freefreedominated.astype('int8'),
                              name='is_freefree'))
                              

category = [('F' if ffd else 'f' if ffc else '_') +
            ('H' if ht else 'C' if cld else '_') +
            ('c' if comp else '_')
            for ffd, ffc, ht, cld, comp in zip(freefreedominated,
                                               freefreecontaminated, hot, cold,
                                               compact)]

cores_merge.add_column(Column(data=category, name='Categories'))

classification = ['HII' if ffd else 
                  'DustyHII' if ffc else
                  'StarlessCore' if (cld and comp) else
                  'HotCore' if (ht and comp) else
                  'ExtendedHotCore' if ht else
                  'ExtendedColdCore' if cld else
                  'UncertainCompact' if comp else
                  'UncertainExtended'
                  for ffd, ffc, ht, cld, comp in zip(freefreedominated,
                                                     freefreecontaminated, hot, cold,
                                                     compact)]
                  
cores_merge.add_column(Column(data=classification, name='Classification'))

for category in np.unique(classification):
    print("{0}: {1}".format(category, np.count_nonzero(np.array(classification)
                                                       == category)))


cores_merge.write(paths.tpath('core_continuum_and_line.ipac'),
                  format='ascii.ipac', overwrite=True)
