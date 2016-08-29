"""
Generally meant to be run from merge_outflow_core_tables.py
"""
import numpy as np
import paths
from astropy.table import Table, Column
from astropy import table
from astropy import units as u
import masscalc

spectral_line_fit_tbl = Table.read(paths.tpath('spectral_lines_and_fits.csv'))

cores_merge = Table.read(paths.tpath('core_continuum_and_line.ipac'), format='ascii.ipac')

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
cores_merge.add_column(Column(brightest_noncld_fluxes, 'BrightestFittedApMeanFlux', unit=u.Jy))
cores_merge.add_column(Column(brightest_fitted_brightness, 'BrightestFittedApMeanBrightness', unit=u.K))

# need to add back in continuum because we're concerned with the *absolute*
# brightness
# moved to other merge jtok_eq = u.brightness_temperature(cores_merge['beam_area'], 225*u.GHz)
# moved to other merge cont_brightness = (u.beam * cores_merge['sum']/cores_merge['npix']).to(u.K, jtok_eq)
cont_brightness = cores_merge['MeanContinuumBrightness']
contincluded_line_brightness = cores_merge['BrightestFittedApMeanBrightness'] + cont_brightness
cores_merge.add_column(Column(contincluded_line_brightness, 'BrightestFittedApMeanBrightnessWithcont', unit=u.K))

temperature_corrected_aperturemass = Column([(masscalc.mass_conversion_factor(20).value
                                          if np.isnan(row['BrightestFittedApMeanBrightnessWithcont'])
                                          else masscalc.mass_conversion_factor(row['BrightestFittedApMeanBrightnessWithcont']).value)
                                         * row['sum']/ppbeam for row in cores_merge],
                                        name='T_corrected_aperturemass',
                                        unit=u.M_sun)
cores_merge.add_column(temperature_corrected_aperturemass)

cores_merge.write(paths.tpath('core_continuum_and_line.ipac'), format='ascii.ipac', overwrite=True)
