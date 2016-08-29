import numpy as np
import paths
from astropy.table import Table, Column
from astropy import table
from astropy import units as u
import masscalc

outflow_tbl = Table.read(paths.tpath("outflow_co_photometry.ipac"), format='ascii.ipac')
core_velo_tbl = Table.read(paths.tpath("core_velocities.ipac"), format="ascii.ipac")
core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')


core_phot_tbl.rename_column('name','SourceID')
cores_merge = table.join(core_velo_tbl, core_phot_tbl,)

brightest_line_flux = np.array([cores_merge[y].data for y in ('peak0','peak1','peak2','peak3')])
peak_line_flux = np.nanmax(brightest_line_flux, axis=0)
peak_line_id = np.nanargmax(brightest_line_flux, axis=0)
brightest_line_name = np.array([cores_merge[('peak0species','peak1species','peak2species','peak3species')[y]][ii] for ii,y in enumerate(peak_line_id)])
cores_merge.add_column(Column(peak_line_flux, name='PeakLineFlux', unit=cores_merge['peak0'].unit))
cores_merge.add_column(Column(brightest_line_name, name='PeakLineSpecies'))

peak_line_brightness = (peak_line_flux*u.Jy).to(u.K, u.brightness_temperature(cores_merge['beam_area'], 220*u.GHz))
cores_merge.add_column(Column(peak_line_brightness, name='PeakLineBrightness'))

temperature_corrected_mass = Column([(masscalc.mass_conversion_factor(20).value
                                      if np.isnan(row['PeakLineBrightness'])
                                      else masscalc.mass_conversion_factor(row['PeakLineBrightness'])).value
                                     * row['peak'] for row in cores_merge],
                                    name='T_corrected_mass',
                                    unit=u.M_sun)
cores_merge.add_column(temperature_corrected_mass)

cores_merge = cores_merge['SourceID',
                          'peak_mass',
                          'T_corrected_mass',
                          'PeakLineBrightness',
                          'PeakLineFlux',
                          'PeakLineSpecies',
                          'peak_col',
                          'mean_velo',
                          'continuum20pct0',
                          'continuum20pct1',
                          'continuum20pct2',
                          'continuum20pct3',
                          'peak',
                          'sum',
                          'peak0',
                          'peak0velo',
                          'peak0freq',
                          'peak0species',
                          'peak1',
                          'peak1velo',
                          'peak1freq',
                          'peak1species',
                          'peak2',
                          'peak2velo',
                          'peak2freq',
                          'peak2species',
                          'peak3',
                          'peak3velo',
                          'peak3freq',
                          'peak3species',
                          'npix',
                          'beam_area',
                          'RA',
                          'Dec',
                          'cont_flux0p2arcsec',
                          'KUbandcont_flux0p2arcsec',
                          'cont_flux0p4arcsec',
                          'KUbandcont_flux0p4arcsec',
                          'cont_flux0p6arcsec',
                          'KUbandcont_flux0p6arcsec',
                          'cont_flux0p8arcsec',
                          'KUbandcont_flux0p8arcsec',
                          'cont_flux1p0arcsec',
                          'KUbandcont_flux1p0arcsec',
                          'cont_flux1p5arcsec',
                          'KUbandcont_flux1p5arcsec',
                         ]
                         

cores_merge.write(paths.tpath('core_continuum_and_line.ipac'), format='ascii.ipac')


### Add columns to the outflow table from the core table ###
newcol = Column([core_phot_tbl['peak_mass'][core_phot_tbl['SourceID'] == name][0]
                 if any(core_phot_tbl['SourceID'] == name) else np.nan
                 for name in outflow_tbl['SourceID']],
                name='CoreMass')
outflow_tbl.add_column(newcol)

newcol = Column([cores_merge['T_corrected_mass'][cores_merge['SourceID'] == name][0]
                 if any(cores_merge['SourceID'] == name) else np.nan
                 for name in outflow_tbl['SourceID']],
                name='TCorrectedCoreMass')
outflow_tbl.add_column(newcol)

#should already be in the file
# newcol = Column([cores_merge['mean_velo'][cores_merge['SourceID'] == name][0]
#                  if any(cores_merge['SourceID'] == name) else np.nan
#                  for name in outflow_tbl['SourceID']],
#                 name='CoreVelocity')
# outflow_tbl.add_column(newcol)

outflow_tbl.write(paths.tpath('outflows_with_cores.ipac'), format='ascii.ipac')

# exec other merge now
with open(paths.apath('merge_spectral_fits_with_photometry.py')) as source_file:
    exec(source_file.read(), globals(), locals())
