import numpy as np
import paths
from astropy.table import Table, Column
from astropy import table
from astropy import units as u
from astropy import coordinates
import masscalc

dendro_velo_tbl = Table.read(paths.tpath("dendro_core_velocities.ipac"), format="ascii.ipac")
dendro_phot_tbl = Table.read(paths.tpath("dendrogram_continuum_catalog.ipac"), format='ascii.ipac')

dendro_phot_tbl.rename_column('_idx','SourceID')
dendro_merge = table.join(dendro_velo_tbl, dendro_phot_tbl,)

brightest_line_flux = np.array([dendro_merge[y].data for y in
                                ('peak0','peak1','peak2','peak3')])
peak_line_flux = np.nanmax(brightest_line_flux, axis=0)
peak_line_id = np.nanargmax(np.nan_to_num(brightest_line_flux), axis=0)
peak_line_id[np.isnan(peak_line_flux)] = -999
brightest_line_name = np.array([dendro_merge[('peak0species', 'peak1species',
                                              'peak2species',
                                              'peak3species')[y]][ii]
                                if y >= 0 else 'NONE'
                                for ii,y in enumerate(peak_line_id)])
dendro_merge.add_column(Column(peak_line_flux, name='PeakLineFlux',
                               unit=dendro_merge['peak0'].unit))
dendro_merge.add_column(Column(brightest_line_name, name='PeakLineSpecies'))

peak_line_beam_area = u.Quantity([dendro_merge['beam{0}area'.format(pid)][ii]
                                  if pid >= 0 else 1.0
                                  for ii,pid in enumerate(peak_line_id)],
                                 u.sr)
peak_line_freq = u.Quantity([dendro_merge['peak{0}freq'.format(pid)][ii]
                             if pid >= 0 else 1.0
                             for ii,pid in enumerate(peak_line_id)],
                            u.GHz)
peak_line_cont = u.Quantity([dendro_merge['continuum20pct{0}'.format(pid)][ii]
                             if pid >= 0 else 0.0
                             for ii,pid in enumerate(peak_line_id)],
                            u.Jy)
tbequiv = u.brightness_temperature(peak_line_beam_area, peak_line_freq)
peak_line_brightness = (peak_line_flux*u.Jy).to(u.K, tbequiv)
dendro_merge.add_column(Column(peak_line_brightness, name='PeakLineBrightness'))
continuum20pct_K = (peak_line_cont.to(u.K, tbequiv))
dendro_merge.add_column(Column(continuum20pct_K, name='PeakLineContinuumBG'))

tcm = Column([(masscalc.mass_conversion_factor(20).to(u.M_sun).value if
               np.isnan(row['PeakLineBrightness']) else
               masscalc.mass_conversion_factor(TK=row['PeakLineBrightness']).to(u.M_sun).value)
              * row['peak_cont_flux'] for row in dendro_merge],
             name='T_corrected_mass', unit=u.M_sun)
temperature_corrected_mass = tcm
dendro_merge.add_column(temperature_corrected_mass)

columns = ['SourceID',
           'peak_cont_mass',
           'T_corrected_mass',
           'PeakLineBrightness',
           'PeakLineFlux',
           'PeakLineSpecies',
           'peak_cont_col',
           'mean_velo',
           'peak_cont_flux',
          ]
columns += sorted([c for c in dendro_merge.colnames if c not in columns])

dendro_merge = dendro_merge[columns]

cat = coordinates.SkyCoord(u.Quantity(dendro_merge['x_cen'].data, u.deg), u.Quantity(dendro_merge['y_cen'].data, u.deg), frame='fk5')
nearest = cat.match_to_catalog_sky(cat, 2)
dendro_merge.add_column(Column(nearest[1], name='nndist'))

dendro_merge.write(paths.tpath('dendro_merge_continuum_and_line.ipac'), format='ascii.ipac')
