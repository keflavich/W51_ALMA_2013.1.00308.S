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

brightest_line_flux = np.array([dendro_merge[y].data for y in ('peak0','peak1','peak2','peak3')])
peak_line_flux = np.nanmax(brightest_line_flux, axis=0)
peak_line_id = np.nanargmax(np.nan_to_num(brightest_line_flux), axis=0)
peak_line_id[np.isnan(peak_line_flux)] = -999
brightest_line_name = np.array([dendro_merge[('peak0species','peak1species','peak2species','peak3species')[y]][ii]
                                if y >= 0 else 'NONE'
                                for ii,y in enumerate(peak_line_id)])
dendro_merge.add_column(Column(peak_line_flux, name='PeakLineFlux', unit=dendro_merge['peak0'].unit))
dendro_merge.add_column(Column(brightest_line_name, name='PeakLineSpecies'))

# WRONG!!!! the beam area is for the continuum data, which is smaller than that for the line data in general
peak_line_brightness = (peak_line_flux*u.Jy).to(u.K, u.brightness_temperature(u.Quantity(np.array(dendro_merge['beam_area']),
                                                                                         u.sr),
                                                                              220*u.GHz))
dendro_merge.add_column(Column(peak_line_brightness, name='PeakLineBrightness'))

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
           'continuum20pct',
           'peak_cont_flux',
          ]
columns += sorted([c for c in dendro_merge.colnames if c not in columns])

dendro_merge = dendro_merge[columns]

cat = coordinates.SkyCoord(u.Quantity(dendro_merge['x_cen'].data, u.deg), u.Quantity(dendro_merge['y_cen'].data, u.deg), frame='fk5')
nearest = cat.match_to_catalog_sky(cat, 2)
dendro_merge.add_column(Column(nearest[1], name='nndist'))

dendro_merge.write(paths.tpath('dendro_merge_continuum_and_line.ipac'), format='ascii.ipac')
