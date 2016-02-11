import numpy as np
import paths
from astropy.table import Table, Column

outflow_tbl = Table.read(paths.tpath("outflow_co_photometry.ipac"), format='ascii.ipac')
core_velo_tbl = Table.read(paths.tpath("core_velocities.ipac"), format="ascii.ipac")
core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')

newcol = Column([core_phot_tbl['peak_mass'][core_phot_tbl['name'] == name][0]
                 if any(core_phot_tbl['name'] == name) else np.nan
                 for name in outflow_tbl['SourceID']],
                name='CoreMass')
outflow_tbl.add_column(newcol)
