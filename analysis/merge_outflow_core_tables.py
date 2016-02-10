import paths
from astropy.table import Table

outflow_tbl = Table.read(paths.tpath("outflow_co_photometry.ipac"), format='ascii.ipac')
core_velo_tbl = Table.read(paths.tpath("core_velocities.ipac"), format="ascii.ipac")
core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
