import numpy as np
from astropy.table import Table, Column
from astropy import constants
from astropy import units as u
from latex_info import latexdict

tbl = Table.read('spw_table.txt', format='ascii')

frq_range = list(zip(tbl['Ch0(MHz)']*u.MHz, (tbl['Ch0(MHz)'] + tbl['ChanWid(kHz)']/1e3*tbl['#Chans'])*u.MHz))

tbl.add_column(Column(data=u.Quantity([min(x) for x in frq_range], u.GHz), name='Minimum Frequency'))
tbl.add_column(Column(data=u.Quantity([max(x) for x in frq_range], u.GHz), name='Maximum Frequency'))
cwnu = 'Channel Width [$\\nu$]'
cwv = 'Channel Width [$v$]'
tbl.rename_column('ChanWid(kHz)', cwnu)
tbl[cwnu].unit = u.kHz
tbl.add_column(Column(data=np.abs((u.Quantity(tbl[cwnu]) /
                                   u.Quantity(tbl['Minimum Frequency']) *
                                   constants.c).to(u.km/u.s)),
                      name=cwv))

tbl = tbl['SpwID', 'Minimum Frequency', 'Maximum Frequency', cwnu, cwv]

latexdict['header_start'] = '\label{tab:spw}'
latexdict['caption'] = 'Spectral Setup'
#latexdict['tablefoot'] = ('\par\nJy-Kelvin gives the conversion factor from Jy '
#                          'to Kelvin given the synthesized beam size and '
#                          'observation frequency.')
latexdict['col_align'] = 'l'*len(tbl.columns)
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
latexdict['units'] = {}

formats = {cwv: lambda x: '{0:0.2f}'.format(x)}

tbl.write('../paper1/spwtable.tex', latexdict=latexdict,
          format='ascii.latex',
          formats=formats,
          overwrite=True)
