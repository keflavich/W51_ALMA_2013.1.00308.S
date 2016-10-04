from astropy.table import Table, Column
from astropy import units as u
from latex_info import latexdict
from astroquery.splatalogue import Splatalogue

from line_to_image_list import line_to_image_list, labeldict

for spwid in (0,1,2,3):
    label_col = Column(data=[labeldict[name] for name,_,_,spw in line_to_image_list
                             if spw==spwid],
                       name='Line Name')
    freq_col = Column(data=[float(freq.strip('GHz')) for _,freq,_,spw in line_to_image_list
                            if spw==spwid],
                      name='Frequency',
                      unit=u.GHz,)

    EU = [set(Splatalogue.query_lines(freq*(1-0.001/3e5),
                                      freq*(1+0.001/3e5))['E_U (K)'])
          for freq in u.Quantity(freq_col)]

    #spw_col = Column(data=[spw for _,_,_,spw in line_to_image_list
    #                       if spw==spwid],
    #                 name='Spectral Window',)

    tbl = Table([label_col, freq_col,])# spw_col])


    latexdict['header_start'] = '\label{{tab:linesspw{0}}}'.format(spwid)
    latexdict['caption'] = 'Spectral Lines in SPW {0}'.format(spwid)
    #latexdict['tablefoot'] = ('\par\nJy-Kelvin gives the conversion factor from Jy '
    #                          'to Kelvin given the synthesized beam size and '
    #                          'observation frequency.')
    latexdict['col_align'] = 'l'*len(tbl.columns)
    #latexdict['tabletype'] = 'longtable'
    #latexdict['tabulartype'] = 'longtable'
    latexdict['units'] = {}

    tbl.write('../cores_and_outflows/linetable{0}.tex'.format(spwid),
              latexdict=latexdict, format='ascii.latex', overwrite=True)
