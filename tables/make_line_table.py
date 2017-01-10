from astropy.table import Table, Column
from astropy import units as u
from latex_info import latexdict
from astroquery.splatalogue import Splatalogue

from line_to_image_list import line_to_image_list, labeldict

methanol_rows = []

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

    methanol_rows += [(labeldict[name], float(freq.strip('GHz')),
                       list(set(
                           Splatalogue.query_lines(u.Quantity(freq)*(1-0.01/3e5),
                                                   u.Quantity(freq)*(1+0.01/3e5),
                                                   line_lists=['JPL','SLAIM'],
                                                   chemical_name='Methanol',
                                                  )['E_U (K)']))[0]
                      )
                      for name,freq,_,spw in line_to_image_list
                      if 'CH3OH' in name[:5] and spw==spwid]

    #spw_col = Column(data=[spw for _,_,_,spw in line_to_image_list
    #                       if spw==spwid],
    #                 name='Spectral Window',)

    tbl = Table([label_col, freq_col,])# spw_col])


    latexdict['header_start'] = '\label{{tab:linesspw{0}}}'.format(spwid)
    latexdict['caption'] = 'Spectral Lines in SPW {0}'.format(spwid)
    latexdict['tablefoot'] = ''
    #latexdict['tablefoot'] = ('\par\nJy-Kelvin gives the conversion factor from Jy '
    #                          'to Kelvin given the synthesized beam size and '
    #                          'observation frequency.')
    latexdict['col_align'] = 'l'*len(tbl.columns)
    #latexdict['tabletype'] = 'longtable'
    #latexdict['tabulartype'] = 'longtable'
    latexdict['units'] = {}

    tbl.write('../paper1/linetable{0}.tex'.format(spwid),
              latexdict=latexdict, format='ascii.latex', overwrite=True)

methanol_table = Table(list(map(list, zip(*methanol_rows))))
methanol_table.rename_column('col0', 'Line Name')
methanol_table.rename_column('col1', 'Frequency')
methanol_table['Frequency'].unit = u.GHz
methanol_table.rename_column('col2', 'E$_U$')
methanol_table['E$_U$'].unit = u.K
methanol_table.sort('E$_U$')

latexdict['col_align'] = 'l'*len(methanol_table.columns)
latexdict['header_start'] = '\label{tab:methanol}'
latexdict['caption'] = '\methanol lines used to determine temperature'
methanol_table.write('../paper1/methanol_table.tex', latexdict=latexdict,
                     format='ascii.latex', overwrite=True)

# try to merge tables into one...
with open('../paper1/linetable.tex', 'w') as outfh:

    outfh.write('\\begin{table*}[htp]\n')
    #outfh.write('\caption{Spectral Lines}\n')
    #outfh.write('\label{tab:lines}\n')

    for ii in range(4):
        with open('../paper1/linetable{0}.tex'.format(ii), 'r') as fh:
            lines = fh.readlines()
        outfh.write("\\begin{minipage}[t]{0.5\\textwidth}\n")
        outfh.write("\\centering\n")
        #outfh.writelines(lines)
        outfh.write(lines[1]+"\n") # caption
        outfh.writelines(list(filter(lambda x: True,#'label' not in x,
                                     lines[2:-1])))
        #outfh.write('\\quad\n')
        outfh.write("\\end{minipage}\n")

    outfh.write('\\end{table*}')
