import line_to_image_list
from astropy import table
from astroquery.splatalogue import Splatalogue
from astropy import units as u

tbl_list = []
tbl_dict = {}

for linename, freq, _, _ in line_to_image_list.line_to_image_list:
    frq = float(freq.strip("GHz"))*u.GHz

    tbl = Splatalogue.query_lines(frq*(1-0.25/3e5), frq*(1+0.25/3e5), line_lists=['SLAIM'], noHFS=True)
    if len(tbl) > 0:
        col = tbl['Lovas/AST Intensity']
        if col.dtype.kind != 'U':
            col = col.astype(str)
            tbl.replace_column(col.name, col)
        tbl_list.append(tbl)
        tbl_dict[linename] = tbl

    if len(tbl) > 1:
        print("Contaminants for {0}: ".format(linename))
        print(tbl)

final_tbl = table.vstack(tbl_list)

final_tbl.write("full_line_table.csv")
