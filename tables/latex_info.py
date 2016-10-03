import numpy as np
from astropy import units as u
from astropy.io import ascii

def exp_to_tex(st):
    if st == 'nan':
        return '-'
    elif 'e' in st:
        pt1,pt2 = st.split('e')
        return "{0}\\ee{{{1:d}}}".format(pt1,int(pt2))
    return st

def format_float(st):
    return exp_to_tex("{0:0.2g}".format(st))


colnames = ['Object Name',
            'Amplitude',
            '$\sigma(Amplitude)$',
            '$V_{LSR}$',
            '$\sigma(V_{LSR})$',
            '$dV$',
            '$\sigma(dV)',
            '$\Omega_{ap}$',]
units = {'Amplitude':u.mJy.to_string(u.format.LatexInline),
         '$\sigma(Amplitude)$':u.mJy.to_string(u.format.LatexInline),
         '$V_{LSR}$':(u.km/u.s).to_string(u.format.LatexInline),
         '$\sigma(V_{LSR})$':(u.km/u.s).to_string(u.format.LatexInline),
         '$dV$':(u.km/u.s).to_string(u.format.LatexInline),
         '$\sigma(dV)':(u.km/u.s).to_string(u.format.LatexInline),
         '$\Omega_{ap}$':u.sr.to_string(u.format.LatexInline),
         '$r_{eff}$':u.arcsec.to_string(u.format.LatexInline),
        }
latexdict = ascii.latex.latexdicts['AA']
latexdict['tabletype'] = 'table*'
latexdict['tablealign'] = 'htp'
latexdict['units'] = units

def rounded(value, error, extra=1):
    """
    Return the value and error both rounded to the error's first digit
    """

    if error == 0:
        return (0,0)

    if hasattr(value, 'unit'):
        value = value.value

    digit = int(np.ceil(-np.log10(error))) + extra
    assert np.round(error, digit) != 0
    return np.round(value, digit), np.round(error, digit)#, digit

def round_to_n(x, n):
    if np.isnan(x):
        return np.nan
    else:
        return round(x, -int(np.floor(np.log10(x))) + (n - 1))

def strip_trailing_zeros(x):
    return x.rstrip('.0')
