import numpy as np
import radio_beam
from astropy import units as u, constants
from astropy.table import Table,Column
from dust_emissivity import blackbody
from HII_model import qlyc_of_tb, Snu
import pylab as pl

# stellar parameters come from Pecaut & Mamajek 2013 via
# http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
# from evla paper/analysis/magnitudes.py:

distance_modulus = 5*np.log10(5400)-5
# 13.661968799114842

# http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
mkabso9v = - 3.20
mko9v = distance_modulus + mkabso9v
# 10.461968799114842

# goldader 1994
ak = 2.6
mko9vext = mko9v + ak
# 13.061968799114842

# for e2e
snu_max = 0.6*u.mJy
beam = radio_beam.Beam(0.34*u.arcsec)
phys_radius = 0.34 * 5400 * u.au
tbmax = snu_max.to(u.K, u.brightness_temperature(beam, 14.5*u.GHz))
qmax = qlyc_of_tb(tbmax, radius=phys_radius).decompose()
print("Upper limit continuum if thin = {0}".format(qmax))

# # plot of... not relevant
# pl.plot(np.logspace(40,50),
#         Snu(Te=8500*u.K, nu=14.5*u.GHz, R=110*u.au, Qlyc=np.logspace(40,50)*u.s**-1,
#             beam=phys_radius, angular_beam=beam.major))
# pl.plot(np.logspace(40,50), [0.6e-3]*50)



tbl = Table.read('pecaut_mamajek_table.txt', format='ascii')

wav = np.linspace(100*u.AA, 10000*u.AA, 100000)
lycfrac = []
for row in tbl:
    if row['Teff'] > 7000:
        bb=blackbody.blackbody_wavelength(wav, row['Teff']*u.K)
        lycfrac.append(bb[wav<912*u.AA].sum()/bb.sum())
    else:
        lycfrac.append(0)

tbl.add_column(Column(name='lycfrac', data=lycfrac))
lyclum = (lycfrac*(10**np.array(tbl['logL'], dtype='float')*u.L_sun)/(constants.h*constants.c/(912*u.AA))).to(u.s**-1)
tbl.add_column(Column(name='lyclum', data=lyclum))

#print(tbl['SpT','Msun','lyclum','logL','Teff'][:40])

closest_type = np.argmin(np.abs(tbl['lyclum'] - qmax))

print("Closest spectral type = {0}".format(tbl['SpT'][closest_type]))
print(tbl['SpT','Msun','lyclum','logL','Teff'][closest_type-1:closest_type+1])
