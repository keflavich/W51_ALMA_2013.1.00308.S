import numpy as np
import pyradex
import pylab as pl
from astropy import units as u
from astropy.utils.console import ProgressBar

# remember that because we've loaded the radex file from fortran, there is only
# one RADEX instance, so we can't play with 2 species simultaneously

R_ch3oh_e = pyradex.Radex(species='e-ch3oh', density=1e5, column=1e15, temperature=50)
te = R_ch3oh_e(density=1e5, column=1e15, temperature=50)
ss_e = R_ch3oh_e.get_synthspec(218*u.GHz, 220.5*u.GHz, npts=4000, linewidth=3*u.km/u.s)

R_ch3oh_a = pyradex.Radex(species='a-ch3oh', density=1e5, column=1e15, temperature=50)
ta = R_ch3oh_a(density=1e5, column=1e15, temperature=50)
ss_a = R_ch3oh_a.get_synthspec(218*u.GHz, 220.5*u.GHz, npts=4000, linewidth=3*u.km/u.s)

print(ta[(ta['frequency'] > 218) & (ta['frequency'] < 220.5)])
print(te[(te['frequency'] > 218) & (te['frequency'] < 220.5)])

frq, prof = ss_e.wcs, ss_e.get_profile()+ss_a.get_profile()

pl.plot(frq, prof)

lines_e = [('4_2','3_1'),
           ('8_0','7_1'),
           ('5_-4','6_-3'),
          ]

lines_a = [('4_-2_0','5_-1_0'),
           ('10_-2_','9_-3_0')
          ]

inds_e = [np.where((te['upperlevel'] == lev[0].encode()) &
                   (te['lowerlevel'] == lev[1].encode()))[0][0]
          for lev in lines_e]

inds_a = [np.where((ta['upperlevel'] == lev[0].encode()) &
                   (ta['lowerlevel'] == lev[1].encode()))[0][0]
          for lev in lines_a]

density = 1e6

tem_taus = {}

for temperature in (50,100,200):
    taus = {x: [] for x in lines_e+lines_a}
    tem_taus[temperature] = taus

    # for speed, split the _e and _a (only one load/interp step each time)
    R_ch3oh_e = pyradex.Radex(species='e-ch3oh', density=density,
                              column=column, temperature=temperature)
    for column in ProgressBar(np.logspace(12,18)):
        te = R_ch3oh_e(density=density, column=column, temperature=temperature)

        for ind, line in zip(inds_e, lines_e)
            taus[line].append(te[ind]['tau'])

    R_ch3oh_a = pyradex.Radex(species='e-ch3oh', density=density,
                              column=column, tamperature=tamperature)
    for column in ProgressBar(np.logspace(12,18)):
        ta = R_ch3oh_a(density=density, column=column, tamperature=tamperature)

        for ind, line in zip(inds_a, lines_a)
            taus[line].append(ta[ind]['tau'])
