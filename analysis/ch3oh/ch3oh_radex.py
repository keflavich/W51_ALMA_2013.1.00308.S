import numpy as np
import pyradex
import pylab as pl
from astropy import units as u
from astropy.utils.console import ProgressBar
import paths

# remember that because we've loaded the radex file from fortran, there is only
# one RADEX instance, so we can't play with 2 species simultaneously

R_ch3oh_e = pyradex.Radex(species='e-ch3oh', density=1e5, column=1e15, temperature=50)
te = R_ch3oh_e(density=1e5, column=1e15, temperature=50)
ss_e = R_ch3oh_e.get_synthspec(218*u.GHz, 220.5*u.GHz, npts=4000, linewidth=3*u.km/u.s)

R_ch3oh_a = pyradex.Radex(species='a-ch3oh', density=1e5, column=1e15, temperature=50)
ta = R_ch3oh_a(density=1e5, column=1e15, temperature=50)
ss_a = R_ch3oh_a.get_synthspec(218*u.GHz, 220.5*u.GHz, npts=4000, linewidth=3*u.km/u.s)

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


print(ta[inds_a])
print(te[inds_e])

density = {'H2': 1e7}

tem_taus = {}
columns = np.logspace(12,np.log10(5e16))

for temperature in (50,100,200):
    taus = {x: [] for x in lines_e+lines_a}
    tem_taus[temperature] = taus

    # for speed, split the _e and _a (only one load/interp step each time)
    R_ch3oh_e = pyradex.Radex(species='e-ch3oh', density=density,
                              column=1e15, temperature=temperature)
    for column in ProgressBar(columns):
        te = R_ch3oh_e(density=density, column=column, temperature=temperature)

        for ind, line in zip(inds_e, lines_e):
            taus[line].append(te[ind]['tau'])

    R_ch3oh_a = pyradex.Radex(species='a-ch3oh', density=density,
                              column=1e15, temperature=temperature)
    for column in ProgressBar(columns):
        ta = R_ch3oh_a(density=density, column=column, temperature=temperature)

        for ind, line in zip(inds_a, lines_a):
            taus[line].append(ta[ind]['tau'])

print(ta[inds_a])
print(te[inds_e])

pl.clf()

for tem, color in zip(tem_taus, ('r','k','b')):
    taus = tem_taus[tem]

    for line, style, width, alpha in zip(taus, ("-","--",":","-.","-"),
                                         (1,1,1,1,2), (1,1,1,1,0.5)):
        if len(taus[line]) == len(columns):
            pl.semilogx(columns, taus[line], label="{0}-{1} {2}".format(line[0],
                                                                        line[1],
                                                                        tem),
                        alpha=alpha,
                        linewidth=width,
                        linestyle=style,
                        color=color)
        else:
            print("{0} has wrong length {1}".format(line, len(taus[line])))
pl.ylim(-0.1, 5)
pl.legend(loc='best')
pl.savefig(paths.fpath('ch3oh_radex_tau_vs_column.png'))

