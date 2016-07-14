import pylab as pl
from astropy import units as u
from astroquery.lamda import Lamda, utils as lamutils
from astropy.utils.console import ProgressBar

ech3oh = Lamda.query('e-ch3oh')


lamutils.ncrit(ech3oh, 256, 253, 100)

crates, transitions, levels = ech3oh
upperenergy = [transitions['E_u(K)'][transitions['Upper'] == upper][0]
               for upper in crates['H2']['Upper']]
ncrits = [lamutils.ncrit(ech3oh, upper, lower, 100)
          for upper, lower in ProgressBar(list(zip(transitions['Upper'],
                                                   transitions['Lower'])))]

Jvals = {level: int(levels[levels['Level']==level]['J'][0].split("_")[0])
         for level in levels['Level']}

pl.figure(4).clf()
pl.semilogy(upperenergy, crates['H2']['C_ij(T=200)'], '.')
pl.xlabel("$E_U$ (K)")
pl.ylabel("$C_{ij}$ at $T=200$ K")

pl.figure(1).clf()
pl.semilogy(transitions['Upper'], u.Quantity(ncrits), '.')
pl.xlabel("Level ID")
pl.ylabel("$n_{crit}$ [cm$^{-3}$]")

pl.figure(3).clf()
pl.semilogy([Jvals[lid] for lid in transitions['Upper']], u.Quantity(ncrits), '.')
pl.xlabel("J")
pl.ylabel("$n_{crit}$ [cm$^{-3}$]")

pl.figure(2).clf()
pl.loglog(transitions['EinsteinA'], u.Quantity(ncrits), '.')
pl.xlabel("$A_{ij}$ [s$^{-1}$]")
pl.ylabel("$n_{crit}$ [cm$^{-3}$]")
