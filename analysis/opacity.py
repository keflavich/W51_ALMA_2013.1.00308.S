import numpy as np
import pylab as pl
import dust_emissivity
from astropy import constants
from astropy import units as u

temperatures = np.linspace(10, 2000, 100)
kappas = []
for temperature in temperatures:
    kappas.append(dust_emissivity.dust.planck_average_kappa_table(temperature*u.K).value)

kappas = np.array(kappas) * u.cm**2/u.g

pl.figure(1).clf()
pl.plot(temperatures, kappas)
pl.xlabel("Temperature")
pl.ylabel("Planck-averaged opacity (cm$^2$/g)")


# Eddington luminosity...

Ledd = (constants.G*20*u.M_sun*constants.c/(4*np.pi*kappas)).to(u.L_sun)
pl.figure(2).clf()
pl.loglog(temperatures, Ledd)
pl.xlabel("Temperature")
pl.ylabel("Dust-only Eddington Luminosity (L$_\\odot$)")
