import numpy as np
from astropy import units as u
from constants import distance

def mass_conversion_factor(TK=20, d=distance.to(u.kpc).value):
    return 14.30 * (np.exp(13.01/TK) - 1)*d**2

def col_conversion_factor(TK=20):
    return 2.19e22 * (np.exp(13.01/TK - 1))

def co21_conversion_factor():
    # eqn 6 of Ginsburg 2009
    return 3.27e18 * u.cm**-2 / (u.K*u.km/u.s)
