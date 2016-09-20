import numpy as np
from astropy import units as u
from constants import distance
from astropy import constants
from dust_emissivity import dust

"""
226.6 GHz is the weighted average of these:

218604.028 * 0.5
220296.833 * 1
230435.532 * 1
233040.032 * 1

because 218 is narrower
"""
centerfreq = 226.6*u.GHz
centerfreq_lb = 225*u.GHz # TODO: FIX THIS

#def mass_conversion_factor(TK=20, d=distance.to(u.kpc).value):
#    return 14.30 * (np.exp(13.01/TK) - 1)*d**2
def mass_conversion_factor(TK=20, d=distance.to(u.kpc)):
    return dust.massofsnu(nu=centerfreq, snu=1*u.Jy, distance=d,
                          temperature=u.Quantity(TK, u.K))

#def col_conversion_factor(TK=20):
#    return 2.19e22 * (np.exp(13.01/TK - 1))
def col_conversion_factor(beamomega, TK=20):
    return dust.colofsnu(nu=centerfreq, snu=1*u.Jy, beamomega=beamomega,
                         temperature=u.Quantity(TK, u.K))

def Jnu(T, nu):
    return (2*constants.h*nu**3 / constants.c**2 *
            np.exp(constants.h*nu/(constants.k_B*T)-1)**-1)

def co21_conversion_factor(Tex, Tbg=2.73*u.K):
    # eqn 6 of Ginsburg 2009
    #return 3.27e18 * u.cm**-2 / (u.K*u.km/u.s)

    mu = 1.098e-19*u.esu*u.cm
    Ju = 2 # for co 2-1
    B = 57.64*u.GHz
    nu = 230.538*u.GHz
    Eu = 16.59608*u.K * constants.k_B

    t1 = 3*constants.h/(8*np.pi**3*mu**2*Ju)
    t2 = (constants.k_B * Tex /(constants.h * B) + 1/3.)
    t3 = np.exp(Eu/(constants.k_B*Tex))
    t4 = (np.exp(constants.h*nu/(constants.k_B*Tex))-1)**-1
    # t5 = integral(tau)
    # use Mangum eqn 78
    t5 = (Jnu(Tex, nu) - Jnu(Tbg, nu))**-1

    outunit = 1/(u.K*u.km/u.s)
    result = (t1*t2*t3*t4*t5*constants.k_B).decompose()
    return result.to(outunit)*u.cm**-2

co_abund = 1e-4
