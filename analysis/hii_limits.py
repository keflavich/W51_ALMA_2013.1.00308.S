import numpy as np
import radio_beam
from astropy import units as u, constants
from astropy.table import Table,Column
from dust_emissivity import blackbody
from HII_model import qlyc_of_tb, Snu, EM_of_T
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

te = 8500 * u.K

# for e2e
beam = radio_beam.Beam(0.34*u.arcsec)

freq = 45 * u.GHz
snu_max = 1.4*u.mJy # measured detection (includes dust)
phys_radius = 0.05 * 5400 * u.au
tbmax = snu_max.to(u.K, u.brightness_temperature(beam, 14.5*u.GHz))
qmax = qlyc_of_tb(tbmax, radius=phys_radius, nu=freq, Te=te).decompose()
print("Upper limit continuum if thin @ 45 GHz = {0}".format(qmax))
EM = EM_of_T(tbmax, nu=freq, Te=te)
dens = ((EM/phys_radius)**0.5).to(u.cm**-3)
print("Corresponding assumed density: {0}".format(dens))
thick_size_limit = (tbmax/te)**0.5 * phys_radius
print("Optically Thick Size Upper Limit: {0}".format(thick_size_limit))
print()

freq = 4.8 * u.GHz
snu_max = 0.4*u.mJy # upper limit, eyeballed
phys_radius = 0.34 * 5400 * u.au # not exact
tbmax = snu_max.to(u.K, u.brightness_temperature(beam, 14.5*u.GHz))
qmax = qlyc_of_tb(tbmax, radius=phys_radius, nu=freq, Te=te).decompose()
print("Upper limit continuum if thin @ 4.8 GHz = {0}".format(qmax))
EM = EM_of_T(tbmax, nu=freq, Te=te)
dens = ((EM/phys_radius)**0.5).to(u.cm**-3)
print("Corresponding assumed density: {0}".format(dens))
thick_size_limit = (tbmax/te)**0.5 * phys_radius
print("Optically Thick Size Upper Limit: {0}".format(thick_size_limit))
print()
 
snu_max = 0.6*u.mJy # upper limit, 2-sigma
phys_radius = 0.34 * 5400 * u.au
tbmax = snu_max.to(u.K, u.brightness_temperature(beam, 14.5*u.GHz))
freq = 14.5*u.GHz
qmax = qlyc_of_tb(tbmax, radius=phys_radius, nu=freq, Te=te).decompose()
print("Upper limit continuum if thin @ 14.5 GHz = {0}".format(qmax))
EM = EM_of_T(tbmax, nu=freq, Te=te)
dens = ((EM/phys_radius)**0.5).to(u.cm**-3)
print("Corresponding assumed density: {0}".format(dens))
thick_size_limit = (tbmax/te)**0.5 * phys_radius
print("Optically Thick Size Upper Limit: {0}".format(thick_size_limit))
print()



# s_peak(225) ~ 15 mJy
# s_thick(45) = 0.6 mJy
# s_measured(45) = 1.4 mJy

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

colorder = ['SpT', 'Teff', 'logT', 'Msun', 'logL', 'lycfrac',
            'lyclum', 'B-V', 'Bt-Vt', 'U-B', 'V-Rc', 'V-Ic', 'V-Ks', 'J-H',
            'H-K', 'Ks-W1', 'logAge', 'b-y', '#SpT', 'M_J', 'M_Ks', 'Mbol',
            'i-z', 'z-Y', 'W1-W2', 'BCv', 'Mv', ]

tbl[colorder].write('pecaut_table_with_lyclums.txt', format='ascii.fixed_width', overwrite=True)

#print(tbl['SpT','Msun','lyclum','logL','Teff'][:40])

closest_type = np.argmin(np.abs(tbl['lyclum'] - qmax))

print("Closest spectral type = {0}".format(tbl['SpT'][closest_type]))
print(tbl['SpT','Msun','lyclum','logL','Teff'][closest_type-1:closest_type+1])


# IMF examinations
from imf import imf
cluster_masses = [250, 500, 1000, 2000]
nsamp = 100
clusters = {mass: [imf.make_cluster(mass, silent=True) for ii in range(nsamp)]
            for mass in cluster_masses}
luminosities = {mass: [imf.lum_of_cluster(cl) for cl in clusters[mass]]
                for mass in cluster_masses}
lycluminosities = {mass: [imf.lyc_of_cluster(cl) for cl in clusters[mass]]
                   for mass in cluster_masses}

for mass in cluster_masses:
    pl.plot(luminosities[mass], lycluminosities[mass], '.', label=mass)
    
pl.plot(np.log10(2e4), np.log10(4e45), 'x')

pl.legend(loc='best')
