import numpy as np
import pyregion
import paths
from astropy.io import fits
import radio_beam
from astropy import constants
from astropy import units as u
from astropy import wcs
import masscalc

fh = fits.open(paths.dpath('12m/continuum/selfcal_allspw_selfcal_3_mfs_deeper_r2.0.image.pbcor.fits'))

mywcs = wcs.WCS(fh[0].header)

beam = radio_beam.Beam.from_fits_header(fh[0].header)

pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5
ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam


for rpath in ('e5bubble_inclusive.reg', 'e5bubble_dustonly.reg'):
    print(rpath)
    reg = pyregion.open(paths.rpath(rpath))

    mask = reg.get_mask(fh[0])
    cutout = fh[0].data*mask

    limits = [0.005, 0.015]

    total_flux_perbeam = cutout[(cutout > limits[0]) & (cutout < limits[1])].sum()
    total_flux = total_flux_perbeam / ppbeam.value
    print("Total flux: {0}".format(total_flux))

    rad = reg[0].coord_list[-1]*u.deg
    rad_pc = (rad*masscalc.distance).to(u.pc, u.dimensionless_angles())

    print("Total mass(20K): {0}".format(total_flux*masscalc.mass_conversion_factor(TK=20)))
    print("Total mass(50K): {0}".format(total_flux*masscalc.mass_conversion_factor(TK=50)))
    print("Assumed radius: {0}".format(rad_pc))
    print("Total density(50K): {0}".format((total_flux*masscalc.mass_conversion_factor(TK=50)/(4/3*np.pi*rad_pc**3) / (2.8*u.Da)).to(u.cm**-3)))

Qlyc = 1e49 * 1/u.s
alpha_b = 3e-13*u.cm**3*u.s**-1
n_h = 5e5*u.cm**-3
rstrom = (3 * Qlyc /(4*np.pi*alpha_b*n_h**2))**(1/3.)
print("Stromgren radius: {0}".format(rstrom))

# sound speed at 8500K
cii = 7.1*u.km/u.s

def r_of_t(Rs, t, cii=cii, ):
    return Rs * (1+7/4. * cii/Rs * t)**(4/7.)

def rdot_of_t(Rs, t, cii=cii, ):
    return Rs * (cii/Rs) * (1+7/4. * cii/Rs * t)**(-3/7.)

print("Radius at t=10^4 yr: {0}".format(r_of_t(rstrom, 1e4*u.yr,).to(u.pc)))
print("Radius at t=10^5 yr: {0}".format(r_of_t(rstrom, 1e5*u.yr,).to(u.pc)))
print("Velocity at t=10^4 yr: {0}".format(rdot_of_t(rstrom, 1e4*u.yr,).to(u.km/u.s)))
print("Velocity at t=10^5 yr: {0}".format(rdot_of_t(rstrom, 1e5*u.yr,).to(u.km/u.s)))


def tfrag(cs, qh, n):
    return 1.56 *u.Myr * (cs/(0.2*u.km/u.s))**(7/11.) * (qh/(1e49*u.s**-1))**(-1/11.) * (n/(1e3*u.cm**-3))**(-5/11.)

cs_50 = ((constants.k_B * 50*u.K / (2.8*u.Da))**0.5).to(u.km/u.s)
print("Fragmentation time: {0}".format(tfrag(cs_50, Qlyc, n_h*2)))

def rfrag(cs, qh, n):
    return 5.8 * u.pc * (cs/(0.2*u.km/u.s))**(4/11.) * (qh/(1e49*u.s**-1))**(1/11.) * (n/(1e3*u.cm**-3))**(-6/11.)

print("Fragmentation size scale: {0}".format(rfrag(cs_50, Qlyc, n_h*2)))
