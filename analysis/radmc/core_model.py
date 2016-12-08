import numpy as np
from astropy import units as u
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter, RadMC3DSource
from yt.utilities.physical_constants import kboltz
import yt
from core_models import broken_powerlaw

x_co = 1.0e-4
x_h2co = 1.0e-9
x_ch3oh = 1e-9
zh2 = 2.8
mu_h2 = yt.YTArray(zh2 * u.Da.to(u.g), 'g')


# Problem setup: pure density field
sz = 32
max_rad = 10000*u.au
rbreak = 1000*u.au
zz,yy,xx = np.indices([sz,sz,sz])
rr = ((zz-sz/2)**2 + (yy-sz/2)**2 + (xx-sz/2)**2)**0.5
rr = rr * max_rad / (sz/2.)
dens = broken_powerlaw(rr, rbreak=rbreak, n0=1e8*u.cm**-3, power=-1.5)


data = dict(density=((dens*u.Da*zh2).to(u.g/u.cm**3), "g/cm**3"))
bbox = np.array([[-max_rad.value,max_rad.value]]*3)
ds = yt.load_uniform_grid(data, dens.shape, length_unit="au", bbox=bbox, nprocs=64)

dust_to_gas = 0.01
def _DustDensity(field, data):
    return dust_to_gas * data['density']
ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3")

def _NumberDensityCH3OH(field, data):
    return (1./mu_h2)*x_h2co*data['density']*x_ch3oh # data['density']#
ds.add_field(("gas", "number_density_CH3OH"), function=_NumberDensityCH3OH, units="cm**-3")

def _NumberDensityH2(field, data):
    return (1./mu_h2)*data['density'] # data['density']#
ds.add_field(("gas", "number_density_H2"), function=_NumberDensityH2, units="cm**-3")

def _GasTemperature(field, data):
    return yt.YTArray(np.ones_like(data['density'])*100., 'K')
ds.add_field(("gas", "temperature"), function=_GasTemperature, units='K')

def _DustTemperature(field, data):
    return yt.YTArray(np.ones_like(data['density'])*100., 'K')
ds.add_field(("gas", "dust_temperature"), function=_DustTemperature, units='K')

writer = RadMC3DWriter(ds)

writer.write_amr_grid()
#writer.write_line_file(("gas", "number_density_CO"), "numberdens_co.inp")
dens_fn = "numberdens_h2.inp" # not h2_numberdens.inp?
writer.write_line_file(("gas", "number_density_H2"), dens_fn)
writer.write_line_file(("gas", "number_density_CH3OH"), "numberdens_ch3oh.inp")
#writer.write_line_file(("gas", "number_density_CH3OH"), "ch3oh_numberdens.inp")
#writer.write_line_file(("gas", "number_density_H2CO"), "numberdens_h2co.inp")
writer.write_dust_file(("gas", "temperature"), "gas_temperature.inp")
writer.write_dust_file(("gas", "dust_density"), "dust_density.inp")
#writer.write_dust_file(("gas", "dust_temperature"), "dust_temperature.inp")

#velocity_fields = [('deposit',"all_cic_velocity_x"), ('deposit',"all_cic_velocity_y"),
#                   ('deposit',"all_cic_velocity_z")]
#writer.write_line_file(velocity_fields, "gas_velocity.inp")

# central star
radius_cm = 1*u.au.to(u.cm)
mass_g = 1*u.M_sun.to(u.g)
position_cm = [0.0, 0.0, 0.0]
temperature_K = 1000.0
star = RadMC3DSource(radius_cm, mass_g, position_cm, temperature_K)

sources_list = [star]
wavelengths_micron = np.logspace(-1.0, 4.0, 1000)

writer.write_source_files(sources_list, wavelengths_micron)

import os
import shutil

shutil.copy('/Users/adam/repos/radmc-3d/version_0.39/python/python_examples/datafiles/dustkappa_silicate.inp', '.')
shutil.copy('/Users/adam/work/jimsims/code/dustopac.inp',
            'dustopac.inp')
shutil.copy('/Users/adam/repos/radmc-3d/version_0.39/python/python_examples/datafiles/molecule_co.inp', '.')
shutil.copy('/Users/adam/LAMDA/e-ch3oh.dat','molecule_ch3oh.inp')

params=dict(istar_sphere=0, itempdecoup=1, lines_mode=3, nphot=1000000,
            nphot_scat=30000, nphot_spec=100000, rto_style=3,
            scattering_mode=0, scattering_mode_max=1, tgas_eq_tdust=1,)

params_string = """
istar_sphere = {istar_sphere}
itempdecoup = {itempdecoup}
lines_mode = {lines_mode}
nphot = {nphot}
nphot_scat = {nphot_scat}
nphot_spec = {nphot_spec}
rto_style = {rto_style}
scattering_mode = {scattering_mode}
scattering_mode_max = {scattering_mode_max}
tgas_eq_tdust = {tgas_eq_tdust}
"""

with open('wavelength_micron.inp', 'w') as fh:
    fh.write("{0}\n".format(len(wavelengths_micron)))
    for nu in wavelengths_micron:
        fh.write("{0}\n".format(nu))

# with open('stars.inp', 'w') as fh:
#     fh.write("2\n")
#     nstars = 1
#     nlam = nfrq
#     fh.write("{0} {1}\n".format(nstars,nlam))
#     rstar = (1*u.au).to(u.cm).value
#     mstar = 15*u.M_sun.to(u.g)
#     x,y,z = 0,0,0
#     fh.write("{0} {1} {2} {3} {4}\n".format(rstar, mstar, x, y, z))
#     for nu in wavelengths_micron:
#         fh.write("{0}\n".format(nu))
#     temperature = 1000 * u.K
#     for nu in wavelengths_micron:
#         fh.write("{0}\n".format(-temperature.to(u.K).value))


with open('radmc3d.inp','w') as f:
    params['lines_mode'] = 50
    f.write(params_string.format(**params))

assert os.system('radmc3d calcpop writepop') == 0

with open('radmc3d.inp','w') as f:
    params['lines_mode'] = 3 # 3 = sobolev (LVG)
    f.write(params_string.format(**params))

# compute the dust temperature
assert os.system('radmc3d mctherm') == 0

with open('lines.inp','w') as fh:
    fh.write("2\n")
    # n_lines
    fh.write("1\n")
    # line name (filename molecule_{linename}.inp must exist)
    # file type (leiden)
    # 0 # expert-only
    # 0 # expert-only
    # 0 for LTE or 1 for nLTE
    fh.write("ch3oh leiden 0 0 1\n")
    fh.write("h2")

lengthscale = (1e4*u.au.to(u.cm))
with open('escprob_lengthscale.inp','w') as fh:
    with open(dens_fn,'r') as h2fh:
        lines = h2fh.read().split("\n")
        nlines = len(lines)
    fh.write("\n".join(lines[0:3]) + "\n")
    # want to write 1 less than nlines
    for ii in range(nlines-1):
        fh.write("{0}\n".format(lengthscale))

import radmc3dPy
# iline: 1 = CO 1-0, 2 = CO 2-1, etc.
# widthkms = full width of output spectrum, divided by linenlam
# linenlam: number of wavelength bins
# linelist
wavelength_center = (218.440063*u.GHz).to(u.um, u.spectral()).value

assert os.system('radmc3d calcpop writepop noscat nostar') == 0

radmc3dPy.image.makeImage(iline=240, widthkms=1, linenlam=40,
                          nostar=True,
                          noscat=True,
                          vkms=0,
                          npix=50, incl=0,
                          lambdarange=[wavelength_center*(1-10/3e5),
                                       wavelength_center*(1+10/3e5)],
                          nlam=40, sizeau=10000)
shutil.move('image.out', 'ch3oh_422-312_image.out')
im = radmc3dPy.image.readImage('ch3oh_422-312_image.out')
im.writeFits('ch3oh_422-312_image.fits', fitsheadkeys={}, dpc=5400,
             coord='19h23m43.963s +14d30m34.56s')





#radmc3dPy.image.makeImage(iline=240, widthkms=5, linenlam=40, nostar=True,
#                          npix=500, incl=0,
#                          lambdarange=[wavelength_center*(1-10/3e5),
#                                       wavelength_center*(1+10/3e5)],
#                          nlam=20, sizeau=10000)
#shutil.move('image.out', 'ch3oh_422-312_image.out')
#im = radmc3dPy.image.readImage('ch3oh_422-312_image.out')
#im.writeFits('ch3oh_422-312_image.fits', fitsheadkeys={}, dpc=5400, coord='19h23m43.963s +14d30m34.56s')
# radmc3dPy.image.makeImage(iline=3, widthkms=5, linenlam=40, nostar=True)
# shutil.move('image.out', 'h2co_303-202_image.out')
# radmc3dPy.image.makeImage(iline=13, widthkms=5, linenlam=40, nostar=True)
# shutil.move('image.out', 'h2co_321-220_image.out')
# radmc3dPy.image.makeImage(iline=2, widthkms=5, linenlam=40, nostar=True)
# shutil.move('image.out', 'co_2-1_image.out')
# radmc3dPy.image.makeImage(iline=1, widthkms=5, linenlam=40, nostar=True)
# shutil.move('image.out', 'co_1-0_image.out')
# #radmc3d image iline 1 widthkms 10 linenlam 40 linelist nostar
