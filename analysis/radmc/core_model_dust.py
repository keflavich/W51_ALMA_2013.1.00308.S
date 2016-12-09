from __future__ import print_function
import os
import shutil
from astroquery import lamda
import radmc3dPy
import matplotlib
import pylab as pl
import numpy as np
from astropy import units as u
from astropy import constants
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter, RadMC3DSource
from yt.utilities.physical_constants import kboltz
import yt
import struct
from core_models import broken_powerlaw

def core_model_dust(outname, x_co=1.0e-4, x_h2co=1.0e-9, x_ch3oh=1e-9, zh2=2.8,
                    sz=16, max_rad=10000*u.au, rbreak=1000*u.au,
                    radius_cm=1*u.au.to(u.cm), mass_g=1*u.M_sun.to(u.g),
                    power=-1.5,
                    luminosity=2e4*u.L_sun,):



    mu_h2 = yt.YTArray(zh2 * u.Da.to(u.g), 'g')


    # Problem setup: pure density field
    zz,yy,xx = np.indices([sz,sz,sz])
    rr = ((zz-(sz-1)/2.)**2 + (yy-(sz-1)/2.)**2 + (xx-(sz-1)/2.)**2)**0.5
    max_velo = 1*u.km/u.s
    velo = max_velo - np.array([(sz-1)/2.-zz, (sz-1/2.)-yy, (sz-1/2.)-xx]) / rr.max() * max_velo

    # now rr has units
    rr = rr * max_rad / (sz/2.)
    dens = broken_powerlaw(rr, rbreak=rbreak, n0=1e8*u.cm**-3, power=power)

    data = {'density': ((dens*u.Da*zh2).to(u.g/u.cm**3), "g/cm**3"),
            'z_velocity': (velo[0].to(u.km/u.s).value, 'km/s'),
            'y_velocity': (velo[1].to(u.km/u.s).value, 'km/s'),
            'x_velocity': (velo[2].to(u.km/u.s).value, 'km/s'),
           }
    bbox = np.array([[-max_rad.value,max_rad.value]]*3)
    ds = yt.load_uniform_grid(data, dens.shape, length_unit="au", bbox=bbox, nprocs=64)

    dust_to_gas = 0.01
    def _DustDensity(field, data):
        return dust_to_gas * data['density']
    ds.add_field(("gas", "dust_density"), function=_DustDensity, units="g/cm**3")

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
    dens_fn = "numberdens_h2.inp" # not h2_numberdens.inp?
    writer.write_line_file(("gas", "number_density_H2"), dens_fn)
    writer.write_dust_file(("gas", "temperature"), "gas_temperature.inp")
    writer.write_dust_file(("gas", "dust_density"), "dust_density.inp")
    #writer.write_dust_file(("gas", "dust_temperature"), "dust_temperature.inp")

    writer.write_line_file([('gas','x_velocity'), ('gas','y_velocity'),
                            ('gas','z_velocity')], "gas_velocity.inp")

    # central star
    position_cm = [0.0, 0.0, 0.0]
    temperature_K = 1000.0
    temperature_K = ((luminosity /
                      (4 * np.pi * (radius_cm*u.cm)**2 * constants.sigma_sb))**0.25
                    ).to(u.K).value
    star = RadMC3DSource(radius_cm, mass_g, position_cm, temperature_K)

    sources_list = [star]
    wavelengths_micron = np.logspace(-1.0, 4.0, 1000)

    writer.write_source_files(sources_list, wavelengths_micron)


    shutil.copy('/Users/adam/repos/radmc-3d/version_0.39/python/python_examples/datafiles/dustkappa_silicate.inp', '.')
    shutil.copy('/Users/adam/work/jimsims/code/dustopac.inp',
                'dustopac.inp')

    params=dict(istar_sphere=0, itempdecoup=0, lines_mode=3, nphot=1000000,
                nphot_scat=30000, nphot_spec=100000, rto_style=3,
                scattering_mode=0, scattering_mode_max=0, tgas_eq_tdust=1,)

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
        params['lines_mode'] = 1 # 3 = sobolev (LVG)
        f.write(params_string.format(**params))

    # compute the dust temperature
    assert os.system('radmc3d mctherm') == 0

    def read_dust_temperature(dust_tem_fn):
        with open(dust_tem_fn, 'rb') as fh:
            ftype, = struct.unpack('=q', fh.read(8))
            precis, = struct.unpack('=q', fh.read(8))
            nrcells, = struct.unpack('=q', fh.read(8))
            nrspec, = struct.unpack('=q', fh.read(8))
            data = np.fromfile(fh, dtype='float64', count=nrcells)
        print(ftype, precis, nrcells, nrspec)

        assert sz * sz * sz == nrcells

        return data.reshape([sz, sz, sz])

    dust_temperature = read_dust_temperature('dust_temperature.bdat')
    shutil.copy('dust_temperature.bdat','dust_temperature_{0}.bdat'.format(outname))

    fig1 = pl.figure(1)
    fig1.clf()
    ax1 = pl.subplot(1,2,1)
    im = ax1.imshow(dust_temperature[:,:,sz/2], cmap='hot')
    pl.colorbar(im, ax=ax1)
    ax2 = pl.subplot(1,2,2)
    ax2.plot(rr.ravel(), dust_temperature.ravel(), '.', alpha=0.25)
    pl.savefig("midplane_dust_temperature_{0}.png".format(outname))

if __name__ == "__main__":

    for power in (-1.5, -2.0, -1.0):
        for lstar in (2e4, 1e4, 5e4, 1e5):
            core_model_dust(outname="sz32_rad1e4au_mstar1msun_rstar1au_lstar{0:0.1e}lsun_power{1}".format(lstar,power),
                            x_co=1.0e-4, x_h2co=1.0e-9, x_ch3oh=1e-9, zh2=2.8, sz=32,
                            max_rad=10000*u.au, rbreak=1000*u.au,
                            radius_cm=1*u.au.to(u.cm), mass_g=1*u.M_sun.to(u.g),
                            power=power, luminosity=lstar*u.L_sun,)
