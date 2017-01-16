from __future__ import print_function
import os
import shutil
import datetime
from astroquery import lamda
import radmc3dPy
import matplotlib
import pylab as pl
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import constants
from yt.analysis_modules.radmc3d_export.api import RadMC3DWriter, RadMC3DSource
from yt.utilities.physical_constants import kboltz
import yt
import struct
from core_models import broken_powerlaw
from convert_to_K import convert_to_K
from get_dust_opacity import get_dust_opacity

def read_dust_temperature(dust_tem_fn, sz=32):
    with open(dust_tem_fn, 'rb') as fh:
        ftype, = struct.unpack('=q', fh.read(8))
        precis, = struct.unpack('=q', fh.read(8))
        nrcells, = struct.unpack('=q', fh.read(8))
        nrspec, = struct.unpack('=q', fh.read(8))
        data = np.fromfile(fh, dtype='float64', count=nrcells)
    print(dust_tem_fn, ftype, precis, nrcells, nrspec)

    # raise an IOError, because the wrong dust temperature file has been found
    if not sz * sz * sz == nrcells:
        raise IOError("Wrong file dimensions")

    return data.reshape([sz, sz, sz])


def core_model_dust(outname, x_co=1.0e-4, x_h2co=1.0e-9, x_ch3oh=1e-9, zh2=2.8,
                    sz=16, max_rad=10000*u.au, rbreak=1000*u.au,
                    radius_cm=1*u.au.to(u.cm), mass_g=1*u.M_sun.to(u.g),
                    n0=5e8*u.cm**-3,
                    power=-1.5,
                    recompute_dusttemperature=True,
                    luminosity=2e4*u.L_sun,):

    if not os.path.exists('dustkappa_mrn5.inp'):
        get_dust_opacity()


    mu_h2 = yt.YTArray(zh2 * u.Da.to(u.g), 'g')


    # Problem setup: pure density field
    zz,yy,xx = np.indices([sz,sz,sz])
    rr = ((zz-(sz-1)/2.)**2 + (yy-(sz-1)/2.)**2 + (xx-(sz-1)/2.)**2)**0.5
    max_velo = 1*u.km/u.s
    velo = max_velo - np.array([(sz-1)/2.-zz, (sz-1/2.)-yy, (sz-1/2.)-xx]) / rr.max() * max_velo

    # now rr has units
    rr = rr * max_rad / (sz/2.)
    dens = broken_powerlaw(rr, rbreak=rbreak, n0=n0, power=power)

    data = {('gas','density'): ((dens*u.Da*zh2).to(u.g/u.cm**3), "g/cm**3"),
            ('gas','z_velocity'): (velo[0].to(u.km/u.s).value, 'km/s'),
            ('gas','y_velocity'): (velo[1].to(u.km/u.s).value, 'km/s'),
            ('gas','x_velocity'): (velo[2].to(u.km/u.s).value, 'km/s'),
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


    get_dust_opacity()

    params=dict(istar_sphere=0, itempdecoup=0, lines_mode=3, nphot=10000,
                nphot_scat=3000, nphot_spec=10000, rto_style=3,
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

    if recompute_dusttemperature:
        # compute the dust temperature
        assert os.system('radmc3d mctherm') == 0

        dust_temperature = read_dust_temperature('dust_temperature.bdat', sz=sz)
        shutil.move('dust_temperature.bdat','dust_temperature_{0}.bdat'.format(outname))
    else:
        try:
            shutil.copy('dust_temperature_{0}.bdat'.format(outname),
                        'dust_temperature.bdat',)
            success = True
        except Exception as ex:
            success = False
            print(ex)
        dust_temperature = read_dust_temperature('dust_temperature.bdat', sz=sz)
        if success:
            os.remove('dust_temperature.bdat')

    if os.path.exists('lines.inp'):
        os.remove('lines.inp')
    assert os.system('radmc3d image npix 50 incl 0 sizeau 10000 noscat pointau 0.0  0.0  0.0 fluxcons lambda 1323 dpc 5400') == 0
    im = radmc3dPy.image.readImage('image.out')
    im.writeFits('dustim1323um_{0}.fits'.format(outname), fitsheadkeys={}, dpc=5400,
                 coord='19h23m43.963s +14d30m34.56s', overwrite=True)

    return dust_temperature

if __name__ == "__main__":
    tmplt = 'dust_temperature_sz32_rad1e4au_mstar1msun_rstar1au_lstar{0:0.1e}lsun_power{1}_ncen{2:0.1e}.bdat'

    max_rad = 10000*u.au

    sz = 32
    zz,yy,xx = np.indices([sz,sz,sz])
    rr = ((zz-(sz-1)/2.)**2 + (yy-(sz-1)/2.)**2 + (xx-(sz-1)/2.)**2)**0.5
    rr = rr * max_rad / (sz/2.)
    rr_u, inds = np.unique(rr.ravel(), return_index=True)
    rr_2d = ((yy[0,:,:]-(sz-1)/2.)**2 + (xx[0,:,:]-(sz-1)/2.)**2)**0.5
    rr_2d = rr_2d * max_rad / (sz/2.)
    rr_2du, inds_2d = np.unique(rr_2d.ravel(), return_index=True)


    fig2 = pl.figure(2)
    fig2.clf()

    fig3 = pl.figure(3)
    fig3.clf()

    linestyles = {-2.0: '--',
                  -1.5: '-',
                  -1.0: ':'}
    colors = {1e4: 'r',
              2e4: 'g',
              5e4: 'b',
              1e5: 'k'}

    # set up the appropriately-sized grid
    lstar = 2e4
    power = -1.5
    ncen = 5e7
    outname = "sz{2}_rad1e4au_mstar1msun_rstar1au_lstar{0:0.1e}lsun_power{1}_ncen{3:0.1e}".format(lstar,power,sz,ncen)
    try:
        core_model_dust(outname=outname,
                        x_co=1.0e-4, x_h2co=1.0e-9, x_ch3oh=1e-9, zh2=2.8, sz=sz,
                        max_rad=max_rad, rbreak=1000*u.au,
                        recompute_dusttemperature=False,
                        radius_cm=1*u.au.to(u.cm), mass_g=1*u.M_sun.to(u.g),
                       )
    except IOError:
        core_model_dust(outname=outname,
                        x_co=1.0e-4, x_h2co=1.0e-9, x_ch3oh=1e-9, zh2=2.8, sz=sz,
                        max_rad=max_rad, rbreak=1000*u.au,
                        recompute_dusttemperature=True,
                        radius_cm=1*u.au.to(u.cm), mass_g=1*u.M_sun.to(u.g),
                       )

    for central_density in (5e7,5e8,):

        for power in (-1.5, -2.0, -1.0):

            for lstar in (2e4, 1e4, 5e4, 1e5):

                if os.path.exists('dust_temperature.bdat'):
                    shutil.move('dust_temperature.bdat', 'dust_temperature_backup_{0}.bdat'.format(datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')))

                dusttem_fn = tmplt.format(lstar, power, central_density)
                outname = "sz{2}_rad1e4au_mstar1msun_rstar1au_lstar{0:0.1e}lsun_power{1}_ncen{3:0.1e}".format(lstar,power,sz,central_density)
                if not os.path.exists(dusttem_fn):
                    core_model_dust(outname=outname,
                                    x_co=1.0e-4, x_h2co=1.0e-9, x_ch3oh=1e-9, zh2=2.8, sz=sz,
                                    max_rad=max_rad, rbreak=1000*u.au,
                                    n0=central_density*u.cm**-3,
                                    radius_cm=1*u.au.to(u.cm), mass_g=1*u.M_sun.to(u.g),
                                    power=power, luminosity=lstar*u.L_sun,)
                else:
                    dust_temperature = read_dust_temperature(dusttem_fn, sz)

                
                dust_image = 'dustim1323um_{0}.fits'.format(outname)
                dust_image_K = 'dustim1323um_{0}_K.fits'.format(outname)
                if not os.path.exists(dust_image_K):
                    shutil.copy(dusttem_fn, 'dust_temperature.bdat')
                    if os.path.exists('lines.inp'):
                        os.remove('lines.inp')
                    os.system('radmc3d image npix 50 incl 0 sizeau 10000 noscat pointau 0.0  0.0  0.0 fluxcons lambda 1323 dpc 5400')
                    im = radmc3dPy.image.readImage('image.out')
                    im.writeFits(dust_image, fitsheadkeys={}, dpc=5400,
                                 coord='19h23m43.963s +14d30m34.56s', overwrite=True)
                    im_K = convert_to_K(dust_image)
                    im_K.writeto(dust_image_K, clobber=True)
                imdata = fits.getdata(dust_image_K)

                fig1 = pl.figure(1)
                fig1.clf()
                ax1 = pl.subplot(1,2,1)
                im = ax1.imshow(dust_temperature[:,:,sz/2], cmap='hot')
                pl.colorbar(im, ax=ax1)
                ax2 = pl.subplot(1,2,2)
                ax2.plot(rr.ravel(), dust_temperature.ravel(), '.', alpha=0.25)
                pl.savefig("midplane_dust_temperature_{0}.png".format(outname))


                ax3 = fig2.gca()
                ax3.plot(rr_u, dust_temperature.flat[inds], linestyle=linestyles[power],
                         color=colors[lstar], alpha=1.0, label="$L={0:0.1e}, \\kappa={1}$".format(lstar, power))

                ax4 = fig3.gca()
                ax4.plot(rr_2du, dust_temperature.flat[inds_2d], linestyle=linestyles[power],
                         color=colors[lstar], alpha=1.0, label="$L={0:0.1e}, \\kappa={1}$".format(lstar, power))


    # overlay "core_model" version
    sz = 16
    lstar = 2e4
    power = -1.5
    zz,yy,xx = np.indices([sz,sz,sz])
    rr = ((zz-(sz-1)/2.)**2 + (yy-(sz-1)/2.)**2 + (xx-(sz-1)/2.)**2)**0.5
    rr = rr * max_rad / (sz/2.)
    rr_u, inds = np.unique(rr.ravel(), return_index=True)
    rr_2d = ((yy[0,:,:]-(sz-1)/2.)**2 + (xx[0,:,:]-(sz-1)/2.)**2)**0.5
    rr_2d = rr_2d * max_rad / (sz/2.)
    rr_2du, inds_2d = np.unique(rr_2d.ravel(), return_index=True)

    # set up the appropriately-sized grid
    outname = "sz{2}_rad1e4au_mstar1msun_rstar1au_lstar{0:0.1e}lsun_power{1}_ncen5e8".format(lstar,power,sz)
    core_model_dust(outname=outname,
                    x_co=1.0e-4, x_h2co=1.0e-9, x_ch3oh=1e-9, zh2=2.8, sz=sz,
                    n0=5e8*u.cm**-3,
                    max_rad=max_rad, rbreak=1000*u.au,
                    recompute_dusttemperature=False,
                    radius_cm=1*u.au.to(u.cm), mass_g=1*u.M_sun.to(u.g),
                   )

    dusttem_fn = outname+'.bdat'
    dust_temperature = read_dust_temperature(dusttem_fn, sz)
    ax3 = fig2.gca()
    ax3.plot(rr_u, dust_temperature.flat[inds], linestyle='-',
             color='m', alpha=0.5, linewidth=2, label="$L={0:0.1e}, \\kappa={1}$, $n=5\\times10^8$".format(lstar, power))

    dust_image = 'dustim1323um_{0}.fits'.format(outname)
    imdata = fits.getdata(dust_image)
    ax4 = fig3.gca()
    ax4.plot(rr_2du, dust_temperature.flat[inds_2d], linestyle='-',
             color='m', alpha=0.5, label="$L={0:0.1e}, \\kappa={1}$, $n=5\\times10^8$".format(lstar, power))

    pl.figure(2)
    pl.legend(loc='best')
    pl.xlabel("Radius (AU)", fontsize=16)
    pl.ylabel("Dust temperature (K)", fontsize=16)
    fig2.savefig("dust_temperature_radial_profile_comparison.png")

    pl.figure(3)
    pl.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
    pl.xlabel("Radius (AU)", fontsize=16)
    pl.ylabel("Dust brightness (Jy/beam?)", fontsize=16)
    fig3.savefig("dust_brightness_radial_profile_comparison.png", bbox_inches='tight')
