from __future__ import print_function
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

# lvg results in inverted populations for methanol and is therefore useless
lvg = False

do_methanol = True
do_mc_therm = False

x_co = 1.0e-4
x_h2co = 1.0e-9
x_ch3oh = 5e-7
zh2 = 2.8
mu_h2 = yt.YTArray(zh2 * u.Da.to(u.g), 'g')


# Problem setup: pure density field
sz = 16
max_rad = 10000*u.au
rbreak = 1000*u.au
zz,yy,xx = np.indices([sz,sz,sz]) * max_rad / (sz/2.)
mid = max_rad * (sz-1.)/sz
rr = ((zz-mid)**2 + (yy-mid)**2 + (xx-mid)**2)**0.5
max_velo = 2*u.km/u.s
velo = max_velo - u.Quantity([mid-zz, mid-yy, mid-xx]) / rr.max() * max_velo

# now rr has units
dens = broken_powerlaw(rr, rbreak=rbreak, n0=5e8*u.cm**-3, power=-1.5)

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

def _NumberDensityCH3OH(field, data):
    return (1./mu_h2)*data['density']*x_ch3oh # data['density']#
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
writer.write_line_file([('gas','x_velocity'), ('gas','y_velocity'),
                        ('gas','z_velocity')], "gas_velocity.inp")

# central star
radius_cm = 1*u.au.to(u.cm)
mass_g = 1*u.M_sun.to(u.g)
position_cm = [0.0, 0.0, 0.0]
temperature_K = 1000.0
luminosity = 2e4*u.L_sun
temperature_K = ((luminosity /
                  (4 * np.pi * (radius_cm*u.cm)**2 * constants.sigma_sb))**0.25
                ).to(u.K).value
star = RadMC3DSource(radius_cm, mass_g, position_cm, temperature_K)

sources_list = [star]
wavelengths_micron = np.logspace(-1.0, 4.0, 1000)

writer.write_source_files(sources_list, wavelengths_micron)

import os
import shutil

#shutil.copy('/Users/adam/repos/radmc-3d/version_0.39/python/python_examples/datafiles/dustkappa_silicate.inp', '.')
import requests
rslt = requests.get('https://hera.ph1.uni-koeln.de/~ossk/Jena/tables/mrn5')
with open('dustkappa_mrn5.inp','w') as fh:
    lines = rslt.content.rstrip().split("\n")
    wav,opac = np.array([list(map(float, line.split())) for line in lines]).T
    beta = np.log(opac[-1]/opac[-2])/np.log(wav[-2]/wav[-1])
    beta = 1.5
    lastwav = 4000.
    const = opac[-1] * wav[-1]**beta
    lastopac = const * lastwav**-beta
    fh.write("{0:10d}\n".format(1))
    fh.write("{0:10d}\n".format(len(lines)+1))
    for line in lines:
        fh.write("{0}\n".format(line))
    fh.write(" {0:9.3e} {1:9.3e}\n".format(lastwav,lastopac))

dust_type = 'mrn5'

with open('dustopac.inp', 'w') as fh:
    fh.write("2               Format number of this file\n")
    fh.write("1               Nr of dust species\n")
    fh.write("============================================================================\n")
    fh.write("1               Way in which this dust species is read\n")
    fh.write("0               0=Thermal grain, 1=Quantum heated\n")
    fh.write("{0}      Extension of name of dustkappa_***.inp file\n".format(dust_type))
    fh.write("----------------------------------------------------------------------------\n")

wav,opac = np.loadtxt('dustkappa_{0}.inp'.format(dust_type), skiprows=2).T
opac_1300 = np.interp(1300, wav, opac)
print("Max column: {0}".format(((dens[:,8,8] * (max_rad/(sz/2.))).sum().to(u.cm**-2))))
print("Max opacity: {0}".format(((dens[:,8,8] * (max_rad/(sz/2.))).sum().to(u.cm**-2) * opac_1300/100 * u.cm**2/u.g * (zh2*u.Da)).decompose()))


shutil.copy('/Users/adam/repos/radmc-3d/version_0.39/python/python_examples/datafiles/molecule_co.inp', '.')
shutil.copy('/Users/adam/LAMDA/e-ch3oh.dat','molecule_ch3oh.inp')

params=dict(istar_sphere=0, itempdecoup=0, lines_mode=3 if lvg else 1, nphot=1000000,
            nphot_scat=30000, nphot_spec=100000, rto_style=3,
            doppcatch=True,
            scattering_mode=0, scattering_mode_max=0,
            tgas_eq_tdust=0, # usually want this true...
           )

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



# Self consistency check
raddata = radmc3dPy.analyze.readData(ddens=True, dtemp=True, gdens=True,
                                     gtemp=True, ispec='h2', binary=False,
                                     molmass=zh2)
raddata.readDustTemp(binary=True)
print("Max log gas density: {0}".format(np.log10(raddata.ndens_mol.max())))
print("Max log gas mass density: {0}".format(np.log10(raddata.rhogas.max())))
print("Max log dust mass density: {0}".format(np.log10(raddata.rhodust.max())))


os.system('radmc3d image npix 50 incl 0 sizeau 10000 noscat  pointau 0.0  0.0  0.0 fluxcons lambda 1325 tracetau')
im = radmc3dPy.image.readImage('image.out')
pl.figure(7).clf()
pl.imshow(im.image[:,:,0])
pl.colorbar()
pl.savefig("optical_depth_1325um.png")



if do_mc_therm:
    with open('radmc3d.inp','w') as f:
        params['lines_mode'] = 3 if lvg else 1 # 3 = sobolev (LVG)
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

fig1 = pl.figure(1)
fig1.clf()
ax1 = pl.subplot(1,2,1)
im = ax1.imshow(dust_temperature[:,:,sz/2], cmap='hot')
pl.colorbar(im, ax=ax1)
ax2 = pl.subplot(1,2,2)
ax2.plot(rr.ravel(), dust_temperature.ravel(), '.', alpha=0.25)
pl.savefig("midplane_dust_temperature.png")



if do_methanol:

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

    # iline: 1 = CO 1-0, 2 = CO 2-1, etc.
    # widthkms = full width of output spectrum, divided by linenlam
    # linenlam: number of wavelength bins
    # linelist
    wavelength_center = (218.440063*u.GHz).to(u.um, u.spectral()).value

    with open('radmc3d.inp','w') as f:
        params['lines_mode'] = 3 if lvg else 1
        f.write(params_string.format(**params))

    tbl1,tbl2,tbl3 = lamda.parse_lamda_datafile('molecule_ch3oh.inp')

    # no need to calculate populations when lines_mode=1
    if lvg:
        assert os.system('radmc3d calcpop writepop noscat') == 0

        shutil.copy('levelpop_ch3oh.bdat', 'levelpop_ch3oh_all.bdat')

        def read_levels(level_fn):
            with open(level_fn, 'rb') as fh:
                ftype, = struct.unpack('=q', fh.read(8))
                nrcells, = struct.unpack('=q', fh.read(8))
                nrlevels_subset, = struct.unpack('=q', fh.read(8))
                nlevels, = struct.unpack('=q', fh.read(8))
                levels = [struct.unpack('=q', fh.read(8)) for ii in range(nlevels)]
                data = np.fromfile(fh, dtype='float64', count=nlevels*nrlevels_subset)
            print(ftype, nrcells, nrlevels_subset, nlevels, levels)

            assert sz * sz * sz == nrlevels_subset

            return data.reshape([sz, sz, sz, nlevels])

        ch3oh_levels = read_levels('levelpop_ch3oh.bdat')

        #for level in (19,90,42,53):
        for trans in (240,241,242,251):

            level_U = int(tbl2['Upper'][tbl2['Transition'] == trans])
            level_L = int(tbl2['Lower'][tbl2['Transition'] == trans])
            level_U_label = tbl3['J'][tbl3['Level'] == level_U][0]
            level_L_label = tbl3['J'][tbl3['Level'] == level_L][0]

            fig2 = pl.figure(2)
            fig2.clf()
            pl.suptitle("{0} - {1}".format(level_U_label, level_L_label))
            ax1 = pl.subplot(2,2,1)
            # trans 240 = 19-13 = 4_22-312
            im = ax1.imshow(ch3oh_levels[:,:,int(sz/2),level_U-1], cmap='hot',
                            norm=matplotlib.colors.LogNorm())
            pl.colorbar(im, ax=ax1)
            ax2 = pl.subplot(1,2,2)
            ax2.semilogy(rr.ravel(), ch3oh_levels[:,:,:,level_U-1].ravel(), '.', alpha=0.25,
                         label=level_U_label)

            ax3 = pl.subplot(2,2,3)
            im = ax3.imshow(ch3oh_levels[:,:,int(sz/2),level_L-1], cmap='hot',
                            norm=matplotlib.colors.LogNorm())
            pl.colorbar(im, ax=ax3)
            ax4 = pl.subplot(1,2,2)
            ax4.semilogy(rr.ravel(), ch3oh_levels[:,:,:,level_L-1].ravel(), '.', alpha=0.25,
                         label=level_L_label)
            pl.legend(loc='best')

            pl.savefig("ch3oh_{0}-{1}_levelpops.png".format(level_U_label,level_L_label))


        # lines_mode = 50 means read from file
        with open('radmc3d.inp','w') as f:
            params['lines_mode'] = 50
            f.write(params_string.format(**params))



    # Debug: tau=1 surface should exist.
    os.system("radmc3d tausurf 1.0 lambda 1300")
    shutil.move('image.out', 'tausurf_1300um.out')
    im = radmc3dPy.image.readImage('tausurf_1300um.out')
    pl.figure(6).clf()
    pl.imshow(im.image[:,:,0])
    pl.savefig("tausurf_1300um.png")
    im.writeFits('tausurf_1300um.fits', fitsheadkeys={}, dpc=5400,
                 coord='19h23m43.963s +14d30m34.56s', overwrite=True)


    radmc3dPy.image.makeImage(nlam=100,
                              lambdarange=[500, 5000],
                              npix=50,
                              writepop=False,
                              noscat=True,
                              doppcatch=True,
                              nostar=False,
                              incl=0,
                              sizeau=10000)
    im = radmc3dPy.image.readImage('image.out')
    pl.figure(3).clf()
    pl.loglog(im.freq, im.image[25,25,0]*(im.freq/im.freq[0])**2,
              label="$\\nu^2$")
    pl.loglog(im.freq, im.image[25,25,:])
    pl.legend(loc='best')
    pl.savefig("SED.png")

    # debug: make an optical depth image
    os.system('radmc3d image npix 50 incl 0 sizeau 10000 noscat  pointau 0.0  0.0  0.0 fluxcons lambdarange 500 5000 tracetau')
    im = radmc3dPy.image.readImage('image.out')
    pl.figure(3).clf()
    pl.loglog(im.freq, im.image[25,25,:], label='center pixel')
    pl.loglog(im.freq, im.image[15,15,:], label='15,15')
    pl.ylabel("Optical Depth")
    pl.legend(loc='best')
    pl.savefig("optical_depth.png")


    radmc3dPy.image.makeImage(iline=240, widthkms=2,
                              linenlam=40,
                              nostar=False,
                              doppcatch=True,
                              noscat=True,
                              vkms=0,
                              writepop=False,
                              npix=50, incl=0,
                              #lambdarange=[wavelength_center*(1-10/3e5),
                              #             wavelength_center*(1+10/3e5)],
                              #nlam=40,
                              sizeau=10000)
    im = radmc3dPy.image.readImage('image.out')
    pl.figure(3).clf()
    pl.plot(im.freq, im.image[25,25,:], label='center pixel')
    pl.plot(im.freq, im.image[15,15,:], label='15,15')
    pl.ylabel("Line?")
    pl.legend(loc='best')
    pl.savefig("ch3oh_422_spectrum_X={0}.png".format(x_ch3oh))




    ###### TEST
    fig4 = pl.figure(4)
    fig5 = pl.figure(5)
    fig4.clf()
    fig5.clf()
    for iline in range(1,250,25):
        os.system("radmc3d image npix 1 incl 0 sizeau 10000 vkms 0 widthkms 2 noscat doppcatch pointau 0.0  0.0  0.0 fluxcons iline {0} > /dev/null".format(iline))
        im = radmc3dPy.image.readImage('image.out')
        fig4.gca().plot(im.freq, im.image.ravel(), label='{0}'.format(iline))
        fig5.gca().plot(np.linspace(-10,10,40), im.image.ravel(), label='{0}'.format(iline))
        if im.image.max() / im.image.min() < 1.1:
            print("peak is <1.1x min for line {0}".format(iline))
    pl.legend(loc='best')
    fig5.savefig('ch3oh_velocity_test_plot.png')





    shutil.move('image.out', 'ch3oh_422-312_image.out')
    im = radmc3dPy.image.readImage('ch3oh_422-312_image.out')
    im.writeFits('ch3oh_422-312_image_X={0}.fits'.format(x_ch3oh),
                 fitsheadkeys={}, dpc=5400,
                 coord='19h23m43.963s +14d30m34.56s', overwrite=True)


    radmc3dPy.image.makeImage(iline=242, widthkms=2,
                              linenlam=40,
                              nostar=False,
                              doppcatch=True,
                              noscat=True,
                              vkms=0,
                              writepop=False,
                              npix=50, incl=0,
                              #lambdarange=[wavelength_center*(1-10/3e5),
                              #             wavelength_center*(1+10/3e5)],
                              #nlam=40,
                              sizeau=10000)

    im = radmc3dPy.image.readImage('image.out')
    pl.figure(3).clf()
    pl.plot(im.freq, im.image[25,25,:], label='center pixel')
    pl.plot(im.freq, im.image[15,15,:], label='15,15')
    pl.ylabel("Line Brightness (Jy)")
    pl.legend(loc='best')
    pl.savefig("ch3oh_808_spectrum_X={0}.png".format(x_ch3oh))

    shutil.move('image.out', 'ch3oh_808-716_image.out')
    im = radmc3dPy.image.readImage('ch3oh_808-716_image.out')
    im.writeFits('ch3oh_808-716_image_X={0}.fits'.format(x_ch3oh), fitsheadkeys={}, dpc=5400,
                 coord='19h23m43.963s +14d30m34.56s', overwrite=True)


    from astropy.io import fits
    from astropy import wcs

    def convert_to_K(radmc_fits_img, distance=5400*u.pc):
        fh = fits.open(radmc_fits_img)
        mywcs = wcs.WCS(fh[0].header)
        pix_area = np.abs(mywcs.celestial.pixel_scale_matrix.diagonal().prod()) * u.deg**2
        conv = u.Jy.to(u.K, equivalencies=u.brightness_temperature(pix_area, mywcs.wcs.crval[2]*u.Hz))
        fh[0].data *= conv
        fh[0].header['BUNIT'] = 'K'
        return fh


    for trans in (240,241,242,251):

        widthkms = 2
        linenlam = 40
        vkms = 0
        radmc3dPy.image.makeImage(iline=trans, widthkms=widthkms,
                                  linenlam=linenlam,
                                  nostar=False,
                                  doppcatch=True,
                                  noscat=True,
                                  vkms=vkms,
                                  writepop=False,
                                  npix=50, incl=0,
                                  #lambdarange=[wavelength_center*(1-10/3e5),
                                  #             wavelength_center*(1+10/3e5)],
                                  #nlam=40,
                                  sizeau=10000)

        level_U = int(tbl2['Upper'][tbl2['Transition'] == trans])
        level_L = int(tbl2['Lower'][tbl2['Transition'] == trans])
        level_U_label = tbl3['J'][tbl3['Level'] == level_U][0]
        level_L_label = tbl3['J'][tbl3['Level'] == level_L][0]

        prefix = "ch3oh_{0}-{1}".format(level_U_label,level_L_label)

        im = radmc3dPy.image.readImage('image.out')

        shutil.move('image.out', '{0}_image_X={1}.out'.format(prefix, x_ch3oh))
        im = radmc3dPy.image.readImage('{0}_image_X={1}.out'.format(prefix, x_ch3oh))
        im.writeFits('{0}_image_X={1}.fits'.format(prefix, x_ch3oh), fitsheadkeys={}, dpc=5400,
                     coord='19h23m43.963s +14d30m34.56s', overwrite=True)

        fh = convert_to_K('{0}_image_X={1}.fits'.format(prefix, x_ch3oh))
        fh.writeto("{0}_image_K_X={1}.fits".format(prefix, x_ch3oh), clobber=True)

        img = fh[0].data
        velo = np.linspace(vkms-widthkms/2, vkms+widthkms/2, linenlam)
        pl.figure(3).clf()
        pl.plot(velo, img[:,25,25], label='center pixel')
        pl.plot(velo, img[:,15,15], label='15,15')
        pl.plot(velo, img[:,5,5], label='5,5')
        pl.ylabel("Line Brightness (K)")
        pl.xlabel("Velocity (km/s)")
        pl.legend(loc='best')
        pl.savefig("{0}_spectrum_X={1}.png".format(prefix, x_ch3oh))


    #radmc3dPy.image.makeImage(iline=240, widthkms=5, linenlam=40, nostar=True,
    #                          npix=500, incl=0,
    #                          lambdarange=[wavelength_center*(1-10/3e5),
    #                                       wavelength_center*(1+10/3e5)],
    #                          nlam=20, sizeau=10000)
    #shutil.move('image.out', 'ch3oh_422-312_image.out')
    #im = radmc3dPy.image.readImage('ch3oh_422-312_image.out')
    #im.writeFits('ch3oh_422-312_image.fits', fitsheadkeys={}, dpc=5400, coord='19h23m43.963s +14d30m34.56s', overwrite=True)
    # radmc3dPy.image.makeImage(iline=3, widthkms=5, linenlam=40, nostar=True)
    # shutil.move('image.out', 'h2co_303-202_image.out')
    # radmc3dPy.image.makeImage(iline=13, widthkms=5, linenlam=40, nostar=True)
    # shutil.move('image.out', 'h2co_321-220_image.out')
    # radmc3dPy.image.makeImage(iline=2, widthkms=5, linenlam=40, nostar=True)
    # shutil.move('image.out', 'co_2-1_image.out')
    # radmc3dPy.image.makeImage(iline=1, widthkms=5, linenlam=40, nostar=True)
    # shutil.move('image.out', 'co_1-0_image.out')
    # #radmc3d image iline 1 widthkms 10 linenlam 40 linelist nostar
