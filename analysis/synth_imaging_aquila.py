"""
Synthetic imaging / simulated observing of the aquila Herschel map at
250um (selected for not-quite resolution matching) observed at the distance and
sensitivity of W51
"""
from astropy.io import fits
from astropy import units as u
from astropy import constants
import paths
import requests
import os
import numpy as np
import dust_emissivity

daquila = 260.
dw51 = 5410.
distance_scaling = daquila/dw51
freq_alma = 226.6*u.GHz
wave_alma = constants.c/(freq_alma)
#wave_herschel = 250*u.um
#freq_herschel = constants.c/wave_herschel
#nu = u.Quantity([freq_herschel, freq_alma])
#flux_ = dust_emissivity.blackbody.modified_blackbody(nu, temperature=20*u.K, beta=1.5)
#alpha = np.log(flux_[0]/flux_[1])/np.log(nu[0]/nu[1])
#flux_scaling = (wave_herschel/wave_alma).decompose().value**alpha

# Herschel 250um has ~18" beam
#in_bm = (2*np.pi*(18.*u.arcsec/206265./2.35)**2)
# but we want to keep constant MJy/sr (surface brightness) out to a larger
# distance, then convert that to Jy/(ALMA beam)
# use 18.2" from Konyves 2015
herschel_resoln = 18.2 # arcsec
fwhm_scale = (8*np.log(2))**0.5
out_bm = (2*np.pi*((herschel_resoln*distance_scaling)*u.arcsec/fwhm_scale)**2)
#MJySr_to_JyBm = (1*u.MJy/u.sr).to(u.Jy/in_bm).value
MJySr_to_JyBm = (1*u.MJy/u.sr).to(u.Jy/out_bm).value

d_aquila = paths.dpath('not_w51/HGBS_aquilaM2_column_density_map.fits')
t_aquila = paths.dpath('not_w51/HGBS_aquilaM2_dust_temperature_map.fits')
aquila_rescaled = paths.dpath('not_w51/aquila04-250-rescaled.fits')

if not os.path.exists(d_aquila):
    result = requests.get('http://www.herschel.fr/Phocea/file.php?class=astimg&file=66/HGBS_aquilaM2_column_density_map.fits.gz')
    with open(d_aquila,'w') as f:
        f.write(result.content)

if not os.path.exists(t_aquila):
    result = requests.get('http://www.herschel.fr/Phocea/file.php?class=astimg&file=66/HGBS_aquilaM2_dust_temperature_map.fits.gz')
    with open(t_aquila,'w') as f:
        f.write(result.content)

# scale data to 1100um
ffile = fits.open(t_aquila)
orig_pix_area = np.abs(ffile[0].header['CDELT1'] * ffile[0].header['CDELT2']) * u.deg**2
ffile[0].header['CRVAL1'] = 290.92402
ffile[0].header['CRPIX1'] = 2339
ffile[0].header['CUNIT1'] = 'deg'
ffile[0].header['CDELT1'] = ffile[0].header['CDELT1'] * distance_scaling
ffile[0].header['CTYPE1'] = 'RA---TAN'
ffile[0].header['CRVAL2'] = 14.512736
ffile[0].header['CRPIX2'] = 3711
ffile[0].header['CUNIT2'] = 'deg'
ffile[0].header['CDELT2'] = ffile[0].header['CDELT2'] * distance_scaling
ffile[0].header['CTYPE2'] = 'DEC--TAN'
ffile[0].header['CRPIX3'] = 1
ffile[0].header['CRVAL3'] = 2.33946806e+11
ffile[0].header['CUNIT3'] = 'Hz'
ffile[0].header['CDELT3'] = 1e9
ffile[0].header['CTYPE3'] = 'FREQ' # MAYBE frequency?
ffile[0].header['CRPIX4'] = 1
ffile[0].header['CRVAL4'] = 1
ffile[0].header['CUNIT4'] = ''
ffile[0].header['CDELT4'] = 0
ffile[0].header['CTYPE4'] = 'STOKES'
ffile[0].header['RESTFRQ'] = 2.33946806e+11
ffile[0].header['BUNIT'] = 'Jy/pixel'
dustcol = fits.getdata(d_aquila)
dusttem = fits.getdata(t_aquila)
surfbright = dust_emissivity.blackbody.modified_blackbody(freq_alma,
                                                          dusttem*u.K,
                                                          beta=1.5,
                                                          column=dustcol*u.cm**-2,) / u.sr
np.testing.assert_approx_equal(dustcol[3687,2105], 5.2438327e+21, )
np.testing.assert_approx_equal(dusttem[3687,2105], 34.599281)
testval = dust_emissivity.blackbody.modified_blackbody(nu=226.6*u.GHz, temperature=34.599281*u.K, column=5.2438327e+21*u.cm**-2)
np.testing.assert_almost_equal(surfbright[3687,2015].value, testval.value)
assert not np.all(surfbright[surfbright==surfbright]==0)
pixel_area = np.abs(ffile[0].header['CDELT1'] * ffile[0].header['CDELT2']) * u.deg**2
fluxdens = surfbright.to(u.MJy/u.sr).value * MJySr_to_JyBm
#ffile[0].data = np.expand_dims(ffile[0].data, axis=0)

# divide by ppbeam for jy/pixel
ppbeam = (out_bm / pixel_area).decompose().value
ffile[0].data = fluxdens / ppbeam

ffile.writeto(aquila_rescaled, clobber=True)
hdr = ffile[0].header

aquila_casa_image = 'aquila_dusttem_to_w51.image'
# doesn't always work: unreliable = don't use.
# rmresult = rmtables([aquila_casa_image])
os.system('rm -rf {0}'.format(aquila_casa_image))
importfits(fitsimage=aquila_rescaled,
           imagename=aquila_casa_image,
           overwrite=True,
           defaultaxes=True,#['RA---TAN','DEC--TAN','FREQUENCY','STOKES'],
           defaultaxesvalues=['2.909234541667E+02deg',
                              '1.451177222222E+01deg',
                              '233.9468GHz','I'],
           # 18" = 1.22 lambda/D
           beam=["{0}deg".format(herschel_resoln/3600.*distance_scaling),
                 "{0}deg".format(herschel_resoln/3600.*distance_scaling),
                 "0deg"],
          )
print("FITS CDELT1={0}, CDELT2={1}".format(hdr['CDELT1'], hdr['CDELT2']))
print("image CDELT1={0[value]}{0[unit]}, CDELT2={1[value]}{1[unit]}"
      .format(imhead(imagename=aquila_casa_image, mode='get', hdkey='CDELT1'),
              imhead(imagename=aquila_casa_image, mode='get',
                     hdkey='CDELT2'),)
     )

#imhead(aquila_casa_image, mode='put', hdkey='CDELT1', hdvalue={'value':hdr['CDELT1'], 'unit':'deg'})
#imhead(aquila_casa_image, mode='put', hdkey='CDELT2', hdvalue={'value':hdr['CDELT2'], 'unit':'deg'})
#imhead(aquila_casa_image, mode='put', hdkey='CRVAL1', hdvalue={'value':hdr['CRVAL1'], 'unit':'deg'})
#imhead(aquila_casa_image, mode='put', hdkey='CRVAL2', hdvalue={'value':hdr['CRVAL2'], 'unit':'deg'})
#imhead(aquila_casa_image, mode='put', hdkey='CRPIX1', hdvalue=hdr['CRPIX1'])
#imhead(aquila_casa_image, mode='put', hdkey='CRPIX2', hdvalue=hdr['CRPIX2'])
#imhead(aquila_casa_image, mode='put', hdkey='CTYPE3', hdvalue='FREQ')
# can convert units and use this
#imhead(aquila_casa_image, mode='put', hdkey='BUNIT', hdvalue='Jy/beam')
#imhead(aquila_casa_image, mode='put', hdkey='BUNIT', hdvalue='Jy/pixel') #WRONG!!!
exportfits(aquila_casa_image, aquila_casa_image+".fits",
           overwrite=True)
hdr = fits.getheader(aquila_casa_image+".fits")
print("CDELT1={0}, CDELT2={1}".format(hdr['CDELT1'], hdr['CDELT2']))

os.system('rm -rf {0}'.format(aquila_casa_image))
# try re-importing (this definitely should not fix any outstanding issues)
importfits(fitsimage=aquila_rescaled,
           imagename=aquila_casa_image,
           overwrite=True,
           # The beam is OBVIOUSLY AND CORRECTLY in the header.
           # but if you leave this out, CASA complains loudly
           beam=["{0}deg".format(herschel_resoln/3600.*distance_scaling),
                 "{0}deg".format(herschel_resoln/3600.*distance_scaling),
                 "0deg"],
          )
ia.open(aquila_casa_image)
print("BUNIT: ",ia.summary()['unit'])
print("brightnessunit: ",ia.brightnessunit())
print("Peak value: {max}  Total: {sum}".format(**ia.statistics()))
ia.close()

#sm.openfromms("w51_contvis_selfcal_0.ms")
sm.openfromms("continuum_7m12m_noflag.ms")
#sm.openfromms("w51_test_onechan.ms")
sm.setvp()
ia.open(aquila_casa_image)
print("After setvp: ")
print("BUNIT: ",ia.summary()['unit'])
print("brightnessunit: ",ia.brightnessunit())
print("Peak value: {max}  Total: {sum}".format(**ia.statistics()))
ia.close()


success = sm.predict(aquila_casa_image)

ia.open(aquila_casa_image)
print("After prediction: ")
print("BUNIT: ",ia.summary()['unit'])
print("brightnessunit: ",ia.brightnessunit())
print("Peak value: {max}  Total: {sum}".format(**ia.statistics()))
ia.close()

# TODO: get these from ASDM_CALWVR and WEATHER
success2 = sm.setnoise(mode='tsys-atm', relhum=60.0, pwv='2mm', tatmos=265.0, )

ia.open(aquila_casa_image)
print("After setnoise: ")
print("BUNIT: ",ia.summary()['unit'])
print("brightnessunit: ",ia.brightnessunit())
print("Peak value: {max}  Total: {sum}".format(**ia.statistics()))
ia.close()

success3 = sm.corrupt()

ia.open(aquila_casa_image)
print("After corrupt: ")
print("BUNIT: ",ia.summary()['unit'])
print("brightnessunit: ",ia.brightnessunit())
print("Peak value: {max}  Total: {sum}".format(**ia.statistics()))
ia.close()


sm.done()
sm.close()

ia.open(aquila_casa_image)
print("After prediction & corruption & close: ")
print("BUNIT: ",ia.summary()['unit'])
print("brightnessunit: ",ia.brightnessunit())
print("Peak value: {max}  Total: {sum}".format(**ia.statistics()))
ia.close()



os.system('rm -rf aquila_dusttem_model.ms')
assert split(vis="continuum_7m12m_noflag.ms", outputvis="aquila_dusttem_model.ms",
             datacolumn='data')

phasecenter = 'J2000 19h23m41.580 +14d30m41.37'

delmod(vis='aquila_dusttem_model.ms')
os.system('rm -rf aquila_dusttem_model_tclean_dirty*')
tclean(vis='aquila_dusttem_model.ms',
       imagename='aquila_dusttem_model_tclean_dirty',
       field='',
       spw='',
       specmode='mfs',
       deconvolver='clark',
       imsize = 1536/2,
       cell= '0.2arcsec',
       weighting = 'natural',
       phasecenter=phasecenter,
       #scales=[0,3,9,27,81],
       robust = 2.0,
       uvrange='0~500m',
       niter = 0,
       threshold = '1.0mJy',
       interactive = False,
       gridder = 'mosaic',
       savemodel='none',
       )
exportfits(imagename='aquila_dusttem_model_tclean_dirty.residual', fitsimage='aquila_dusttem_model_tclean_dirty.image.fits',  dropdeg=True, overwrite=True)


os.system('rm -rf aquila_dusttem_model_tclean_clean*')
tclean(vis='aquila_dusttem_model.ms',
       imagename='aquila_dusttem_model_tclean_clean',
       field='',
       spw='',
       specmode='mfs',
       deconvolver='clark',
       imsize = 1536/2,
       cell= '0.2arcsec',
       weighting = 'natural',
       phasecenter=phasecenter,
       #scales=[0,3,9,27,81],
       robust = 2.0,
       uvrange='0~500m',
       niter = 50000,
       threshold = '0.5mJy',
       interactive = False,
       gridder = 'mosaic',
       savemodel='none',
       )
exportfits(imagename='aquila_dusttem_model_tclean_clean.image', fitsimage='aquila_dusttem_model_tclean_clean.image.fits',  dropdeg=True, overwrite=True)
exportfits(imagename='aquila_dusttem_model_tclean_clean.model', fitsimage='aquila_dusttem_model_tclean_clean.model.fits',  dropdeg=True, overwrite=True)

os.system('rm -rf aquila_dusttem_model_tclean_msclean*')
tclean(vis='aquila_dusttem_model.ms',
       imagename='aquila_dusttem_model_tclean_msclean',
       field='',
       spw='',
       specmode='mfs',
       deconvolver='multiscale',
       imsize = [1536,1536],
       cell= '0.1arcsec',
       weighting = 'uniform',
       phasecenter=phasecenter,
       scales=[0,3,9,27],
       robust = 2.0,
       niter = 50000,
       uvrange='0~500m',
       threshold = '0.5mJy',
       interactive = False,
       gridder = 'mosaic',
       savemodel='none',
       )
exportfits(imagename='aquila_dusttem_model_tclean_msclean.image', fitsimage='aquila_dusttem_model_tclean_msclean.image.fits',  dropdeg=True, overwrite=True)
exportfits(imagename='aquila_dusttem_model_tclean_msclean.model', fitsimage='aquila_dusttem_model_tclean_msclean.model.fits',  dropdeg=True, overwrite=True)
