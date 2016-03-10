"""
Synthetic imaging / simulated observing of the NGC 1333/Perseus Herschel map at
250um (selected for resolution matching) observed at the distance and
sensitivity of W51

This time, it is made artificially 10x brighter in a region that looks
kind of similar to W51
"""
from astropy.io import fits
from astropy import units as u
from astropy import constants
import paths
import requests
import os
import numpy as np
import dust_emissivity

dperseus = 140.
dw51 = 5410.
distance_scaling = dperseus/dw51
freq_alma = 225*u.GHz
wave_alma = constants.c/(freq_alma)
wave_herschel = 250*u.um
freq_herschel = constants.c/wave_herschel
nu = u.Quantity([freq_herschel, freq_alma])
flux = dust_emissivity.blackbody.modified_blackbody(nu, temperature=20*u.K, beta=1.5)
alpha = np.log(flux[0]/flux[1])/np.log(nu[0]/nu[1])
flux_scaling = (wave_herschel/wave_alma).decompose().value**alpha
flux_scaling *= 100

# Herschel 250um has ~18" beam
#in_bm = (2*np.pi*(18.*u.arcsec/206265./2.35)**2)
# but we want to keep constant MJy/sr (surface brightness) out to a larger
# distance, then convert that to Jy/(ALMA beam)
out_bm = (2*np.pi*((18.*distance_scaling)*u.arcsec/2.35)**2)
#MJySr_to_JyBm = (1*u.MJy/u.sr).to(u.Jy/in_bm).value
MJySr_to_JyBm = (1*u.MJy/u.sr).to(u.Jy/out_bm).value

im_perseus = paths.dpath('perseus04-250.fits.gz')
perseus_rescaled = paths.dpath('perseus04-250-rescaled_2.fits')

if not os.path.exists(im_perseus):
    result = requests.get('http://www.herschel.fr/cea/gouldbelt/en/archives/perseus04/perseus04-250.fits.gz')
    with open(im_perseus,'w') as f:
        f.write(result.content)

# scale data to 1100um
ffile = fits.open(im_perseus)
ffile[0].header['CRVAL1'] = 290.92402
ffile[0].header['CRPIX1'] = 580
ffile[0].header['CUNIT1'] = 'deg'
ffile[0].header['CDELT1'] = ffile[0].header['CDELT1'] * distance_scaling
ffile[0].header['CTYPE1'] = 'RA---TAN'
ffile[0].header['CRVAL2'] = 14.512736
ffile[0].header['CRPIX2'] = 1442
ffile[0].header['CUNIT2'] = 'deg'
ffile[0].header['CDELT2'] = ffile[0].header['CDELT2'] * distance_scaling
ffile[0].header['CTYPE2'] = 'DEC--TAN'
ffile[0].header['CRPIX3'] = 1
ffile[0].header['CRVAL3'] = 2.33946806e+11
ffile[0].header['CUNIT3'] = 'Hz'
ffile[0].header['CDELT3'] = 1e9
ffile[0].header['CTYPE3'] = 'FREQ'
ffile[0].header['CRPIX4'] = 1
ffile[0].header['CRVAL4'] = 1
ffile[0].header['CUNIT4'] = ''
ffile[0].header['CDELT4'] = 0
ffile[0].header['CTYPE4'] = 'STOKES'
ffile[0].header['RESTFRQ'] = 2.33946806e+11
ffile[0].header['BUNIT'] = 'Jy/beam'
ffile[0].data *= flux_scaling * MJySr_to_JyBm
#ffile[0].data = np.expand_dims(ffile[0].data, axis=0)
ffile.writeto(perseus_rescaled, clobber=True)
hdr = ffile[0].header

perseus_casa_image = 'perseus_250_to_w51_2.image'
# doesn't always work: unreliable = don't use.
# rmresult = rmtables([perseus_casa_image])
os.system('rm -rf {0}'.format(perseus_casa_image))
importfits(fitsimage=perseus_rescaled,
           imagename=perseus_casa_image,
           overwrite=True,
           defaultaxes=True,#['RA---TAN','DEC--TAN','FREQUENCY','STOKES'],
           defaultaxesvalues=['2.909234541667E+02deg',
                              '1.451177222222E+01deg',
                              '233.9468GHz','I'],
           # 18" = 1.22 lambda/D
           beam=["{0}deg".format(18/3600.*distance_scaling),
                 "{0}deg".format(18/3600.*distance_scaling),
                 "0deg"],
          )
print("FITS CDELT1={0}, CDELT2={1}".format(hdr['CDELT1'], hdr['CDELT2']))
print("image CDELT1={0[value]}{0[unit]}, CDELT2={1[value]}{1[unit]}"
      .format(imhead(imagename=perseus_casa_image, mode='get', hdkey='CDELT1'),
              imhead(imagename=perseus_casa_image, mode='get',
                     hdkey='CDELT2'),)
     )

#imhead(perseus_casa_image, mode='put', hdkey='CDELT1', hdvalue={'value':hdr['CDELT1'], 'unit':'deg'})
#imhead(perseus_casa_image, mode='put', hdkey='CDELT2', hdvalue={'value':hdr['CDELT2'], 'unit':'deg'})
#imhead(perseus_casa_image, mode='put', hdkey='CRVAL1', hdvalue={'value':hdr['CRVAL1'], 'unit':'deg'})
#imhead(perseus_casa_image, mode='put', hdkey='CRVAL2', hdvalue={'value':hdr['CRVAL2'], 'unit':'deg'})
#imhead(perseus_casa_image, mode='put', hdkey='CRPIX1', hdvalue=hdr['CRPIX1'])
#imhead(perseus_casa_image, mode='put', hdkey='CRPIX2', hdvalue=hdr['CRPIX2'])
#imhead(perseus_casa_image, mode='put', hdkey='CTYPE3', hdvalue='FREQ')
# can convert units and use this
#imhead(perseus_casa_image, mode='put', hdkey='BUNIT', hdvalue='Jy/beam')
#imhead(perseus_casa_image, mode='put', hdkey='BUNIT', hdvalue='Jy/pixel') #WRONG!!!
exportfits(perseus_casa_image, perseus_casa_image+".fits",
           overwrite=True)
hdr = fits.getheader(perseus_casa_image+".fits")
print("CDELT1={0}, CDELT2={1}".format(hdr['CDELT1'], hdr['CDELT2']))

os.system('rm -rf {0}'.format(perseus_casa_image))
# try re-importing (this definitely should not fix any outstanding issues)
importfits(fitsimage=perseus_rescaled,
           imagename=perseus_casa_image,
           overwrite=True,
           # The beam is OBVIOUSLY AND CORRECTLY in the header.
           # but if you leave this out, CASA complains loudly
           beam=["{0}deg".format(18/3600.*distance_scaling),
                 "{0}deg".format(18/3600.*distance_scaling),
                 "0deg"],
          )
ia.open(perseus_casa_image)
print("BUNIT: ",ia.summary()['unit'])
ia.close()

#sm.openfromms("w51_contvis_selfcal_0.ms")
sm.openfromms("continuum_7m12m_noflag.ms")
#sm.openfromms("w51_test_onechan.ms")
sm.setvp()
success = sm.predict(perseus_casa_image)
# TODO: get these from ASDM_CALWVR and WEATHER
success2 = sm.setnoise(mode='tsys-atm', relhum=60.0, pwv='2mm', tatmos=265.0, )
success3 = sm.corrupt()
sm.done()
sm.close()

# problem:
# plotms(vis='continuum_7m12m_noflag.ms', xaxis='uvdist', ydatacolumn='model')

base_name = 'perseus_250_2'

os.system('rm -rf {0}_model.ms'.format(base_name))
assert split(vis="continuum_7m12m_noflag.ms", outputvis="{0}_model.ms".format(base_name),
             datacolumn='data')

phasecenter = 'J2000 19h23m41.580 +14d30m41.37'

os.system('rm -rf {0}_model_tclean_dirty*'.format(base_name))
tclean(vis='{0}_model.ms'.format(base_name),
       imagename='{0}_model_tclean_dirty'.format(base_name),
       field='',
       spw='',
       specmode='mfs',
       deconvolver='clark',
       imsize = [1536,1536],
       cell= '0.1arcsec',
       weighting = 'uniform',
       phasecenter=phasecenter,
       #scales=[0,3,9,27,81],
       robust = -2.0,
       niter = 0,
       threshold = '1.0mJy',
       interactive = False,
       gridder = 'mosaic',
       savemodel='none',
       )
exportfits(imagename='{0}_model_tclean_dirty.residual'.format(base_name), fitsimage='{0}_model_tclean_dirty.image.fits'.format(base_name),  dropdeg=True, overwrite=True)

# dirtyimage = '{0}_model_tclean_dirty.residual'.format(base_name)
# assert ia.open(dirtyimage)
# assert ia.calcmask(mask=dirtyimage+" > 0.1", name='dirty_mask_100mJy')
# assert ia.close()
# makemask(mode='copy', inpimage=dirtyimage,
#          inpmask=dirtyimage+":dirty_mask_100mJy", output='dirty_100mJy.mask',
#          overwrite=True)
# exportfits('dirty_100mJy.mask', 'dirty_100mJy.mask.fits', dropdeg=True, overwrite=True)
# 
# os.system('rm -rf {0}_model_tclean_clean_masked*'.format(base_name))
# tclean(vis='{0}_model.ms'.format(base_name),
#        imagename='{0}_model_tclean_clean_masked'.format(base_name),
#        field='',
#        spw='',
#        specmode='mfs',
#        deconvolver='clark',
#        imsize = [1536,1536],
#        cell= '0.1arcsec',
#        weighting = 'uniform',
#        phasecenter=phasecenter,
#        #scales=[0,3,9,27,81],
#        robust = -2.0,
#        niter = 50000,
#        threshold = '7.0mJy',
#        interactive = False,
#        gridder = 'mosaic',
#        savemodel='none',
#        mask='dirty_100mJy.mask',
#        )
# exportfits(imagename='{0}_model_tclean_clean.image'.format(base_name), fitsimage='{0}_model_tclean_clean.image.fits'.format(base_name),  dropdeg=True, overwrite=True)
# exportfits(imagename='{0}_model_tclean_clean.model'.format(base_name), fitsimage='{0}_model_tclean_clean.model.fits'.format(base_name),  dropdeg=True, overwrite=True)


os.system('rm -rf {0}_model_tclean_clean*'.format(base_name))
tclean(vis='{0}_model.ms'.format(base_name),
       imagename='{0}_model_tclean_clean'.format(base_name),
       field='',
       spw='',
       specmode='mfs',
       deconvolver='clark',
       imsize = [1536,1536],
       cell= '0.1arcsec',
       weighting = 'uniform',
       phasecenter=phasecenter,
       #scales=[0,3,9,27,81],
       robust = -2.0,
       niter = 50000,
       threshold = '2.0mJy',
       interactive = False,
       gridder = 'mosaic',
       savemodel='none',
       )
exportfits(imagename='{0}_model_tclean_clean.image'.format(base_name), fitsimage='{0}_model_tclean_clean.image.fits'.format(base_name),  dropdeg=True, overwrite=True)
exportfits(imagename='{0}_model_tclean_clean.model'.format(base_name), fitsimage='{0}_model_tclean_clean.model.fits'.format(base_name),  dropdeg=True, overwrite=True)

os.system('rm -rf {0}_model_tclean_msclean*'.format(base_name))
tclean(vis='{0}_model.ms'.format(base_name),
       imagename='{0}_model_tclean_msclean'.format(base_name),
       field='',
       spw='',
       specmode='mfs',
       deconvolver='multiscale',
       imsize = [1536,1536],
       cell= '0.1arcsec',
       weighting = 'uniform',
       phasecenter=phasecenter,
       scales=[0,3,9,27],
       robust = -2.0,
       niter = 50000,
       threshold = '2.0mJy',
       interactive = False,
       gridder = 'mosaic',
       savemodel='none',
       )
exportfits(imagename='{0}_model_tclean_msclean.image'.format(base_name), fitsimage='{0}_model_tclean_msclean.image.fits'.format(base_name),  dropdeg=True, overwrite=True)
exportfits(imagename='{0}_model_tclean_msclean.model'.format(base_name), fitsimage='{0}_model_tclean_msclean.model.fits'.format(base_name),  dropdeg=True, overwrite=True)
