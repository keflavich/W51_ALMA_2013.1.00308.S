"""
Synthetic imaging / simulated observing of the NGC 1333/Perseus Herschel map at
250um (selected for resolution matching) observed at the distance and
sensitivity of W51
"""
from astropy.io import fits
import paths
import requests
import os
import numpy as np

dperseus = 140.
dw51 = 5410.
distance_scaling = dperseus/dw51
flux_scaling = (250/1100.)**3.5

im_perseus = paths.dpath('perseus04-250.fits.gz')
perseus_rescaled = paths.dpath('perseus04-250-rescaled.fits')

if not os.path.exists(im_perseus):
    result = requests.get('http://www.herschel.fr/cea/gouldbelt/en/archives/perseus04/perseus04-250.fits.gz')
    with open(im_perseus,'w') as f:
        f.write(result.content)

# scale data to 1100um
ffile = fits.open(im_perseus)
ffile[0].header['CRVAL1'] = 290.92402
ffile[0].header['CRPIX1'] = 1100
ffile[0].header['CUNIT1'] = 'deg'
ffile[0].header['CDELT1'] = ffile[0].header['CDELT1'] * distance_scaling
ffile[0].header['CTYPE1'] = 'RA---TAN'
ffile[0].header['CRVAL2'] = 14.512736
ffile[0].header['CRPIX2'] = 1553
ffile[0].header['CUNIT2'] = 'deg'
ffile[0].header['CDELT2'] = ffile[0].header['CDELT2'] * distance_scaling
ffile[0].header['CTYPE2'] = 'DEC--TAN'
ffile[0].header['CRPIX3'] = 1
ffile[0].header['CRVAL3'] = 220.0e9
ffile[0].header['CUNIT3'] = 'Hz'
ffile[0].header['CDELT3'] = 1e9
ffile[0].header['RESTFRQ'] = 220.0e9
ffile[0].data *= flux_scaling
ffile[0].data = np.expand_dims(ffile[0].data, axis=0)
ffile.writeto(perseus_rescaled, clobber=True)
hdr = ffile[0].header

perseus_casa_image = 'perseus_250_to_w51.image'
# doesn't always work: unreliable = don't use.
# rmresult = rmtables([perseus_casa_image])
os.system('rm -rf {0}'.format(perseus_casa_image))
importfits(fitsimage=perseus_rescaled,
           imagename=perseus_casa_image,
           overwrite=True,
           defaultaxes=True,
           defaultaxesvalues=['2.909234541667E+02deg',
                              '1.451177222222E+01deg',
                              '220GHz','I'],
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
#imhead(perseus_casa_image, mode='put', hdkey='BUNIT', hdvalue='Jy/beam') #WRONG!!!
imhead(perseus_casa_image, mode='put', hdkey='BUNIT', hdvalue='Jy/pixel') #WRONG!!!
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
           #beam=["{0}deg".format(18/3600.*distance_scaling),
           #      "{0}deg".format(18/3600.*distance_scaling),
           #      "0deg"],
          )

sm.openfromms("w51_contvis_selfcal_0.ms")
sm.openfromms("continuum_7m12m_noflag.ms")
sm.setvp()
success = sm.predict(perseus_casa_image)
sm.done()
sm.close()

#im.open("continuum_7m12m_noflag.ms")
#im.defineimage(nx=hdr['NAXIS1'],
#               cellx='{0}arcsec'.format(hdr['CDELT2']*3600),
#               mode='mfs'),
#
##im.makemodelfromsd(sdimage=perseus_casa_image, modelimage='modelimage_'+perseus_casa_image)
#im.ft(model=perseus_casa_image)
#im.close()

#ft(vis='continuum_7m12m_noflag.ms',
#   model=perseus_casa_image, usescratch=True)


#inbright = distance_scaling
#indirection = "J2000 19h23m41.509 +14d30m41.46"
#incell = "{0}arcsec".format(np.abs(ffile[0].header['CDELT1']*3600) * dperseus/dw51)
#incenter = '220GHz'
#inwidth = '2GHz'

