import os
if not os.path.exists('test_frequency.image.fits'):
    os.system('rm -rf w51_test_small.ms')
    concat(vis=['w51_test_small_spw3.ms','w51_test_small_spw7.ms'],
           concatvis='w51_test_small.ms')

    os.system('rm -rf test_frequency.*')
    clean(vis='w51_test_small.ms', imagename="test_frequency",
          field="", spw='', mode='channel', outframe='LSRK',
          interpolation='linear', imagermode='mosaic', interactive=False,
          niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
          weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
          pbcor=False, usescratch=True, robust=1.0)
    exportfits('test_frequency.image', 'test_frequency.image.fits', dropdeg=True, overwrite=True)
    exportfits('test_frequency.model', 'test_frequency.model.fits', dropdeg=True, overwrite=True)
    for suffix in ('image','model','flux','psf','residual'):
        os.system('rm -rf test_frequency.{0}'.format(suffix))

import numpy as np
import spectral_cube
from astropy import units as u
cube = spectral_cube.SpectralCube.read('test_frequency.image.fits').with_spectral_unit(u.Hz)
cont = cube.min(axis=0)
contsub = cube-cont
# temporary hack for issue #251
contsub._data = contsub._data.value
contsub = contsub.with_mask(contsub>15*u.mJy)
ppbeam = np.abs((cube.beam.sr / (cube.wcs.pixel_scale_matrix[0,0]*cube.wcs.pixel_scale_matrix[1,1]*u.deg**2)).decompose())
hdu = contsub.hdu
# this scaling may be necessary if setjy is used
# I found that the cube version was about 70x too low, maybe worse...
#hdu.data *= ppbeam # because apparently CASA divides by this?
hdu.writeto('test_contsub_cube.fits', clobber=True)
header = contsub.header
importfits('test_contsub_cube.fits', 'test_contsub_cube.image',
           overwrite=True, defaultaxes=T,
           defaultaxesvalues=[header['CRVAL1'], header['CRVAL2'],
                              header['CRVAL3'], 'I'])

os.system('rm -rf w51_test_small_linecubesub.ms')
split('w51_test_small.ms', 'w51_test_small_linecubesub.ms', datacolumn='data')

im.open('w51_test_small_linecubesub.ms')
from astropy import coordinates
c = coordinates.SkyCoord(header['CRVAL1'], header['CRVAL2'], unit=('deg','deg'), frame='fk5')
im.defineimage(nx=cube.shape[1],
               cellx='{0}arcsec'.format(header['CDELT2']*3600),
               mode='channel',
               nchan=cube.shape[0],
               phasecenter='J2000 {0} {1}'.format(c.ra.to_string(),
                                                  c.dec.to_string())),
os.system('rm -rf model_contsubcube')
im.makemodelfromsd(sdimage='test_contsub_cube.image', modelimage='model_contsubcube')
# apparently images can't be applied across frequency, but models can
#im.ft(model='test_continuum_min.image')
im.ft(model='model_contsubcube')
im.close()
# setjy does some strange scaling
#setjy(vis='w51_test_small_linecubesub.ms', model='model_test', usescratch=True,
#      field='', standard='manual')
#
# ft claims there is no frequency overlap between the image and the data
#ft(vis='w51_test_small_linecubesub.ms', model='test_continuum_min.image', usescratch=True, nterms=1)

uvsub('w51_test_small_linecubesub.ms')

os.system('rm -rf test_frequency_linecubesub.*')
clean(vis='w51_test_small_linecubesub.ms', imagename="test_frequency_linecubesub",
      field="", spw='', mode='channel', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_frequency_linecubesub.image', 'test_frequency_linecubesub.image.fits', dropdeg=True, overwrite=True)
exportfits('test_frequency_linecubesub.model', 'test_frequency_linecubesub.model.fits', dropdeg=True, overwrite=True)
for suffix in ('image','model','flux','psf','residual'):
    os.system('rm -rf test_frequency_linecubesub.{0}'.format(suffix))
