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
cube = spectral_cube.SpectralCube.read('test_frequency.image.fits')
cont = cube.min(axis=0)
ppbeam = np.abs((cube.beam.sr / (cube.wcs.pixel_scale_matrix[0,0]*cube.wcs.pixel_scale_matrix[1,1]*u.deg**2)).decompose())
hdu = cont.hdu
hdu.data *= ppbeam # because apparently CASA divides by this?
hdu.writeto('test_continuum_min.fits', clobber=True)
header = cont.header
importfits('test_continuum_min.fits', 'test_continuum_min.image',
           overwrite=True, defaultaxes=T,
           defaultaxesvalues=[header['CRVAL1'], header['CRVAL2'],
                              header['RESTFRQ'], 'I'])

os.system('rm -rf w51_test_small_imcont.ms')
split('w51_test_small.ms', 'w51_test_small_imcont.ms', datacolumn='data')

im.open('w51_test_small_imcont.ms')
from astropy import coordinates
c = coordinates.SkyCoord(header['CRVAL1'], header['CRVAL2'], unit=('deg','deg'), frame='fk5')
im.defineimage(nx=cube.shape[1], cellx='{0}arcsec'.format(header['CDELT2']*3600), phasecenter='J2000 {0} {1}'.format(c.ra.to_string(), c.dec.to_string())),
im.makemodelfromsd(sdimage='test_continuum_min.image', modelimage='model_test')
im.close()
setjy(vis='w51_test_small_imcont.ms', model='model_test', usescratch=True, field='', standard='manual')

#ft(vis='w51_test_small_imcont.ms', model='test_continuum_min.image', usescratch=True, nterms=1)

uvsub('w51_test_small_imcont.ms')

os.system('rm -rf test_frequency_mincontsub.*')
clean(vis='w51_test_small_imcont.ms', imagename="test_frequency_mincontsub",
      field="", spw='', mode='channel', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_frequency_mincontsub.image', 'test_frequency_mincontsub.image.fits', dropdeg=True, overwrite=True)
exportfits('test_frequency_mincontsub.model', 'test_frequency_mincontsub.model.fits', dropdeg=True, overwrite=True)
for suffix in ('image','model','flux','psf','residual'):
    os.system('rm -rf test_frequency_mincontsub.{0}'.format(suffix))
