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

import spectral_cube
cube = SpectralCube.read('test_frequency.image.fits')
cont = cube.min(axis=0)
cont.hdu.write('test_continuum_min.fits')
importfits('test_continuum_min.fits', 'test_continuum_min.image', overwrite=True)

os.system('rm -rf w51_test_small_imcont.ms')
split('w51_test_small.ms', 'w51_test_small_imcont.ms', datacolumn='data')
ft(vis='w51_test_small_imcont.ms', model='test_continuum_min.image', usescratch=True, nterms=0)
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
