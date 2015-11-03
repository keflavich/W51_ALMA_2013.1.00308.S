import os

os.system('rm -rf test_mfs.*')
clean(vis = ['w51_test_small_spw3.ms','w51_test_small_spw7.ms'],
  imagename = "test_mfs",
  field = "",
  spw = '',
  mode = 'mfs',
  outframe = 'LSRK',
  interpolation = 'linear',
  imagermode='mosaic',
  interactive = False,
  niter = 10000,
  threshold = '20.0mJy', # got to 11 mJy, but it showed really bad instabilities.  15 was still too deep
  imsize = [512,512],
  cell = '0.052arcsec',
  phasecenter='J2000 19h23m43.905 +14d30m28.08',
  weighting = 'briggs',
  pbcor=False,
  robust = 0.0)
exportfits('test_mfs.image', 'test_mfs.image.fits')

os.system('rm -rf test_frequency.*')
clean(vis = ['w51_test_small_spw3.ms','w51_test_small_spw7.ms'],
  imagename = "test_frequency",
  field = "",
  spw = '',
  mode = 'channel',
  outframe = 'LSRK',
  interpolation = 'linear',
  imagermode='mosaic',
  interactive = False,
  niter = 10000,
  threshold = '20.0mJy', # got to 11 mJy, but it showed really bad instabilities.  15 was still too deep
  imsize = [512,512],
  cell = '0.052arcsec',
  weighting = 'briggs',
  phasecenter='J2000 19h23m43.905 +14d30m28.08',
  pbcor=False,
  robust = 0.0)
exportfits('test_frequency.image', 'test_frequency.image.fits')
