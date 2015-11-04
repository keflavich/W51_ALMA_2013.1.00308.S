import os

concat(vis=['w51_test_small_spw3.ms','w51_test_small_spw7.ms'],
       concatvis='w51_test_small.ms')
flagdata(vis='w51_test_small.ms',
         mode='manual', timerange='2015/08/27/01:53:32.3~2015/08/27/01:53:56.5',
         antenna='DV08&DV09')

os.system('rm -rf test_mfs.*')
clean(vis = 'w51_test_small.ms',
  imagename = "test_mfs",
      multiscale=[0,3,6,12],
  field = "",
  spw = '',
  mode = 'mfs',
  outframe = 'LSRK',
  interpolation = 'linear',
  imagermode='mosaic',
  interactive = False,
  niter = 10000,
  threshold = '2.0mJy', # got to 11 mJy, but it showed really bad instabilities.  15 was still too deep
  imsize = [512,512],
  cell = '0.052arcsec',
  phasecenter='J2000 19h23m43.905 +14d30m28.08',
  weighting = 'briggs',
      usescratch=True,
  pbcor=False,
  robust = 0.0)
exportfits('test_mfs.image', 'test_mfs.image.fits', dropdeg=True, overwrite=True)

os.system('rm -rf test_frequency.*')
clean(vis = 'w51_test_small.ms',
  imagename = "test_frequency",
  #    multiscale=[0,3,6,12],
  field = "",
  spw = '',
  mode = 'channel',
  outframe = 'LSRK',
  interpolation = 'linear',
  imagermode='mosaic',
  interactive = False,
  niter = 10000,
  threshold = '2.0mJy', # got to 11 mJy, but it showed really bad instabilities.  15 was still too deep
  imsize = [512,512],
  cell = '0.052arcsec',
  weighting = 'briggs',
  phasecenter='J2000 19h23m43.905 +14d30m28.08',
  pbcor=False,
      usescratch=True,
  robust = 1.0)
exportfits('test_frequency.image', 'test_frequency.image.fits', dropdeg=True, overwrite=True)

import spectral_cube
cube = spectral_cube.SpectralCube.read('test_frequency.image.fits')
mean = cube.mean(axis=0)
mean.hdu.writeto('test_frequency.mean.fits', clobber=True)
median = cube.median(axis=0)
median.hdu.writeto('test_frequency.median.fits', clobber=True)
from astropy.io import fits
mfs = fits.getdata('test_mfs.image.fits')

print("Stats:")
print("MFS:    sigma={0:0.5f}".format(mfs[:200,:200].std()))
print("median: sigma={0:0.5f}".format(median.value[:200,:200].std()))
print("mean:   sigma={0:0.5f}".format(mean.value[:200,:200].std()))
