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
cube = spectral_cube.SpectralCube.read('test_frequency.image.fits')
cont = cube[1:-1,:,:].min(axis=0)
scube = cube - cont
hdu = scube.hdu
hdu.data = hdu.data.reshape((1,) + hdu.data.shape)
hdu.header['CRVAL1'] = cube.wcs.wcs.crval[0]
hdu.header['CRVAL2'] = cube.wcs.wcs.crval[1]
hdu.header['CRPIX4'] = 1
hdu.header['CRVAL4'] = 1
hdu.header['CDELT4'] = 1
hdu.header['CTYPE4'] = 'STOKES'
hdu.writeto('test_frequency_minsub.fits', clobber=True)
importfits('test_frequency_minsub.fits', 'test_frequency_minsub.image', overwrite=True)

os.system('rm -rf w51_test_small_linecubesub.ms')
split('w51_test_small.ms', 'w51_test_small_linecubesub.ms', datacolumn='data')
ft(vis='w51_test_small_linecubesub.ms', model='test_frequency_minsub.image', usescratch=True, nterms=1)
raise ValueError("ft does not work.")
uvsub('w51_test_small_linecubesub.ms')

os.system('rm -rf test_mfs_linecubesub.*')
clean(vis='w51_test_small_linecubesub.ms', imagename="test_mfs_linecubesub",
      field="", spw='', mode='mfs', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_mfs_linecubesub.image', 'test_mfs_linecubesub.image.fits', dropdeg=True, overwrite=True)
exportfits('test_mfs_linecubesub.model', 'test_mfs_linecubesub.model.fits', dropdeg=True, overwrite=True)
for suffix in ('image','model','flux','psf','residual'):
    os.system('rm -rf test_mfs_linecubesub.{0}'.format(suffix))
