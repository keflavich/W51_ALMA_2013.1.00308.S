import time
t0 = time.time()

os.system('rm -rf w51_test_small.ms')
os.system('rm -rf w51_test_small.ms.flagversions')
assert split(vis='w51_spw3_continuum_flagged.split',
      outputvis='w51_test_small.ms',
      #field=','.join([str(x-4) for x in (31,32,33,39,40,24,25)]),
      field='28', # 32-4
      spw='',
      datacolumn='data',
     )

print("Done splitting")

solint = 'int'
threshold = '20.0mJy'
multiscale = [0,5,15,45]
#multiscale = []

clearcal(vis='w51_test_small.ms')
#flagmanager(vis='w51_test_small.ms', versionname='flagdata_1', mode='restore')
myimagebase = "test_mfs_dirty"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small.ms', imagename=myimagebase, field="", spw='',
      mode='mfs', outframe='LSRK', interpolation='linear', imagermode='mosaic',
      interactive=False, niter=0, threshold=threshold, imsize=[512,512],
      cell='0.06arcsec', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      weighting='briggs', usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = "test_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small.ms', imagename=myimagebase, field="", spw='',
      mode='mfs', outframe='LSRK', interpolation='linear', imagermode='mosaic',
      multiscale=multiscale,
      interactive=False, niter=10000, threshold=threshold, imsize=[512,512],
      cell='0.06arcsec', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      weighting='briggs', usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables('phase.cal')
gaincal(vis='w51_test_small.ms', caltable="phase.cal", field="", solint='inf',
        calmode="p", refant="", gaintype="G")

#plotcal(caltable="phase.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis='w51_test_small.ms', mode='save', versionname='backup')
applycal(vis="w51_test_small.ms", field="", gaintable=["phase.cal"],
         interp="linear")
flagmanager(vis='w51_test_small.ms', mode='restore', versionname='backup')
os.system('rm -rf w51_test_small_selfcal.ms')
os.system('rm -rf w51_test_small_selfcal.ms.flagversions')
split(vis="w51_test_small.ms", outputvis="w51_test_small_selfcal.ms",
      datacolumn="corrected")

myimagebase="test_selfcal_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small_selfcal.ms', imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables('phase_2.cal')
gaincal(vis="w51_test_small_selfcal.ms", caltable="phase_2.cal", field="",
        solint=solint, calmode="p", refant="", gaintype="G")
#plotcal(caltable="phase_2.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)


flagmanager(vis='w51_test_small_selfcal.ms', mode='save', versionname='backup')
applycal(vis="w51_test_small_selfcal.ms", field="", gaintable=["phase_2.cal"],
         interp="linear")
flagmanager(vis='w51_test_small_selfcal.ms', mode='restore', versionname='backup')
os.system('rm -rf w51_test_small_selfcal_2.ms')
os.system('rm -rf w51_test_small_selfcal_2.ms.flagversions')
split(vis="w51_test_small_selfcal.ms", outputvis="w51_test_small_selfcal_2.ms",
      datacolumn="corrected")

myimagebase = "test_selfcal_2_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small_selfcal_2.ms', imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables("phase_3.cal")
gaincal(vis="w51_test_small_selfcal_2.ms", caltable="phase_3.cal", field="",
        solint=solint, calmode="p", refant="", gaintype="G")
#plotcal(caltable="phase_3.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis='w51_test_small_selfcal_2.ms', mode='save', versionname='backup')
applycal(vis="w51_test_small_selfcal_2.ms", field="", gaintable=["phase_2.cal"],
         interp="linear")
flagmanager(vis='w51_test_small_selfcal_2.ms', mode='restore', versionname='backup')
os.system('rm -rf w51_test_small_selfcal_3.ms')
os.system('rm -rf w51_test_small_selfcal_3.ms.flagversions')
split(vis="w51_test_small_selfcal_2.ms", outputvis="w51_test_small_selfcal_3.ms",
      datacolumn="corrected")

myimagebase = "test_selfcal_3_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small_selfcal_3.ms', imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables("phase_4.cal")
gaincal(vis="w51_test_small_selfcal_3.ms", caltable="phase_4.cal", field="",
        solint=solint, calmode="p", refant="", gaintype="G")
#plotcal(caltable="phase_4.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)


rmtables("ampphase.cal")
gaincal(vis="w51_test_small_selfcal_3.ms", caltable="ampphase.cal", field="",
        solint=solint, solnorm=True, calmode="ap", refant="", gaintype="G")

flagmanager(vis='w51_test_small_selfcal_3.ms', mode='save', versionname='backup')
applycal(vis="w51_test_small_selfcal_3.ms", field="", gaintable=["phase_3.cal"],
         interp="linear")
flagmanager(vis='w51_test_small_selfcal_3.ms', mode='restore', versionname='backup')
os.system('rm -rf w51_test_small_selfcal_4.ms')
os.system('rm -rf w51_test_small_selfcal_4.ms.flagversions')
split(vis="w51_test_small_selfcal_3.ms", outputvis="w51_test_small_selfcal_4.ms",
      datacolumn="corrected")

myimagebase = "test_selfcal_4ampphase_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small_selfcal_4.ms', imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)




os.system('rm -rf w51_test_small_multifield.ms')
os.system('rm -rf w51_test_small_multifield.ms.flagversions')
assert split(vis='w51_spw3_continuum_flagged.split',
      outputvis='w51_test_small_multifield.ms',
      field=','.join([str(x-4) for x in (31,32,33,39,40,24,25)]),
      #field='28', # 32-4
      spw='',
      datacolumn='data',
     )

ft(vis="w51_test_small_multifield.ms", model="test_selfcal_4ampphase_mfs.model",)
rmtables("phase_multifield.cal")
gaincal(vis="w51_test_small_multifield.ms", caltable="phase_multifield.cal",
        field="", solint='inf', solnorm=True, calmode="p", refant="",
        gaintype="G")



myimagebase = "test_multifield_mfs_dirty"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small_multifield.ms', imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=0, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


myimagebase = "test_multifield_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small_multifield.ms', imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


flagmanager(vis='w51_test_small_multifield.ms', mode='save', versionname='backup')
applycal(vis="w51_test_small_multifield.ms", field="",
         gaintable=["phase_multifield.cal"], interp="linear")
flagmanager(vis='w51_test_small_multifield.ms', mode='restore', versionname='backup')



myimagebase = "test_multifield_selfcal_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis='w51_test_small_multifield.ms', imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

import numpy as np
from astropy.io import fits
print("Stats (mfs):")
slc = slice(80,200), slice(80,200)
sigma, peak = (fits.getdata('test_mfs_dirty.image.fits')[slc].std(),     np.nanmax(fits.getdata('test_mfs_dirty.image.fits')))
print("dirty:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_mfs.image.pbcor.fits')[slc].std(),           np.nanmax(fits.getdata('test_mfs.image.pbcor.fits')))
print("clean:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_selfcal_mfs.image.pbcor.fits')[slc].std(),   np.nanmax(fits.getdata('test_selfcal_mfs.image.pbcor.fits')))
print("selfcal:           peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_selfcal_2_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('test_selfcal_2_mfs.image.pbcor.fits')))
print("selfcal2:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_selfcal_3_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('test_selfcal_3_mfs.image.pbcor.fits')))
print("selfcal3:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_selfcal_4ampphase_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('test_selfcal_4ampphase_mfs.image.pbcor.fits')))
print("selfcal4 ampphase: peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_multifield_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('test_multifield_mfs.image.pbcor.fits')))
print("multifield:        peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_multifield_selfcal_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('test_multifield_selfcal_mfs.image.pbcor.fits')))
print("multifield_selfcal peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
print("Completed in {0}s".format(time.time()-t0))
