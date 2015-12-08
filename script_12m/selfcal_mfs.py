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
os.system('rm -rf test_mfs_dirty.*')
#flagmanager(vis='w51_test_small.ms', versionname='flagdata_1', mode='restore')
clean(vis='w51_test_small.ms', imagename="test_mfs_dirty", field="", spw='',
      mode='mfs', outframe='LSRK', interpolation='linear', imagermode='mosaic',
      interactive=False, niter=0, threshold=threshold, imsize=[512,512],
      cell='0.06arcsec', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      weighting='briggs', usescratch=True, pbcor=False, robust=-2.0)
exportfits('test_mfs_dirty.image', 'test_mfs_dirty.image.fits', dropdeg=True, overwrite=True)
exportfits('test_mfs_dirty.model', 'test_mfs_dirty.model.fits', dropdeg=True, overwrite=True)
exportfits('test_mfs_dirty.residual', 'test_mfs_dirty.residual.fits', dropdeg=True, overwrite=True)

os.system('rm -rf test_mfs.*')
clean(vis='w51_test_small.ms', imagename="test_mfs", field="", spw='',
      mode='mfs', outframe='LSRK', interpolation='linear', imagermode='mosaic',
      multiscale=multiscale,
      interactive=False, niter=10000, threshold=threshold, imsize=[512,512],
      cell='0.06arcsec', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      weighting='briggs', usescratch=True, pbcor=False, robust=-2.0)
exportfits('test_mfs.image', 'test_mfs.image.fits', dropdeg=True, overwrite=True)
exportfits('test_mfs.model', 'test_mfs.model.fits', dropdeg=True, overwrite=True)
exportfits('test_mfs.residual', 'test_mfs.residual.fits', dropdeg=True, overwrite=True)

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

os.system('rm -rf test_selfcal_mfs.*')
clean(vis='w51_test_small_selfcal.ms', imagename="test_selfcal_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('test_selfcal_mfs.image', 'test_selfcal_mfs.image.fits', dropdeg=True, overwrite=True)
exportfits('test_selfcal_mfs.model', 'test_selfcal_mfs.model.fits', dropdeg=True, overwrite=True)
exportfits('test_selfcal_mfs.residual', 'test_selfcal_mfs.residual.fits', dropdeg=True, overwrite=True)

os.system("rm -rf phase_2.cal")
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

os.system('rm -rf test_selfcal_2_mfs.*')
clean(vis='w51_test_small_selfcal_2.ms', imagename="test_selfcal_2_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('test_selfcal_2_mfs.image', 'test_selfcal_2_mfs.image.fits', dropdeg=True, overwrite=True)
exportfits('test_selfcal_2_mfs.model', 'test_selfcal_2_mfs.model.fits', dropdeg=True, overwrite=True)
exportfits('test_selfcal_2_mfs.residual', 'test_selfcal_2_mfs.residual.fits', dropdeg=True, overwrite=True)

os.system("rm -rf phase_3.cal")
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

os.system('rm -rf test_selfcal_3_mfs.*')
clean(vis='w51_test_small_selfcal_3.ms', imagename="test_selfcal_3_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('test_selfcal_3_mfs.image', 'test_selfcal_3_mfs.image.fits', dropdeg=True, overwrite=True)
exportfits('test_selfcal_3_mfs.model', 'test_selfcal_3_mfs.model.fits', dropdeg=True, overwrite=True)
exportfits('test_selfcal_3_mfs.residual', 'test_selfcal_3_mfs.residual.fits', dropdeg=True, overwrite=True)

os.system("rm -rf phase_4.cal")
gaincal(vis="w51_test_small_selfcal_3.ms", caltable="phase_4.cal", field="",
        solint=solint, calmode="p", refant="", gaintype="G")
#plotcal(caltable="phase_4.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)


os.system("rm -rf ampphase.cal")
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

os.system('rm -rf test_selfcal_4_mfs.*')
clean(vis='w51_test_small_selfcal_4.ms', imagename="test_selfcal_4ampphase_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('test_selfcal_4ampphase_mfs.image', 'test_selfcal_4ampphase_mfs.image.fits', dropdeg=True, overwrite=True)
exportfits('test_selfcal_4ampphase_mfs.model', 'test_selfcal_4ampphase_mfs.model.fits', dropdeg=True, overwrite=True)
exportfits('test_selfcal_4ampphase_mfs.residual', 'test_selfcal_4ampphase_mfs.residual.fits', dropdeg=True, overwrite=True)



os.system('rm -rf w51_test_small_multifield.ms')
os.system('rm -rf w51_test_small_multifield.ms.flagversions')
assert split(vis='w51_spw3_continuum_flagged.split',
      outputvis='w51_test_small_multifield.ms',
      field=','.join([str(x-4) for x in (31,32,33,39,40,24,25)]),
      #field='28', # 32-4
      spw='',
      datacolumn='data',
     )


os.system('rm -rf test_multifield_mfs.*')
clean(vis='w51_test_small_multifield.ms', imagename="test_multifield_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('test_multifield_mfs.image', 'test_multifield_mfs.image.fits', dropdeg=True, overwrite=True)
exportfits('test_multifield_mfs.model', 'test_multifield_mfs.model.fits', dropdeg=True, overwrite=True)
exportfits('test_multifield_mfs.residual', 'test_multifield_mfs.residual.fits', dropdeg=True, overwrite=True)


flagmanager(vis='w51_test_small_multifield.ms', mode='save', versionname='backup')
applycal(vis="w51_test_small_multifield.ms", field="",
         gaintable=["phase.cal","phase_2.cal","phase_3.cal"], interp="linear")
flagmanager(vis='w51_test_small_multifield.ms', mode='restore', versionname='backup')



os.system('rm -rf test_multifield_selfcal_mfs.*')
clean(vis='w51_test_small_multifield.ms', imagename="test_multifield_selfcal_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('test_multifield_selfcal_mfs.image', 'test_multifield_selfcal_mfs.image.fits', dropdeg=True, overwrite=True)
exportfits('test_multifield_selfcal_mfs.model', 'test_multifield_selfcal_mfs.model.fits', dropdeg=True, overwrite=True)
exportfits('test_multifield_selfcal_mfs.residual', 'test_multifield_selfcal_mfs.residual.fits', dropdeg=True, overwrite=True)

import numpy as np
from astropy.io import fits
print("Stats (mfs):")
slc = slice(80,200), slice(80,200)
sigma, peak = (fits.getdata('test_mfs_dirty.image.fits')[slc].std(),     np.nanmax(fits.getdata('test_mfs_dirty.image.fits')))
print("dirty:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_mfs.image.fits')[slc].std(),           np.nanmax(fits.getdata('test_mfs.image.fits')))
print("clean:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_selfcal_mfs.image.fits')[slc].std(),   np.nanmax(fits.getdata('test_selfcal_mfs.image.fits')))
print("selfcal:           peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_selfcal_2_mfs.image.fits')[slc].std(), np.nanmax(fits.getdata('test_selfcal_2_mfs.image.fits')))
print("selfcal2:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_selfcal_3_mfs.image.fits')[slc].std(), np.nanmax(fits.getdata('test_selfcal_3_mfs.image.fits')))
print("selfcal3:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_selfcal_4ampphase_mfs.image.fits')[slc].std(), np.nanmax(fits.getdata('test_selfcal_4ampphase_mfs.image.fits')))
print("selfcal4 ampphase: peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_multifield_mfs.image.fits')[slc].std(), np.nanmax(fits.getdata('test_multifield_mfs.image.fits')))
print("multifield:        peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('test_multifield_selfcal_mfs.image.fits')[slc].std(), np.nanmax(fits.getdata('test_multifield_selfcal_mfs.image.fits')))
print("multifield_selfcal peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
