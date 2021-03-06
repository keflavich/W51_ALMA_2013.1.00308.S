if not os.path.exists('w51_progressive_test_small.ms'):
    os.system('cp -r w51_test_small.ms w51_progressive_test_small.ms')
#    split(vis='w51_spw3_continuum_flagged.ms',
#          outputvis='w51_progressive_test_small.ms',
#          field='31,32,33,39,40,24,25',
#          spw='',
#          datacolumn='data',
#         )

assert os.path.exists('w51_progressive_test_small.ms')

clearcal(vis='w51_progressive_test_small.ms')
os.system('rm -rf progressive_test_mfs_dirty.*')
flagmanager(vis='w51_progressive_test_small.ms', versionname='flagdata_1', mode='restore')
clean(vis='w51_progressive_test_small.ms', imagename="progressive_test_mfs_dirty", field="", spw='',
      mode='mfs', outframe='LSRK', interpolation='linear', imagermode='mosaic',
      interactive=False, niter=0, imsize=[512,512],
      cell='0.06arcsec', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      weighting='briggs', usescratch=True, pbcor=False, robust=-2.0)
exportfits('progressive_test_mfs_dirty.image', 'progressive_test_mfs_dirty.image.fits', dropdeg=True, overwrite=True)
os.system('rm -rf progressive_test_mfs.*')
clean(vis='w51_progressive_test_small.ms', imagename="progressive_test_mfs", field="", spw='',
      mode='mfs', outframe='LSRK', interpolation='linear', imagermode='mosaic',
      interactive=False, niter=1000, threshold='50.0mJy', imsize=[512,512],
      cell='0.06arcsec', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      weighting='briggs', usescratch=True, pbcor=False, robust=-2.0)
exportfits('progressive_test_mfs.image', 'progressive_test_mfs.image.fits', dropdeg=True, overwrite=True)

gaincal(vis='w51_progressive_test_small.ms', caltable="phase.cal", field="", solint="30s",
        calmode="p", refant="", gaintype="G")
#plotcal(caltable="phase.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis='w51_progressive_test_small.ms', mode='save', versionname='backup')
applycal(vis="w51_progressive_test_small.ms", field="", gaintable=["phase.cal"],
         interp="linear", uvrange='400~2000', applymode='calonly')
flagmanager(vis='w51_progressive_test_small.ms', mode='restore', versionname='backup')
os.system('rm -rf w51_progressive_test_small_selfcal.ms')
os.system('rm -rf w51_progressive_test_small_selfcal.ms.flagversions')
split(vis="w51_progressive_test_small.ms", outputvis="w51_progressive_test_small_selfcal.ms",
      datacolumn="corrected")

os.system('rm -rf progressive_test_selfcal_mfs.*')
clean(vis='w51_progressive_test_small_selfcal.ms', imagename="progressive_test_selfcal_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=1000, threshold='50.0mJy', imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('progressive_test_selfcal_mfs.image', 'progressive_test_selfcal_mfs.image.fits', dropdeg=True, overwrite=True)

os.system("rm -rf phase_2.cal")
gaincal(vis="w51_progressive_test_small_selfcal.ms", caltable="phase_2.cal", field="",
        solint="30s", calmode="p", refant="", gaintype="G")
#plotcal(caltable="phase_2.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)


flagmanager(vis='w51_progressive_test_small_selfcal.ms', mode='save', versionname='backup')
applycal(vis="w51_progressive_test_small_selfcal.ms", field="", gaintable=["phase_2.cal"],
         interp="linear", uvrange='200~2000', applymode='calonly')
flagmanager(vis='w51_progressive_test_small_selfcal.ms', mode='restore', versionname='backup')
os.system('rm -rf w51_progressive_test_small_selfcal_2.ms')
os.system('rm -rf w51_progressive_test_small_selfcal_2.ms.flagversions')
split(vis="w51_progressive_test_small_selfcal.ms", outputvis="w51_progressive_test_small_selfcal_2.ms",
      datacolumn="corrected")

os.system('rm -rf progressive_test_selfcal_2_mfs.*')
clean(vis='w51_progressive_test_small_selfcal_2.ms', imagename="progressive_test_selfcal_2_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=1000, threshold='50.0mJy', imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('progressive_test_selfcal_2_mfs.image', 'progressive_test_selfcal_2_mfs.image.fits', dropdeg=True, overwrite=True)

os.system("rm -rf phase_3.cal")
gaincal(vis="w51_progressive_test_small_selfcal_2.ms", caltable="phase_3.cal", field="",
        solint="30s", calmode="p", refant="", gaintype="G")
#plotcal(caltable="phase_3.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis='w51_progressive_test_small_selfcal_2.ms', mode='save',
            versionname='backup')
applycal(vis="w51_progressive_test_small_selfcal_2.ms", field="",
         gaintable=["phase_3.cal"],
         interp="linear", applymode='calonly')
flagmanager(vis='w51_progressive_test_small_selfcal.ms', mode='restore',
            versionname='backup')
os.system('rm -rf w51_progressive_test_small_selfcal_3.ms')
os.system('rm -rf w51_progressive_test_small_selfcal_3.ms.flagversions')
split(vis="w51_progressive_test_small_selfcal_2.ms",
      outputvis="w51_progressive_test_small_selfcal_3.ms",
      datacolumn="corrected")

os.system('rm -rf progressive_test_selfcal_3_mfs.*')
clean(vis='w51_progressive_test_small_selfcal_3.ms',
      imagename="progressive_test_selfcal_3_mfs",
      field="", spw='', mode='mfs', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=1000, threshold='50.0mJy', imsize=[512,512], cell='0.06arcsec',
      phasecenter='J2000 19h23m43.905 +14d30m28.08', weighting='briggs',
      usescratch=True, pbcor=False, robust=-2.0)
exportfits('progressive_test_selfcal_3_mfs.image',
           'progressive_test_selfcal_3_mfs.image.fits', dropdeg=True,
           overwrite=True)


from astropy.io import fits
print("Stats (mfs):")
print("dirty:    peak={1:0.5f} sigma={0:0.5f}".format(fits.getdata('progressive_test_mfs_dirty.image.fits')[:200,:200].std(),     fits.getdata('progressive_test_mfs_dirty.image.fits').max()))
print("clean:    peak={1:0.5f} sigma={0:0.5f}".format(fits.getdata('progressive_test_mfs.image.fits')[:200,:200].std(),           fits.getdata('progressive_test_mfs.image.fits').max()))
print("selfcal:  peak={1:0.5f} sigma={0:0.5f}".format(fits.getdata('progressive_test_selfcal_mfs.image.fits')[:200,:200].std(),   fits.getdata('progressive_test_selfcal_mfs.image.fits').max()))
print("selfcal2: peak={1:0.5f} sigma={0:0.5f}".format(fits.getdata('progressive_test_selfcal_2_mfs.image.fits')[:200,:200].std(), fits.getdata('progressive_test_selfcal_2_mfs.image.fits').max()))
print("selfcal3: peak={1:0.5f} sigma={0:0.5f}".format(fits.getdata('progressive_test_selfcal_3_mfs.image.fits')[:200,:200].std(), fits.getdata('progressive_test_selfcal_3_mfs.image.fits').max()))
