uvcontsub(vis='w51_test_small.ms', fitspw='0:1~20, 0:70~73', solint='int')
os.system('rm -rf test_frequency_contsub.*')
clean(vis='w51_test_small.ms.contsub', imagename="test_frequency_contsub",
      field="", spw='', mode='channel', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_frequency_contsub.image', 'test_frequency_contsub.image.fits', dropdeg=True, overwrite=True)
exportfits('test_frequency_contsub.model', 'test_frequency_contsub.model.fits', dropdeg=True, overwrite=True)

os.system('rm -rf w51_test_small_sub.ms')
split('w51_test_small.ms', 'w51_test_small_sub.ms', datacolumn='data')
clearcal('w51_test_small_sub.ms', addmodel=True)
tb.open('w51_test_small.ms.contsub')
model=tb.getcol("MODEL_DATA")
tb.close()
tb.open('w51_test_small_sub.ms', nomodify=False)
tb.putcol("MODEL_DATA", model)
tb.flush()
tb.close()

uvsub('w51_test_small_sub.ms')

os.system('rm -rf test_frequency_linesub.*')
clean(vis='w51_test_small_sub.ms', imagename="test_frequency_linesub",
      field="", spw='', mode='channel', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_frequency_linesub.image', 'test_frequency_linesub.image.fits', dropdeg=True, overwrite=True)
exportfits('test_frequency_linesub.model', 'test_frequency_linesub.model.fits', dropdeg=True, overwrite=True)

os.system('rm -rf test_mfs_linesub.*')
clean(vis='w51_test_small_sub.ms', imagename="test_mfs_linesub",
      field="", spw='', mode='mfs', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_mfs_linesub.image', 'test_mfs_linesub.image.fits', dropdeg=True, overwrite=True)
exportfits('test_mfs_linesub.model', 'test_mfs_linesub.model.fits', dropdeg=True, overwrite=True)



# SECOND ITERATION
os.system('rm -rf w51_test_small_contsub2.ms')
split('w51_test_small.ms', 'w51_test_small_contsub2.ms', datacolumn='data')
clearcal('w51_test_small_contsub2.ms', addmodel=True)
tb.open('w51_test_small_sub.ms')
model=tb.getcol("MODEL_DATA")
tb.close()
tb.open('w51_test_small_contsub2.ms', nomodify=False)
tb.putcol("MODEL_DATA", model)
tb.flush()
tb.close()

uvsub('w51_test_small_contsub2.ms')
os.system('rm -rf test_frequency_contsub2.*')
clean(vis='w51_test_small_contsub2.ms', imagename="test_frequency_contsub2",
      field="", spw='', mode='channel', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_frequency_contsub2.image', 'test_frequency_contsub2.image.fits', dropdeg=True, overwrite=True)
exportfits('test_frequency_contsub2.model', 'test_frequency_contsub2.model.fits', dropdeg=True, overwrite=True)

os.system('rm -rf w51_test_small_sub2.ms')
split('w51_test_small.ms', 'w51_test_small_sub2.ms', datacolumn='data')
clearcal('w51_test_small_sub2.ms', addmodel=True)
tb.open('w51_test_small_contsub2.ms')
model=tb.getcol("MODEL_DATA")
tb.close()
tb.open('w51_test_small_sub2.ms', nomodify=False)
tb.putcol("MODEL_DATA", model)
tb.flush()
tb.close()

uvsub('w51_test_small_sub2.ms')

os.system('rm -rf test_mfs_linesub2.*')
clean(vis='w51_test_small_sub2.ms', imagename="test_mfs_linesub2",
      field="", spw='', mode='mfs', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_mfs_linesub2.image', 'test_mfs_linesub2.image.fits', dropdeg=True, overwrite=True)
exportfits('test_mfs_linesub2.model', 'test_mfs_linesub2.model.fits', dropdeg=True, overwrite=True)


# THIRD ITERATION
os.system('rm -rf w51_test_small_contsub3.ms')
split('w51_test_small.ms', 'w51_test_small_contsub3.ms', datacolumn='data')
clearcal('w51_test_small_contsub3.ms', addmodel=True)
tb.open('w51_test_small_sub2.ms')
model=tb.getcol("MODEL_DATA")
tb.close()
tb.open('w51_test_small_contsub3.ms', nomodify=False)
tb.putcol("MODEL_DATA", model)
tb.flush()
tb.close()

uvsub('w51_test_small_contsub3.ms')
os.system('rm -rf test_frequency_contsub3.*')
clean(vis='w51_test_small_contsub3.ms', imagename="test_frequency_contsub3",
      field="", spw='', mode='channel', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_frequency_contsub3.image', 'test_frequency_contsub3.image.fits', dropdeg=True, overwrite=True)
exportfits('test_frequency_contsub3.model', 'test_frequency_contsub3.model.fits', dropdeg=True, overwrite=True)

os.system('rm -rf w51_test_small_sub3.ms')
split('w51_test_small.ms', 'w51_test_small_sub3.ms', datacolumn='data')
clearcal('w51_test_small_sub3.ms', addmodel=True)
tb.open('w51_test_small_contsub3.ms')
model=tb.getcol("MODEL_DATA")
tb.close()
tb.open('w51_test_small_sub3.ms', nomodify=False)
tb.putcol("MODEL_DATA", model)
tb.flush()
tb.close()

uvsub('w51_test_small_sub3.ms')

os.system('rm -rf test_mfs_linesub3.*')
clean(vis='w51_test_small_sub3.ms', imagename="test_mfs_linesub3",
      field="", spw='', mode='mfs', outframe='LSRK',
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold='50.0mJy', imsize=[512,512], cell='0.052arcsec',
      weighting='briggs', phasecenter='J2000 19h23m43.905 +14d30m28.08',
      pbcor=False, usescratch=True, robust=1.0)
exportfits('test_mfs_linesub3.image', 'test_mfs_linesub3.image.fits', dropdeg=True, overwrite=True)
exportfits('test_mfs_linesub3.model', 'test_mfs_linesub3.model.fits', dropdeg=True, overwrite=True)
