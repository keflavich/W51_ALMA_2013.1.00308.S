"""
Look for pointing offsets between 7m and different 12m pointings in the UV
overlap range (~15-40m)
"""
vis = 'continuum_7m12m.ms'

field='w51'
multiscale = [0,5,15]
imsize = [256,256]
cell = '1arcsec'
phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"

for obsid in range(6):

    # there's no obvious way to determine whether it is 7m or 12m
    # you can specify obs#, but that doesn't tell you...
    # would be nice to have some non-interactive way to determine that
    myimagebase = 'continuum_7m12m_obs{0}'.format(obsid)
    os.system('rm -rf {0}.*'.format(myimagebase))
    tclean(vis=vis, imagename=myimagebase, field=field, spw='',
           observation=str(obsid),
           outframe='LSRK', interpolation='linear', gridder='mosaic',
           scales=multiscale,
           interactive=False, niter=1000,
           threshold='100mJy', imsize=imsize, specmode='mfs',
           pblimit=0.1, cell=cell, phasecenter=phasecenter, weighting='briggs',
           robust=-2.0, uvrange='15~40m',
           savemodel='modelcolumn',
          )
    exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True,
               overwrite=True)
    exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True,
               overwrite=True)

myimagebase = 'continuum_7m12m_allobs'.format(obsid)
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis, imagename=myimagebase, field=field, spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       scales=multiscale,
       interactive=False, niter=1000,
       threshold='100mJy', imsize=imsize, specmode='mfs',
       pblimit=0.1, cell=cell, phasecenter=phasecenter, weighting='briggs',
       robust=-2.0, uvrange='15~40m',
       savemodel='modelcolumn',
      )
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True,
           overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True,
           overwrite=True)
