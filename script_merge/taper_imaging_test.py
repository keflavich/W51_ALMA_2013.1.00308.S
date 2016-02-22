import os

vis0 = 'w51_continuum_7m12m_contvis_selfcal_0.ms'
multiscale = [0,5,15,45,135]
imsize = [3072,3072]
cell = '0.05arcsec'
phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"

myimagebase = "w51_continuum_7m12m_contvis_taper_lt50m"
# threshold = 50mJy with no other restrictions -> infinite divergence
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       scales=multiscale, interactive=False, niter=10000,
       threshold='100mJy', imsize=imsize, specmode='mfs',
       pblimit=0.5, cell=cell, phasecenter=phasecenter, weighting='briggs',
       robust=-2.0, uvrange='0~50m',
      )
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)


multiscale = [0,5,15]
myimagebase = "w51_continuum_7m12m_contvis_taper_50to500m"
# threshold = 50mJy with no other restrictions -> infinite divergence
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       scales=multiscale, interactive=False, niter=10000,
       threshold='100mJy', imsize=imsize, specmode='mfs',
       pblimit=0.5, cell=cell, phasecenter=phasecenter, weighting='briggs',
       robust=-2.0, uvrange='50~500m',
      )
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)



myimagebase = "w51_continuum_7m12m_contvis_taper_gt500m"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       interactive=False, niter=100000,
       threshold='10mJy', imsize=imsize, specmode='mfs',
       pblimit=0.5, cell=cell, phasecenter=phasecenter, weighting='briggs',
       robust=-2.0, uvrange='500~5000m',
      )
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)



myimagebase = "w51_continuum_7m12m_contvis_taper_gt300m"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       interactive=False, niter=100000,
       threshold='5mJy', imsize=imsize, specmode='mfs',
       pblimit=0.5, cell=cell, phasecenter=phasecenter, weighting='briggs',
       robust=-2.0, uvrange='300~5000m',
       savemodel='modelcolumn',
      )
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)

