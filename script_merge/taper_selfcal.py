import os

field=''
solint='int'
vis0 = 'w51_continuum_7m12m_contvis_selfcal_0.ms'
multiscale = [0,5,15,45,135]
imsize = [3072,3072]
cell = '0.05arcsec'
phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"

contvis='continuum_7m12m.ms'
vis0 = 'w51_continuum_7m12m_contvis_selfcal_0.ms'
vis1 = 'w51_continuum_7m12m_contvis_selfcal_1.ms'
vis2 = 'w51_continuum_7m12m_contvis_selfcal_2.ms'

os.system('rm -rf {0}'.format(vis0))
os.system('rm -rf {0}.flagversions'.format(vis0))
assert split(vis=contvis,
      outputvis=vis0,
      #field=','.join([str(x-4) for x in (31,32,33,39,40,24,25)]),
      field='w51', # 32-4
      spw='',
      datacolumn='data',
     )

myimagebase = "w51_continuum_7m12m_contvis_taper_gt300m_selfcaliter0"
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

gaincal(vis=vis0, caltable="phase_longbaselines_0.cal", field=field,
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5)
applycal(vis=vis0, field="", gaintable=["phase_longbaselines_0.cal"],
         interp="linear", applymode='calonly', calwt=False)

split(vis=vis0, outputvis=vis1,
      datacolumn="corrected")

myimagebase = "w51_continuum_7m12m_contvis_taper_gt300m_selfcaliter1"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis1, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       interactive=False, niter=100000,
       threshold='5mJy', imsize=imsize, specmode='mfs',
       pblimit=0.5, cell=cell, phasecenter=phasecenter, weighting='briggs',
       robust=-2.0, uvrange='300~5000m',
       savemodel='modelcolumn',
      )
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)

gaincal(vis=vis1, caltable="phase_longbaselines_1.cal", field=field,
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5)
applycal(vis=vis1, field="", gaintable=["phase_longbaselines_1.cal"],
         interp="linear", applymode='calonly', calwt=False)

split(vis=vis1, outputvis=vis2,
      datacolumn="corrected")

myimagebase = "w51_continuum_7m12m_contvis_taper_gt300m_selfcaliter2"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis2, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       interactive=False, niter=100000,
       threshold='5mJy', imsize=imsize, specmode='mfs',
       pblimit=0.5, cell=cell, phasecenter=phasecenter, weighting='briggs',
       robust=-2.0, uvrange='300~5000m',
       savemodel='modelcolumn',
      )
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
