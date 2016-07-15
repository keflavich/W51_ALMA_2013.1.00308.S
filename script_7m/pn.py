velocity_range = 0,120
velocity_res = 1.0
vis='calibrated_7m.ms'
line, restfreq, velocity_res, spw = ('PN5-4', "234.93569GHz",1.2, 3)
nchans = int((velocity_range[1]-velocity_range[0])/velocity_res)

for spw in (3,7,11,15):
    output='W51_B6_7m.{0}.spw{1}'.format(line,spw)
    tclean(vis=vis,
           imagename=output+".tclean",
           field='w51',
           spw='{0}'.format(spw),
           gridder='mosaic',
           specmode='cube',
           start='{0}km/s'.format(velocity_range[0]),
           width='{0}km/s'.format(velocity_res), interpolation='linear',
           nchan=nchans,
           restfreq=restfreq,
           veltype='radio',
           outframe='LSRK',
           interactive=F,
           niter=1000,
           imsize=[512,512],
           cell='0.5arcsec',
           weighting='briggs',
           robust=0.5,
           phasecenter='',
           threshold='15mJy',
           savemodel='none',
          )
    myimagebase = output+".tclean"
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True, dropdeg=True)
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', overwrite=True, dropdeg=True)
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', overwrite=True, dropdeg=True)

    clean(vis=vis,
           imagename=output+".clean",
           field='w51',
           spw='',
           imagermode='mosaic',
           mode='velocity',
           start='{0}km/s'.format(velocity_range[0]),
           width='{0}km/s'.format(velocity_res), interpolation='linear',
           nchan=nchans,
           restfreq=restfreq,
           veltype='radio',
           outframe='LSRK',
           interactive=F,
           niter=1000,
           imsize=[512,512],
           cell='0.5arcsec',
           weighting='briggs',
           robust=0.5,
           phasecenter='',
           threshold='15mJy',
          )
    myimagebase = output+".clean"
    exportfits(imagename=myimagebase+'.image', fitsimage=myimagebase+'.image.fits', overwrite=True, dropdeg=True)

# full spw0
tclean(vis=vis,
       imagename='full_spw0_7m'+".tclean",
       field='w51',
       spw='0,4,8,12',
       gridder='mosaic',
       specmode='cube',
       veltype='radio',
       outframe='LSRK',
       interactive=F,
       niter=100,
       imsize=[256,256],
       cell='0.5arcsec',
       weighting='briggs',
       robust=0.5,
       phasecenter='',
       threshold='15mJy',
       savemodel='none',
      )
myimagebase = 'full_spw0_7m.tclean'
exportfits(imagename=myimagebase+'.image', fitsimage=myimagebase+'.image.fits', overwrite=True, dropdeg=True)
