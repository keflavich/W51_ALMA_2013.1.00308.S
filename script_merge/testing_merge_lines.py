
linesub='w51_concat_7m12m.spw0.merge.contsub' #result of uvcontsub contains only the selected uvdata by the task

rmtables(['h2co_line_test.ms'])
split(vis=linesub,
      outputvis='h2co_line_test.ms',
      field="w51",
      spw='0:3006~3036',
      width=5,
      datacolumn='data',)
listobs(vis=linesub, listfile=linesub+".listfile", overwrite=True)
listobs(vis='h2co_line_test.ms', listfile='h2co_line_test.ms'+".listfile",
        overwrite=True)

myimagebase="h2co_line_test_tclean"
rmtables(["h2co_line_test.{0}".format(x)
          for x in ('psf', 'weight', 'sumwt', 'model', 'residual', 'image',
                    'image.pbcor')])
tclean(vis='h2co_line_test.ms',
       imagename = myimagebase,
       field = "w51",
       spw = '',
       specmode = 'cube',
       restfreq = '218.222GHz',
       outframe = 'LSRK',
       interpolation = 'linear',
       gridder='mosaic',
       interactive = False,
       niter = 0,
       threshold = '10mJy', #req rms 5.85 mJy, 38 antennas, 33.8min tos, pwv 2 (for EB Xb4b),0.2 arcsec res, 0.728 MHz BW gives 1.8mJy!
       imsize = [960,960], # make it much bigger to avoid edge effects
       cell = '0.15arcsec',  #synth beam expected to be 0.2 arcsec, so 0.2/3= 0.06 arcsec cell
       weighting = 'briggs',
       #mask='auto-pb',
       pblimit=0.5,
       )
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits',
           overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.image',fitsimage=myimagebase+'.image.fits',
           overwrite=True) # export the PB image


myimagebase="h2co_line_test_clean"
rmtables(["h2co_line_test.{0}".format(x)
          for x in ('psf', 'weight', 'flux', 'model', 'residual', 'image',
                    'image.pbcor')])
clean(vis='h2co_line_test.ms',
      imagename = myimagebase,
      field = "w51",
      spw = '',
      mode = 'velocity',
      restfreq = '218.222GHz',
      outframe = 'LSRK',
      interpolation = 'linear',
      imagermode='mosaic',
      interactive = False,
      niter = 0,
      threshold = '10mJy', #req rms 5.85 mJy, 38 antennas, 33.8min tos, pwv 2 (for EB Xb4b),0.2 arcsec res, 0.728 MHz BW gives 1.8mJy!
      imsize = [960,960], # make it much bigger to avoid edge effects
      cell = '0.15arcsec',  #synth beam expected to be 0.2 arcsec, so 0.2/3= 0.06 arcsec cell
      weighting = 'briggs',
      minpb=0.5,
      )
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits',
           overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.image',fitsimage=myimagebase+'.image.fits',
           overwrite=True) # export the PB image
