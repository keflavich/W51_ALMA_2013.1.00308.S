
imexts = ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual',
          '.flux.pbcoverage']

for field in [(x-4) for x in (31,32,33,39,40,24,25)]:

    myimagebase = "adjacent_field_test_{0}_mosaic".format(field)
    clean(vis='w51_spw3_continuum_flagged.split',
          imagename=myimagebase,
          field=str(field), spw='', mode='mfs', outframe='LSRK',
          interpolation='linear', imagermode='mosaic',
          interactive=False, niter=10000, threshold='50mJy', imsize=[512,512],
          cell='0.06arcsec',
          weighting='briggs', usescratch=True, pbcor=False, robust=-2.0)

    impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
            outfile=myimagebase+'.image.pbcor', overwrite=True)

    exportfits(myimagebase+".image", myimagebase+".image.fits", dropdeg=True,
               overwrite=True)
    exportfits(myimagebase+".image.pbcor",
               myimagebase+".image.pbcor.fits", dropdeg=True,
               overwrite=True)

    for ext in imexts:
        rmtables(myimagebase+ext)


    myimagebase = "adjacent_field_test_{0}_csclean".format(field)
    clean(vis='w51_spw3_continuum_flagged.split',
          imagename=myimagebase,
          field=str(field), spw='', mode='mfs', outframe='LSRK',
          interpolation='linear', imagermode='csclean',
          interactive=False, niter=10000, threshold='50mJy', imsize=[512,512],
          cell='0.06arcsec',
          weighting='briggs', usescratch=True, pbcor=False, robust=-2.0)

    impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
            outfile=myimagebase+'.image.pbcor', overwrite=True)

    exportfits(myimagebase+".image", myimagebase+".image.fits", dropdeg=True,
               overwrite=True)
    exportfits(myimagebase+".image.pbcor",
               myimagebase+".image.pbcor.fits", dropdeg=True,
               overwrite=True)

    for ext in imexts:
        rmtables(myimagebase+ext)
