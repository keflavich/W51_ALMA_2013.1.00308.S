
imexts = ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf', '.residual',
          '.flux.pbcoverage']

fieldoffset = -4
fieldlist = (24,25,31,32,33,39,40,)
fieldmap = {x:x+fieldoffset for x in fieldlist}
inputvis = 'w51_spw3_continuum_flagged.split'
fieldoffset = 0
fieldmap = {x:ii for ii,x in enumerate(sorted(fieldlist))}
inputvis = "w51_test_small_multifield.ms"

for field in fieldlist:

    assert os.path.exists(inputvis)

    myimagebase = "adjacent_field_test_{0}_mosaic".format(field)
    for ext in imexts:
        if os.path.exists(myimagebase+ext):
            rmtables(myimagebase+ext)

    clean(vis=inputvis,
          imagename=myimagebase,
          field=str(fieldmap[field]), spw='', mode='mfs', outframe='LSRK',
          interpolation='linear', imagermode='mosaic',
          interactive=False, niter=10000, threshold='50mJy', imsize=[768,768],
          cell='0.06arcsec',
          weighting='briggs', usescratch=True, pbcor=False, robust=-2.0)

    impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
            outfile=myimagebase+'.image.pbcor', overwrite=True)

    exportfits(myimagebase+".model", myimagebase+".model.fits", dropdeg=True,
               overwrite=True)
    exportfits(myimagebase+".image", myimagebase+".image.fits", dropdeg=True,
               overwrite=True)
    exportfits(myimagebase+".image.pbcor",
               myimagebase+".image.pbcor.fits", dropdeg=True,
               overwrite=True)

    for ext in imexts:
        if os.path.exists(myimagebase+ext):
            rmtables(myimagebase+ext)


    myimagebase = "adjacent_field_test_{0}_csclean".format(field)
    for ext in imexts:
        if os.path.exists(myimagebase+ext):
            rmtables(myimagebase+ext)
    clean(vis=inputvis,
          imagename=myimagebase,
          field=str(fieldmap[field]), spw='', mode='mfs', outframe='LSRK',
          interpolation='linear', imagermode='csclean',
          interactive=False, niter=10000, threshold='50mJy', imsize=[768,768],
          cell='0.06arcsec',
          weighting='briggs', usescratch=True, pbcor=False, robust=-2.0)

    impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
            outfile=myimagebase+'.image.pbcor', overwrite=True)

    exportfits(myimagebase+".model", myimagebase+".model.fits", dropdeg=True,
               overwrite=True)
    exportfits(myimagebase+".image", myimagebase+".image.fits", dropdeg=True,
               overwrite=True)
    exportfits(myimagebase+".image.pbcor",
               myimagebase+".image.pbcor.fits", dropdeg=True,
               overwrite=True)

    for ext in imexts:
        rmtables(myimagebase+ext)
