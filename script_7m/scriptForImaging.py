
# continuum mapping
os.system('rm -rf w51-cont.*')
# Just ignoring all line emission contamination...
clean(vis="calibrated.ms",
      imagename="w51-cont",
      field="w51*",
      spw="1,2,3",
      mode="mfs",
      niter=1000,gain=0.1,threshold="0.0mJy",
      psfmode="clark",imagermode="mosaic",ftmachine="mosaic",
      interactive=True,
      outframe="BARY",veltype="optical",
      imsize=128,cell="1.2arcsec",
      phasecenter="",
      stokes="I",
      weighting="briggs",robust=0.5,
      pbcor=False)


# continuum subtraction 
# (is diffcult because so many lines...)
split('calibrated.ms',spw='0',outputvis='calibrated-spw0.ms')
split('calibrated.ms',spw='1',outputvis='calibrated-spw1.ms')
split('calibrated.ms',spw='2',outputvis='calibrated-spw2.ms')
split('calibrated.ms',spw='3',outputvis='calibrated-spw3.ms')

uvcontsub(
    vis                 = 'calibrated-spw0.ms',  
    field               =     'w51*',       
    fitspw              = '0:90~1000;1680~2300;2980~3180;3470~3950')
uvcontsub(
    vis                 = 'calibrated-spw1.ms',  
    field               =     'w51*',       
    fitspw              = '0:1750~2750;3350~3750')
uvcontsub(
    vis                 = 'calibrated-spw2.ms',  
    field               =     'w51*',       
    fitspw              = '0:460~1260;1930~3180')
uvcontsub(
    vis                 = 'calibrated-spw3.ms',  # too hard - too many lines
    field               =     'w51*',       
    fitspw              = '0:510~770;3270~3370')


# line mapping: H2CO
os.system('rm -rf w51-H2CO.*')
clean(vis="calibrated-spw0.ms.contsub",
      imagename="w51-H2CO",
      field="w51*",
      spw="0",
      mode="velocity",
      niter=1000,gain=0.1,threshold="0.0mJy",
      psfmode="clark",imagermode="mosaic",ftmachine="mosaic",
      interactive=True,
      nchan=50,start="40km/s",
      width="1km/s",outframe="LSRK",veltype="radio",
      imsize=128,cell="1.2arcsec",
      phasecenter="",
      restfreq="218.22219GHz",stokes="I",
      weighting="briggs",robust=0.5,
      pbcor=False)

# line mapping: C18O
os.system('rm -rf w51-C18O.*')
clean(vis="calibrated-spw1.ms.contsub",
      imagename="w51-C18O",
      field="w51*",
      spw="0",
      mode="velocity",
      niter=1000,gain=0.1,threshold="0.0mJy",
      psfmode="clark",imagermode="mosaic",ftmachine="mosaic",
      interactive=True,
      nchan=50,start="40km/s",
      width="1km/s",outframe="LSRK",veltype="radio",
      imsize=128,cell="1.2arcsec",
      phasecenter="",
      restfreq="219.5603582GHz",stokes="I",
      weighting="briggs",robust=0.5,
      pbcor=False)

# line mapping: CO(2-1)
os.system('rm -rf w51-CO.*')
clean(vis="calibrated-spw2.ms.contsub",
      imagename="w51-CO",
      field="w51*",
      spw="0",
      mode="velocity",
      niter=1000,gain=0.1,threshold="0.0mJy",
      psfmode="clark",imagermode="mosaic",ftmachine="mosaic",
      interactive=True,
      nchan=50,start="40km/s",
      width="1km/s",outframe="LSRK",veltype="radio",
      imsize=128,cell="1.2arcsec",
      phasecenter="",
      restfreq="230.53800GHz",stokes="I",
      weighting="briggs",robust=0.5,
      pbcor=False)




myimagebase = 'w51-cont'
impbcor(imagename=myimagebase+'.image',
        pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor',
        overwrite=True) # perform PBcorr

exportfits(imagename=myimagebase+'.image.pbcor',
           fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',
           fitsimage=myimagebase+'.flux.fits') # export the corrected image

myimagebase = 'w51-H2CO'
impbcor(imagename=myimagebase+'.image',
        pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor',
        overwrite=True) # perform PBcorr

exportfits(imagename=myimagebase+'.image.pbcor',
           fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',
           fitsimage=myimagebase+'.flux.fits') # export the corrected image

myimagebase = 'w51-C18O'
impbcor(imagename=myimagebase+'.image',
        pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor',
        overwrite=True) # perform PBcorr

exportfits(imagename=myimagebase+'.image.pbcor',
           fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',
           fitsimage=myimagebase+'.flux.fits') # export the corrected image

myimagebase = 'w51-CO'
impbcor(imagename=myimagebase+'.image',
        pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor',
        overwrite=True) # perform PBcorr

exportfits(imagename=myimagebase+'.image.pbcor',
           fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',
           fitsimage=myimagebase+'.flux.fits') # export the corrected image

