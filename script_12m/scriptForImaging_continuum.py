import os
field='4~40' # science field(s). For a mosaic, select all mosaic fields. DO
              # NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT
              # WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
phasecenter=''
cell='.15arcsec' # cell size for imaging.
imsize = [960,960] # size of image in pixels.

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean.

weighting = 'briggs'
robust=0.0
niter=5000
threshold = '0.0mJy'

spws = {0: '0,4',
        1: '1,5',
        2: '2,6',
        3: '3,7',
       }
nchans_total = {0: 3840, 1: 3840, 2: 3840, 3: 3840}
ncubes_per_window = 20
finalvis='w51_concat.ms.split.cal'
# don't know how to contsub yet
linevis = finalvis#+'.contsub'

for spwnum in '3210':
    print "# running clean on all lines in spw{0}".format(spwnum)
    spw = spws[int(spwnum)]
    nchans_total_thiscube = nchans_total[int(spwnum)]
    nchans_per_cube = nchans_total_thiscube/ncubes_per_window
    inputvis = linevis
    os.system("rm -rf w51_cont_spw{0}_hires.*".format(spwnum))
    imagename = "w51_cont_spw{0}_hires".format(spwnum)
    clean(vis = 'w51_concat.ms.split.cal',
          imagename = imagename,
          field = "w51",
          spw = spw,
          mode = 'mfs',
          outframe = 'LSRK',
          interpolation = 'linear',
          imagermode='mosaic',
          interactive = False,
          niter = 10000,
          threshold = '20.0mJy', # got to 11 mJy, but it showed really bad instabilities.  15 was still too deep
          imsize = [2560,2560],
          cell = '0.052arcsec',
          weighting = 'briggs',
          pbcor=False,
          robust = 0.0)

          
    if not os.path.exists(imagename+".image.pbcor"):
        myimagebase = imagename
        impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
        exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
        exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)


