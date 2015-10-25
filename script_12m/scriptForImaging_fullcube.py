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
linevis = finalvis+'.contsub'

for spwnum in '3210':
    print "# running clean on all lines in spw{0}".format(spwnum)
    spw = spws[int(spwnum)]
    nchans_total_thiscube = nchans_total[int(spwnum)]
    nchans_per_cube = nchans_total_thiscube/ncubes_per_window
    inputvis = linevis
    for ii in range(ncubes_per_window):
        start = nchans_per_cube*ii
        end = nchans_per_cube*(ii+1)
        output = 'piece_of_full_W51_cube.spw{0}.channels{1}to{2}'.format(spwnum, start, end)
        #---------------------------------------------------
        # LINE IMAGING (MOSAIC MODE)
        if not os.path.exists(output+".image"):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            clean(vis = inputvis,
                  imagename = output,
                  field = field,
                  spw = spw,
                  imagermode = 'mosaic',
                  mode = 'channel',
                  width = 1,
                  start = start,
                  nchan = nchans_per_cube,
                  chaniter = True,
                  veltype = 'radio',
                  outframe = 'LSRK',
                  interactive = F,
                  niter = 2000,
                  imsize = imsize,
                  cell = cell,
                  psfmode='clark',
                  weighting = weighting,
                  phasecenter = phasecenter,
                  robust = robust,
                  threshold = threshold,
                  pbcor = T,
                  usescratch= T)

          
        if not os.path.exists(output+".image.pbcor"):
            myimagebase = output
            impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
            exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
            exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)


