import os

if not os.path.exists('w51_pointing32.split.cal'):
    split(
          vis                 = 'w51_concat.ms.split.cal', #  Name of input measurement set
          outputvis           = 'w51_pointing32.split.cal', #  Name of output measurement set
          datacolumn          =     'data',       #  Which data column(s) to split out
          field               =       '32',        #  Select field using ID(s) or name(s)
    )

field='w51'
phasecenter=''
cell='.10arcsec' # cell size for imaging.
imsize = [256,256] # size of image in pixels.

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
finalvis='w51_pointing32.split.cal'
# don't know how to contsub yet
linevis = finalvis#+'.contsub'

for spwnum in '0123':
    print "# running clean on all lines in spw{0}".format(spwnum)
    spw = spws[int(spwnum)]
    nchans_total_thiscube = nchans_total[int(spwnum)]
    nchans_per_cube = nchans_total_thiscube/ncubes_per_window
    inputvis = linevis
    for ii in range(ncubes_per_window):
        start = nchans_per_cube*ii
        end = nchans_per_cube*(ii+1)
        output = 'piece_of_w51pointing32_cube.spw{0}.channels{1}to{2}'.format(spwnum, start, end)
        #---------------------------------------------------
        # LINE IMAGING (MOSAIC MODE)
        if not os.path.exists(output+".image"):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            clean(vis = inputvis,
                  imagename = output,
                  field = field,
                  spw = spw,
                  imagermode = 'csclean',
                  mode = 'channel',
                  width = 1,
                  start = start,
                  nchan = nchans_per_cube,
                  chaniter = True,
                  veltype = 'radio',
                  outframe = 'LSRK',
                  interactive = F,
                  niter = niter,
                  imsize = imsize,
                  cell = cell,
                  psfmode='clark',
                  weighting = weighting,
                  phasecenter = phasecenter,
                  robust = robust,
                  threshold = threshold,
                  pbcor = T,
                  usescratch= F)

          
        if not os.path.exists(output+".image.pbcor"):
            myimagebase = output
            impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
            exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
            exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)


