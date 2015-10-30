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
nchans_total = {ii: 3840 for ii in range(8)}
ncubes_per_window = 20
finalvis='w51_pointing32.split.cal'
# don't know how to contsub yet
linevis = finalvis#+'.contsub'

cell={x:'.07' for x in (0,1,2,3)}
cell.update({x:'0.05arcsec' for x in (4,5,6,7)}) # cell size for imaging.
imsize = [256,256] # size of image in pixels.
imsize={x:[512,512] for x in (0,1,2,3)}
imsize.update({x:[512,512] for x in (4,5,6,7)})

freqs = {0: 218604.028,
         1: 220296.833,
         2: 230435.532,
         3: 233040.032,}
deltafreq = {0:-0.122070,
             1:-0.488281,
             2: 0.488281,
             3: 0.488281,
            }

for spwnum in '1230':
    print "# running clean on all lines in spw{0}".format(spwnum)
    spw = spws[int(spwnum)] # this tries to merge incompatible channels
    #spw = int(spwnum)
    nchans_total_thiscube = nchans_total[int(spwnum)]
    nchans_per_cube = nchans_total_thiscube/ncubes_per_window
    inputvis = linevis
    for ii in range(ncubes_per_window):
        start = nchans_per_cube*ii
        end = nchans_per_cube*(ii+1)
        output = 'piece_of_w51pointing32_cube_bothSB.spw{0}.channels{1}to{2}'.format(spwnum, start, end)
        #---------------------------------------------------
        # LINE IMAGING (MOSAIC MODE)
        if not os.path.exists(output+".image"):
            stfrq = freqs[int(spwnum)] + start*deltafreq[int(spwnum)]
            endfrq = freqs[int(spwnum)] + end*deltafreq[int(spwnum)]
            print("Imaging {0} with freqs {1} MHz-{2} MHz,"
                  " {3} channels, df={4}MHz".format(output, stfrq, endfrq,
                                                    nchans_per_cube,
                                                    deltafreq[int(spwnum)]))
            os.system('rm -rf ' + output + '*')
            clean(vis = inputvis,
                  imagename = output,
                  field = field,
                  spw = str(spw),
                  imagermode = 'csclean',
                  mode = 'frequency',
                  width = "{0}MHz".format(deltafreq[int(spwnum)]),
                  start = "{0}MHz".format(stfrq),
                  nchan = nchans_per_cube,
                  chaniter = True,
                  veltype = 'radio',
                  outframe = 'LSRK',
                  interactive = F,
                  niter = niter,
                  imsize = imsize[int(spwnum)],
                  cell = cell[int(spwnum)],
                  psfmode='clark',
                  weighting = weighting,
                  phasecenter = phasecenter,
                  robust = robust,
                  threshold = threshold,
                  pbcor = F,
                  usescratch= F)

          
        if not (os.path.exists(output+".image.pbcor.fits") or os.path.exists(output+".image.fits")):
            myimagebase = output
            if os.path.exists(output+".flux"):
                impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
                exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True,dropdeg=True)
                exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True,dropdeg=True)
            else:
                exportfits(imagename=myimagebase+'.image', fitsimage=myimagebase+'.image.fits', overwrite=True,dropdeg=True)


for spwnum in '12356704':
    print "# running clean on all lines in spw{0}".format(spwnum)
    #spw = spws[int(spwnum)] # this tries to merge incompatible channels
    spw = int(spwnum)
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
                  spw = str(spw),
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
                  imsize = imsize[spw],
                  cell = cell[spw],
                  psfmode='clark',
                  weighting = weighting,
                  phasecenter = phasecenter,
                  robust = robust,
                  threshold = threshold,
                  pbcor = F,
                  usescratch= F)

          
        if not (os.path.exists(output+".image.pbcor.fits") or os.path.exists(output+".image.fits")):
            myimagebase = output
            if os.path.exists(output+".flux"):
                impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
                exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True,dropdeg=True)
                exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True,dropdeg=True)
            else:
                exportfits(imagename=myimagebase+'.image', fitsimage=myimagebase+'.image.fits', overwrite=True,dropdeg=True)
