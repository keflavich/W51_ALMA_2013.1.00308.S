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
robust=0.5
threshold = '0.0mJy'

spws = {0: '0,4',
        1: '1,5',
        2: '2,6',
        3: '3,7',
       }
nchans_total = {0: 3840, 1: 3840, 2: 3840, 3: 3840}
frange = {0: [218135.2792,218575.868],
          1: [218421.83396,220268.545],
          2: [230435.532, 232310.53104],
          3: [233040.032,234915.03104],
         }
fstep = {0:130, # kHz
         1:500, # kHz
         2:500, # kHz
         3:500, # kHz
        }
nchans_total = {ii: int(np.abs(np.diff(frange[ii])/fstep[ii]*1000)[0])
                for ii in frange}

ncubes_per_window = 20
finalvis='w51_concat.ms.split.cal'
# don't know how to contsub yet
linevis = finalvis#+'.contsub'


for spwnum in '3210':
    swpnum = int(spwnum)
    spw = spws[spwnum]

    print "# running cvel on all lines in spw{0}".format(spwnum)
    cvel(vis='w51_concat.ms.split.cal',
         outputvis='w51_concat.spw{0}.cvel'.format(spwnum),
         passall=False, field=field spw=spw, selectdata=True, timerange='',
         array='', antenna='', scan='', mode='frequency',
         nchan=nchans_total[spwnum], start='{0}MHz'.format(frange[spwnum][0]),
         width=fstep[spwnum], interpolation='linear', phasecenter='',
         restfreq='', outframe='', veltype='radio', hanning=False,)

    print "# running clean on all lines in spw{0}".format(spwnum)
    nchans_total_thiscube = nchans_total[spwnum]
    nchans_per_cube = int(nchans_total_thiscube/ncubes_per_window)
    inputvis = linevis
    for ii in range(ncubes_per_window):
        start = nchans_per_cube*ii
        end = nchans_per_cube*(ii+1)
        if end > nchans_total_thiscube:
            end = nchans_total_thiscube
        output = 'piece_of_full_W51_cube.spw{0}.channels{1}to{2}'.format(spwnum, start, end)
        #---------------------------------------------------
        # LINE IMAGING (MOSAIC MODE)
        if not os.path.exists(output+".image"):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            clean(vis = inputvis,
                  imagename = output,
                  field = '',
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
                  pbcor = F,
                  usescratch= F)

          
        if not os.path.exists(output+".image.pbcor"):
            myimagebase = output
            impbcor(imagename=myimagebase+'.image',
                    pbimage=myimagebase+'.flux',
                    outfile=myimagebase+'.image.pbcor', overwrite=True)
            exportfits(imagename=myimagebase+'.image.pbcor',
                       fitsimage=myimagebase+'.image.pbcor.fits',
                       overwrite=True,
                       dropdeg=True)
            exportfits(imagename=myimagebase+'.flux',
                       fitsimage=myimagebase+'.flux.fits', overwrite=True,
                       dropdeg=True)
