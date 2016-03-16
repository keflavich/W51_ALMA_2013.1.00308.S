import numpy as np

field='4~40' # science field(s). For a mosaic, select all mosaic fields. DO
             # NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT
             # WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
phasecenter='J2000 19h23m43.957 +14d30m34.67'
cell='0.05arcsec' # cell size for imaging.
imsize = [256,256] # size of image in pixels.

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean.

weighting = 'uniform'
robust=-2.0
threshold = '30.0mJy'

spws = {0: '0,4',
        1: '1,5',
        2: '2,6',
        3: '3,7',
       }
nchans_total = {0: 3840, 1: 3840, 2: 3840, 3: 3840}
frange = {0: [218136., 218575.],
          1: [218422., 220268.],
          2: [230436., 232310.],
          3: [233041., 234915.],
         }
fstep = {0:130., # kHz
         1:500., # kHz
         2:500., # kHz
         3:500., # kHz
        }
nchans_total = {ii: int(np.abs(np.diff(frange[ii])/fstep[ii]*1000.)[0])
                for ii in frange}

ncubes_per_window = 20


for spwnum in '3201':
    spwnum = int(spwnum)
    spw = spws[spwnum]

    concatvis = 'w51_concat.spw{0}.merge'.format(spwnum)
    if not os.path.exists(concatvis):
        print "# running cvel on all lines in spw{0}".format(spwnum)
        cvelvises = []
        for ss in spw.split(","):
            ss = int(ss)
            cvelvis = 'w51_concat.spw{0}.cvel'.format(ss)
            cvelvises.append(cvelvis)
            cvel(vis='w51_concat.ms.split.cal',
                 outputvis=cvelvis,
                 passall=False, field=field, spw=str(ss), selectdata=True,
                 timerange='', array='', antenna='', scan='', mode='frequency',
                 nchan=nchans_total[spwnum],
                 start='{0}MHz'.format(frange[spwnum][0]),
                 width='{0}kHz'.format(fstep[spwnum]), interpolation='linear',
                 phasecenter='', restfreq='', outframe='', veltype='radio',
                 hanning=False,)
        concat(vis=cvelvises,
               concatvis=concatvis,)
    else:
        print "Already cvel'd spw {0} to {1}".format(spwnum, concatvis)

    print "# running clean on all lines in spw{0}".format(spwnum)
    nchans_total_thiscube = nchans_total[spwnum]
    nchans_per_cube = int(nchans_total_thiscube/ncubes_per_window)
    for ii in range(ncubes_per_window):
        start = nchans_per_cube*ii
        end = nchans_per_cube*(ii+1)
        if end > nchans_total_thiscube:
            end = nchans_total_thiscube
        output = 'piece_of_full_W51e2cutout_cube_uniform.spw{0}.channels{1}to{2}'.format(spwnum, start, end)

        # Channel-based gridding has major bugs when dealing with CVEL'd data
        # It is therefore necessary to compute the frequency gridding by hand
        startfreq = "{0}GHz".format(frange[spwnum][0]/1e3 + start * fstep[spwnum]/1e6)
        width = "{0}kHz".format(fstep[spwnum])


        # LINE IMAGING (MOSAIC MODE)
        if not os.path.exists(output+".image"):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            tclean(vis = concatvis,
                   imagename = output,
                   field = '',
                   spw = '', # there should be only one
                   gridder = 'mosaic',
                   deconvolver = 'clark',
                   specmode = 'cube',
                   width = width,
                   start = startfreq,
                   nchan = nchans_per_cube,
                   veltype = 'radio',
                   outframe = 'LSRK',
                   interactive = F,
                   niter = 1000,
                   imsize = imsize,
                   cell = cell,
                   weighting = weighting,
                   phasecenter = phasecenter,
                   robust = robust,
                   threshold = threshold,
                   pblimit=0.4)

          
        myimagebase = output
        # I've given up on primary beam correction, at least for now
        exportfits(imagename=myimagebase+'.image',
                   fitsimage=myimagebase+'.image.fits',
                   overwrite=True,
                   dropdeg=True)
