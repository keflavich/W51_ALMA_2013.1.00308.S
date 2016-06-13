"""


whitespace is for symmetry with unmerged
"""
finalvis12m='calibrated_12m.ms'
finalvis7m='calibrated_7m.ms'

import numpy as np
from astropy.io import fits

field='w51'
phasecenter=''
cell='.05arcsec' # cell size for imaging.
imsize = [3200,3200] # size of image in pixels.

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean.

weighting = 'briggs'
robust=-0.5
threshold = '100.0mJy'

spws_12m = {0: '0,4',
            1: '1,5',
            2: '2,6',
            3: '3,7',
           }
spws_7m =  {0: '0',
            1: '1',
            2: '2',
            3: '3',
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

ncubes_per_window = 40


for spwnum in '2013':
    spwnum = int(spwnum)


    concatvis = 'w51_concat_7m12m.spw{0}.merge'.format(spwnum)
    if not os.path.exists(concatvis):
        print "# running cvel on all lines in spw{0}".format(spwnum)
        cvelvises = []
        spw = spws_12m[spwnum]
        for ss in spw.split(","):
            ss = int(ss)
            cvelvis12m = 'w51_12m.spw{0}.cvel'.format(ss)
            cvelvises.append(cvelvis12m)
            if not os.path.exists(cvelvis12m):
                print("cveling {0}".format(cvelvis12m))
                cvel(vis=finalvis12m,
                     outputvis=cvelvis12m,
                     passall=False, field=field, spw=str(ss), selectdata=True,
                     timerange='', array='', antenna='', scan='', mode='frequency',
                     nchan=nchans_total[spwnum],
                     start='{0}MHz'.format(frange[spwnum][0]),
                     width='{0}kHz'.format(fstep[spwnum]), interpolation='linear',
                     phasecenter='', restfreq='', outframe='LSRK', veltype='radio',
                     hanning=False,)
            else:
                print("skipping {0}".format(cvelvis12m))
        spw = spws_7m[spwnum]
        for ss in spw.split(","):
            ss = int(ss)
            cvelvis7m = 'w51_7m.spw{0}.cvel'.format(ss)
            cvelvises.append(cvelvis7m)
            if not os.path.exists(cvelvis7m):
                print("cveling {0}".format(cvelvis7m))
                cvel(vis=finalvis7m,
                     outputvis=cvelvis7m,
                     passall=False, field=field, spw=str(ss), selectdata=True,
                     timerange='', array='', antenna='', scan='', mode='frequency',
                     nchan=nchans_total[spwnum],
                     start='{0}MHz'.format(frange[spwnum][0]),
                     width='{0}kHz'.format(fstep[spwnum]), interpolation='linear',
                     phasecenter='', restfreq='', outframe='LSRK', veltype='radio',
                     hanning=False,)
            else:
                print("skipping {0}".format(cvelvis7m))
        concat(vis=cvelvises,
               concatvis=concatvis,)
    else:
        print "Already cvel'd spw {0} to {1}".format(spwnum, concatvis)

    print "# running clean on all lines in spw{0}".format(spwnum)
    nchans_total_thiscube = nchans_total[spwnum]
    nchans_per_cube = int(nchans_total_thiscube/ncubes_per_window)
    for ii in range(ncubes_per_window):
        # add 1 channel at start and 1 at end because tclean mistreats these channels
        start = nchans_per_cube*ii -1
        if start <= 0:
            start = 0
        end = nchans_per_cube*(ii+1) +1
        if end > nchans_total_thiscube:
            end = nchans_total_thiscube
        output = 'piece_of_full_W51_7m12m_cube_hires.spw{0}.channels{1}to{2}'.format(spwnum, start, end)

        # Channel-based gridding has major bugs when dealing with CVEL'd data
        # It is therefore necessary to compute the frequency gridding by hand
        startfreq = "{0}GHz".format(frange[spwnum][0]/1e3 + start * fstep[spwnum]/1e6)
        width = "{0}kHz".format(fstep[spwnum])


        # LINE IMAGING (MOSAIC MODE)
        if (not (os.path.exists(output+".image.fits") or
                 os.path.exists(output+".image.pbcor.fits"))
            or ('reclean' in locals() and reclean)):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            tclean(vis = concatvis,
                   imagename = output,
                   field = '',
                   spw = '', # there should be only one
                   gridder='mosaic',
                   specmode = 'cube',
                   width = width,
                   start = startfreq,
                   nchan = nchans_per_cube + 2, # 1 channel at either end for buffer
                   veltype = 'radio',
                   outframe = 'LSRK',
                   deconvolver='clark',
                   interactive = F,
                   niter = 500000, # huge niter: forcibly go to the threshold
                   # in principle, at least, this might help smooth over the
                   # band-edge issues
                   imsize = imsize,
                   cell = cell,
                   weighting = weighting,
                   phasecenter = phasecenter,
                   robust = robust,
                   threshold = threshold,
                   savemodel='none',
                   )

              
            myimagebase = output
            exportfits(myimagebase+'.image', myimagebase+'.image.fits',
                       dropdeg=True, overwrite=True)
            impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
                    outfile=myimagebase+'.image.pbcor', overwrite=True)
            exportfits(myimagebase+'.image.pbcor',
                       myimagebase+'.image.pbcor.fits', dropdeg=True,
                       overwrite=True)

            for suffix in ('psf', 'weight', 'sumwt', 'pb', 'model', 'residual',
                           'mask', 'image'):
                os.system('rm -rf {0}.{1}'.format(myimagebase, suffix))
        else:
            print("Skipped {0} because it was already completed".format(output))

        # don't crop.  Let's do that later so we stop losing data.
        #if (os.path.exists(output+".image.fits") and
        #    os.path.exists(output+".image.pbcor.fits")):
        #    print("Cropping {0}".format(output))
        #    # crop out the junk channels
        #    fh = fits.open(myimagebase+".image.pbcor.fits")
        #    ds = 1 if start > 0 else 0
        #    de = -1 if end < nchans_total_thiscube else 0
        #    newoutput = 'piece_of_full_W51_7m12m_cube_hires.spw{0}.channels{1}to{2}'.format(spwnum, start+ds, end+de)
        #    if fh[0].data.shape[0] > nchans_per_cube:
        #        fh[0].data = fh[0].data[ds:fh[0].data.shape[0]+de,:,:]
        #        fh[0].header['CRPIX3'] = fh[0].header['CRPIX3'] - ds
        #    else:
        #        print("Size of dim 0 was {0} and desired was {1} so no cropping was done"
        #              .format(fh[0].data.shape[0], nchans_per_cube))

        #    if (len(fh) == 2 and fh[1].data.size > fh[0].data.shape[0]) or len(fh) == 1:
        #        # get the beams
        #        print("Resizing beams from {0} to {1}".format(fh[1].data.size, fh[0].data.shape[0]))
        #        beams = fits.open(myimagebase+".image.fits")[1]
        #        beams.data = beams.data[ds:fh[0].data.shape[0]+de]
        #        if len(fh) == 2:
        #            fh[1] = beams
        #        else:
        #            fh.append(beams)

        #    fh.writeto(newoutput+".image.pbcor.fits", clobber=True)

        #    if os.path.exists(newoutput+".image.pbcor.fits") and fits.getdata(newoutput+".image.pbcor.fits").shape[0] == nchans_per_cube:
        #        print("would delete cubes here but skipping for safety.")
        #        ## a bit dangerous... but I already did it, oh well.  Only way to keep file sizes in check.
        #        #os.system('rm -rf {output}.image.fits'.format(output=output))
        #        #os.system('rm -rf {output}.image.pbcor.fits'.format(output=output))
