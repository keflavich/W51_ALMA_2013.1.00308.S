finalvis7m='calibrated_7m.ms'

import numpy as np

field='w51'
phasecenter=''
cell='.6arcsec' # cell size for imaging.
imsize = [256,256] # size of image in pixels.

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean.

weighting = 'briggs'
robust=1.0
threshold = '100.0mJy'

spws_7m = {0: '0',
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


for spwnum in '3201':
    spwnum = int(spwnum)

    output = 'fullcube_7m_spw{0}'.format(spwnum)

    print "Imaging {0}".format(output)
    os.system('rm -rf ' + output + '*')
    tclean(vis = finalvis7m,
           imagename = output,
           field = '',
           spw = '', # there should be only one
           gridder='mosaic',
           specmode = 'cube',
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
