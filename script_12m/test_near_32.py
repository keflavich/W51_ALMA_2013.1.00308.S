import os

field='w51'
phasecenter=''

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean.

weighting = 'briggs'
robust=0.0
threshold = '0.0mJy'

spws = {0: '0,4',
        1: '1,5',
        2: '2,6',
        3: '3,7',
       }
nchans_total = {ii: 3840 for ii in range(8)}
ncubes_per_window = 20
finalvis='w51_concat.ms.split.cal'
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

output = 'test_near_32_clean'
os.system('rm -rf {0}*'.format(output))
print("Cleaning")
clean(vis = finalvis,
      imagename = output,
      field = field,
      spw = '3:942~946,7:1002~1006',
      imagermode = 'csclean',
      mode = 'mfs',
      chaniter = True,
      veltype = 'radio',
      outframe = 'LSRK',
      interactive = F,
      niter = 5000,
      imsize = [512,512],
      cell = '0.07arcsec',
      psfmode='clark',
      weighting = weighting,
      phasecenter = 32,
      robust = robust,
      threshold = threshold,
      pbcor = F,
      usescratch= F)

print("Exporting")
exportfits(output+".image", output+".image.fits", dropdeg=True)

output = 'test_near_32_clean_freq'
os.system('rm -rf {0}*'.format(output))
print("Cleaning")
clean(vis = finalvis,
      imagename = output,
      field = field,
      spw = '3,7',
      imagermode = 'csclean',
      mode = 'frequency',
      start = '233.50GHz',
      width = '0.5MHz',
      nchan = 100,
      chaniter = True,
      veltype = 'radio',
      outframe = 'LSRK',
      interactive = F,
      niter = 5000,
      imsize = [512,512],
      cell = '0.07arcsec',
      psfmode='clark',
      weighting = weighting,
      phasecenter = 32,
      robust = robust,
      threshold = threshold,
      pbcor = F,
      usescratch= F)

print("Exporting")
exportfits(output+".image", output+".image.fits", dropdeg=True)
