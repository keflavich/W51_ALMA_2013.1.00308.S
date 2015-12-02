"""
Configuration of the calibrated 12m and 7m data
"""

finalvis12m='calibrated_12m.ms'
finalvis7m='calibrated_7m.ms'

import numpy as np

field='4~40' # science field(s). For a mosaic, select all mosaic fields. DO
             # NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT
             # WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
# backup if needed because of tclean bug phasecenter = "J2000 19:23:41.585 +14:30:41.00"
phasecenter=''
cell='.15arcsec' # cell size for imaging.
imsize = [960,960] # size of image in pixels.

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean.

weighting = 'briggs'
robust=0.5
threshold = '10.0mJy'

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
