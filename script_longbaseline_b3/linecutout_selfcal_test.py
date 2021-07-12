import numpy as np
import os
import glob
import datetime
import sys
sys.path.append('.')
from source_ids import sources, sources_fmtd, source_field_mapping

def makefits(myimagebase, cleanup=True):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image

    if cleanup:
        for suffix in ('psf', 'weight', 'sumwt', 'pb', 'model', 'residual',
                       'mask', 'image', 'workdirectory'):
            os.system('rm -rf {0}.{1}'.format(myimagebase, suffix))




mslist = ['calibrated_final.ms']

for ms in mslist:
    listobs(ms, listfile=ms+'.listobs', overwrite=True)

for sourcename, coordinate in sources_fmtd.items():

    for spw,spws in enumerate([(0,4,8,12,16), (1,5,9,13,17), (2,6,10,14,18), (3,7,11,15,19)]):

        for suffix, niter in (('clarkclean1000', 1000), ):

            for robust in (0.0, ):

                imagename = 'W51{3}_only.B3.robust{2}.spw{0}.{1}'.format(spw, suffix, robust, sourcename)


sourcename = 'w51n'
coordinate = sources_fmtd[sourcename]

for datacolumn in ('data','corrected'):
    imagename = 'W51{3}_only.selfcaltest.B3.robust{2}.spw{0}.{1}.{4}'.format(spw, suffix, robust, sourcename, datacolumn)

    if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
        if fits.getheader("{0}.image.pbcor.fits".format(imagename))['RADESYS'] == 'ICRS':
            print("Skipping completed file {0}".format(imagename))
            continue
        else:
            print("Redoing {0} because it's in fk5.".format(imagename))

    print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
    tclean(vis=mslist,
           imagename=imagename,
           spw=",".join(['{0}'.format(ss) for ss in spws]),
           field=source_field_mapping[sourcename],
           specmode='cube',
           outframe='LSRK',
           threshold='15mJy',
           imsize=[128, 128],
           cell=['0.008arcsec'],
           niter=niter,
           #cycleniter=-1, # -1 is default
           #cyclefactor=0.0001, # set very small: try to prevent major cycles
           phasecenter=coordinate,
           deconvolver='clark',
           gridder='standard',
           weighting='briggs',
           robust=robust,
           pbcor=True,
           pblimit=0.2,
           savemodel='none',
           #chanchunks=1,
           datacolumn=datacolumn,
           #parallel=True,
           interactive=False)
    makefits(imagename)

