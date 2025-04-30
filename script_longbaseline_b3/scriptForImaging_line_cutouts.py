from astropy.io import fits
import numpy as np
import os
import glob
import datetime
import sys
sys.path.append('.')
from source_ids import sources_fmtd, source_field_mapping

def makefits(myimagebase, cleanup=True):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    # exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
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

        for suffix, niter in (('clarkclean1e5', int(1e5)), ):

            for robust in (0.5, ):

                imagename = 'W51{3}_only.B3.robust{2}.spw{0}.{1}.1536.bigpix'.format(spw, suffix, robust, sourcename)


                if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
                    if fits.getheader("{0}.image.pbcor.fits".format(imagename))['RADESYS'] == 'ICRS':
                        print("Skipping completed file {0}".format(imagename))
                        continue
                    else:
                        print("Redoing {0} because it's in fk5.".format(imagename))

                # notes from 2025-04-26:
                # 128 -> 1536
                # 1e4 -> 1e5 iter
                # 0.008 -> 0.01 cell size [minor axis was sampled by 4 diagonal pixels before]
                # cyclefactor not forced to skip major cycles
                # threshold 15 -> 2.5mJy (roughly 5-sigma)
                print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
                tclean(vis=mslist,
                       imagename=imagename,
                       datacolumn='data',
                       spw=",".join(['{0}'.format(ss) for ss in spws]),
                       field=source_field_mapping[sourcename],
                       specmode='cube',
                       outframe='LSRK',
                       threshold='2.5mJy',
                       imsize=[1536, 1536],
                       cell=['0.02arcsec'],
                       niter=niter,
                       cycleniter=-1, # -1 is default
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
                       #parallel=True,
                       interactive=False)
                makefits(imagename)


"""
Work log 2025-04-26

cd /red/adamginsburg/w51
rsync -ravpu --partial --progress /orange/adamginsburg/w51/2017.1.00293.S/calibrated_final.ms.tar .
tar -xvf calibrated_final.ms.tar


"""