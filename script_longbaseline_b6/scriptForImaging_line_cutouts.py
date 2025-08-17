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


# ./alma/2015.1.01596.S/archive/2015.1.01596.S/science_goal.uid___A001_X2d6_X10/group.uid___A001_X2d6_X11/member.uid___A001_X2d6_X12/script/w51e2_w51n.ms

mslist = ['w51e2_w51n.ms']

# 0  1  2  3  4  5  6  7  8  9
# 10 11 12 13 14 15 16 17 18 19
# 20 21 22 23 24 25 26 27 28 29
spwlist = [(0, 10, 20), (1, 11, 21), (2, 12, 22), (3, 13, 23), (4, 14, 24), (5, 15, 25), (6, 16, 26), (7, 17, 27), (8, 18, 28), (9, 19, 29)]

for ms in mslist:
    listobs(ms, listfile=ms+'.listobs', overwrite=True)

for sourcename, coordinate in sources_fmtd.items():

    for spw,spws in enumerate(spwlist):

        for suffix, niter in (('clarkclean1e5', int(1e5)), ):

            for robust in (0.5, ):

                imagename = '{3}_only.B6.robust{2}.spw{0}.{1}.1536.bigpix'.format(spw, suffix, robust, sourcename)


                if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
                    if fits.getheader("{0}.image.pbcor.fits".format(imagename))['RADESYS'] == 'ICRS':
                        print("Skipping completed file {0}".format(imagename))
                        continue
                    else:
                        print("Redoing {0} because it's in fk5.".format(imagename))

                # these are very big pixels for B6....
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

for sourcename, coordinate in sources_fmtd.items():

    for spw,spws in enumerate(spwlist):

        for suffix, niter in (('clarkclean1e5', int(1e5)), ):

            for robust in (0.5, ):

                imagename = '{3}_only.B6.robust{2}.spw{0}.{1}.1536.smallpix'.format(spw, suffix, robust, sourcename)


                if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
                    if fits.getheader("{0}.image.pbcor.fits".format(imagename))['RADESYS'] == 'ICRS':
                        print("Skipping completed file {0}".format(imagename))
                        continue
                    else:
                        print("Redoing {0} because it's in fk5.".format(imagename))

                # these are very big pixels for B6....
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
                       cell=['0.01arcsec'],
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