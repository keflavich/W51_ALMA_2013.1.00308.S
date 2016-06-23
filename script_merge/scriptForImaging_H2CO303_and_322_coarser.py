"""
Re-image the 303 and 322 cubes with coarser intrinsic resolution; it has a real
effect.

Have to downsample (average) and then regrid the .ms first
"""
import os
from line_to_image_list import line_to_image_list
ltd = dict(((x[0], x[1:]) for x in line_to_image_list))

for line in ('H2CO303_202', 'H2CO322_221'):
    vis12m = "{line}_W51_B6_12m.cvel.ms".format(line=line)
    concatvis = "{line}_W51_B6_cvel_merge.ms".format(line=line)
    vis12m_out = "{line}_W51_B6_12m.coarse.cvel.ms".format(line=line)
    concatvis_out = "{line}_W51_B6_cvel_merge.coarse.ms".format(line=line)

    velocity_range = 25,95
    velocity_res=1.2 # original is 0.5, regridding to 1.2
    nchans = int((velocity_range[1]-velocity_range[0])/velocity_res)

    restfreq = ltd[line][0]

    if not os.path.exists(vis12m_out):
        print("Regridding {0}".format(vis12m))
        split(vis=vis12m,
              outputvis=vis12m+".split",
              width=2,
              datacolumn='data',
             )
        cvel2(vis=vis12m+".split",
              outputvis=vis12m_out,
              datacolumn='data',
              field='w51',
              mode='velocity',
              spw='',
              nchan=nchans,
              start='{0}km/s'.format(velocity_range[0]),
              width='{0}km/s'.format(velocity_res),
              interpolation='linear',
              phasecenter='',
              restfreq=restfreq,
              outframe='LSRK',
              veltype='radio',
              hanning=False,)

    if not os.path.exists(concatvis_out):
        print("Regridding {0}".format(concatvis))
        split(vis=concatvis,
              outputvis=concatvis+".split",
              width=2,
              datacolumn='data',
             )
        cvel2(vis=concatvis+".split",
              outputvis=concatvis_out,
              datacolumn='data',
              field='w51',
              mode='velocity',
              spw='',
              nchan=nchans,
              start='{0}km/s'.format(velocity_range[0]),
              width='{0}km/s'.format(velocity_res),
              interpolation='linear',
              phasecenter='',
              restfreq=restfreq,
              outframe='LSRK',
              veltype='radio',
              hanning=False,)


    output = 'W51_b6_7M_12M.{0}.coarser'.format(line)
    if not os.path.exists(output+".image.pbcor.fits"):
        os.system('rm -rf ' + output + '*/')
        tclean(vis=concatvis_out,
               imagename=output,
               field='w51',
               spw='',
               gridder='mosaic',
               specmode='cube',
               start='{0}km/s'.format(velocity_range[0]),
               width='{0}km/s'.format(velocity_res), interpolation='linear',
               nchan=nchans,
               restfreq=restfreq,
               veltype='radio',
               outframe='LSRK',
               interactive=F,
               niter=1000,
               imsize=[3200,3200],
               cell='.05arcsec',
               weighting='briggs',
               robust=0.5,
               phasecenter='',
               threshold='15mJy',
               savemodel='none',
              )
        myimagebase = output
        impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
        exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True, dropdeg=True)
        exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', overwrite=True, dropdeg=True)
        exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', overwrite=True, dropdeg=True)

        for suffix in ('pb', 'weight', 'sumwt', 'psf', 'model', 'mask',
                       'image', 'residual'):
            os.system('rm -rf {0}.{1}'.format(output, suffix))

    output = 'W51_b6_12M.{0}.coarser'.format(line)
    if not os.path.exists(output+".image.pbcor.fits"):
        os.system('rm -rf ' + output + '*/')
        tclean(vis=vis12m_out,
               imagename=output,
               field='w51',
               spw='',
               gridder='mosaic',
               specmode='cube',
               start='{0}km/s'.format(velocity_range[0]),
               width='{0}km/s'.format(velocity_res), interpolation='linear',
               nchan=nchans,
               restfreq=restfreq,
               veltype='radio',
               outframe='LSRK',
               interactive=F,
               niter=1000,
               imsize=[3200,3200],
               cell='.05arcsec',
               weighting='briggs',
               robust=0.5,
               phasecenter='',
               threshold='15mJy',
               savemodel='none',
              )
        myimagebase = output
        impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
        exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True, dropdeg=True)
        exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', overwrite=True, dropdeg=True)
        exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', overwrite=True, dropdeg=True)

        for suffix in ('pb', 'weight', 'sumwt', 'psf', 'model', 'mask',
                       'image', 'residual'):
            os.system('rm -rf {0}.{1}'.format(output, suffix))
