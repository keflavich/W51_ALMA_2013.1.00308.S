"""
Joint imaging of 7m and 12m data
"""
import os
from line_to_image_list import line_to_image_list

calpath = './'

vistemplate = 'calibrated_{0}.ms'
vis_7m='calibrated_7m.ms'
vis_tc='calibrated_12m.ms'

velocity_range = 25,95
velocity_res = 1.0

# It's not clear how cvel will handle overlapping SPWs

for line, restfreq, velocity_res in line_to_image_list:

    outms_template = "{line}_W51_B6_{array}.cvel.ms"
    concatvis = "{line}_W51_B6_cvel_merge.ms".format(line=line)

    nchans = int((velocity_range[1]-velocity_range[0])/velocity_res)

    # cvel the data first
    for array in ('7m','12m'):
        outputvis = outms_template.format(line=line, array=array)
        if not os.path.exists(outputvis):
            cvel2(vis=os.path.join(calpath, vistemplate.format(array)),
                  outputvis=outputvis,
                  datacolumn='data',
                  field='w51',
                  mode='velocity',
                  nchan=nchans,
                  start='{0}km/s'.format(velocity_range[0]),
                  width='{0}km/s'.format(velocity_res),
                  interpolation='linear',
                  phasecenter='',
                  restfreq=restfreq,
                  outframe='LSRK',
                  veltype='radio',
                  hanning=False,)

    if not os.path.exists(concatvis):
        concat(vis=[outms_template.format(line=line, array=array) for array in ('7m','12m')],
               concatvis=concatvis)

    assert os.path.exists(concatvis)

    output = 'W51_b6_7M_12M.{0}'.format(line)

    if not os.path.exists(output+".image.pbcor.fits"):
        os.system('rm -rf ' + output + '*/')
        tclean(vis=concatvis,
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
               uvtaper=True,
               innteraper=['10arcsec'],
              )
        myimagebase = output
        impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
        exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True, dropdeg=True)
        exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', overwrite=True, dropdeg=True)
        exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', overwrite=True, dropdeg=True)

        for suffix in ('pb', 'weight', 'sumwt', 'psf', 'model', 'mask',
                       'image', 'residual'):
            os.system('rm -rf {0}.{1}'.format(output, suffix))

