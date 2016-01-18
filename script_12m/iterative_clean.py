import os

imsize = [3072,3072]
cell = '0.05arcsec'
vis4 = 'w51_contvis_selfcal_4.ms'
phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"

if not os.path.exists('dirty_50mJy.mask'):
    myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_dirty"
    os.system('rm -rf {0}.*'.format(myimagebase))
    tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
           deconvolver='clark', gridder='mosaic', outframe='LSRK',
           pblimit=0.4, interpolation='linear',
           interactive=False, niter=0,
           threshold='15mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
           weighting='briggs', savemodel='modelcolumn', robust=-2.0)
    exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
    exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

    dirtyimage = myimagebase+".residual"
    ia.open(dirtyimage)
    ia.calcmask(mask=dirtyimage+" > 0.05", name='dirty_mask_50mJy')
    ia.close()
    makemask(mode='copy', inpimage=dirtyimage,
             inpmask=dirtyimage+":dirty_mask_50mJy", output='dirty_50mJy.mask')
    exportfits('dirty_50mJy.mask', 'dirty_50mJy.mask.fits', dropdeg=True, overwrite=True)

# iterclean...
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_20mJy_mask_iter0"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='hogbom', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='20mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       mask='dirty_50mJy.mask',
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

residualimage = myimagebase+".residual"
ia.open(residualimage)
ia.calcmask(mask=residualimage+" > 0.02", name='clean_mask_20mJy')
ia.close()
makemask(mode='copy', inpimage=residualimage,
         inpmask=residualimage+":clean_mask_20mJy", output='clean_20mJy.mask')
exportfits('clean_20mJy.mask', 'clean_20mJy.mask.fits', dropdeg=True, overwrite=True)

myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_10mJy_mask_iter1"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='hogbom', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='10mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       mask='dirty_20mJy.mask',
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

residualimage = myimagebase+".residual"
ia.open(residualimage)
ia.calcmask(mask=residualimage+" > 0.01", name='clean_mask_10mJy')
ia.close()
makemask(mode='copy', inpimage=residualimage,
         inpmask=residualimage+":clean_mask_10mJy", output='clean_10mJy.mask')
exportfits('clean_10mJy.mask', 'clean_10mJy.mask.fits', dropdeg=True, overwrite=True)


myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_1mJy_mask_iter2"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='hogbom', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='1mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       mask='dirty_10mJy.mask',
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)
