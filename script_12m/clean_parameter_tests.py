
imsize = [3072,3072]
cell = '0.05arcsec'
vis4 = 'w51_contvis_selfcal_4.ms'
phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"

myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='15mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

# 5 mJy clark diverges horribly (but retry anyway)
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_5mJy"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='5mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

# 5 mJy clark diverges horribly (so try 10)
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_10mJy"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='10mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


# 10mJy hogbom
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_10mJy_hogbom"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='hogbom', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='10mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

# retry 5 mjy with hogbom
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_5mJy_hogbom"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='hogbom', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='5mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

# 1 mJy hogbom?
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_1mJy"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='hogbom', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='10mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_scales_0_5"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       scales = [0,5],
       deconvolver='multiscale', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='10mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_scales_0_5_15"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       scales = [0,5,15],
       deconvolver='multiscale', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='10mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)







###### CELL SIZE EXPERIMENTS ######

cell = '0.04arcsec'
imsize = [3840,3840]
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_smallerpix"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='15mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


cell = '0.06arcsec'
imsize = [2560,2560]
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_biggerpix"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='15mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


cell = '0.03arcsec'
imsize = [5120,5120]
myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper_tinypix"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='15mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

