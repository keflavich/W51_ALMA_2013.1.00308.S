import os
# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'

#field='w51n'
#cleanmask='cleanmask_north.crtf'
#field='w51e2'
#cleanmask='cleanmask_e2.crtf'

#skip this
# contimagename = field+'.spw0thru19.20000.robust0.0.thr0.1mJy.mfs.I'
# 
# if not os.path.exists(contimagename+'.image.tt0.pbcor'):
#     for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
#         if os.path.exists(contimagename+ext):
#             rmtables(contimagename+ext)
# 
#     assert os.path.exists(cleanmask)
# 
#     tclean(vis=contvis,
#            imagename=contimagename,
#            field=field,
#            specmode='mfs',
#            deconvolver='mtmfs',
#            nterms=2,
#            imsize = 20000,
#            scales=[0,6,18],
#            cell= '0.005arcsec',
#            weighting = 'briggs',
#            robust = 0.0,
#            niter = 100000,
#            threshold = '0.1mJy',
#            interactive = False,
#            gridder = 'standard',
#            pbcor = True,
#            savemodel='none',
#            mask=cleanmask)


# Set the ms and continuum image name.
contvis = 'calibrated_final_cont_'+field+'selfcal.ms'
if not os.path.exists(contvis):
    split("calibrated_final_cont.ms", contvis, field=field)


contimagename = field+'.spw0thru19.20000.robust0.0.thr0.075mJy.mfs.I'
threshold='0.075mJy' # noise ~0.015 mJy

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists(cleanmask)

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           scales=[0,6,18],
           savemodel='none', mask=cleanmask)
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           scales=[0,6,18],
           savemodel='modelcolumn', mask=cleanmask)

    gaincal(vis=contvis, field=field, caltable=field+'_b3_selfcal_phase1_T.cal', calmode="p",
            gaintype="T", solint="inf", solnorm=True)

cals = [field+'_b3_selfcal_phase1_T.cal']

applycal(vis=contvis,
         gaintable=cals,
         gainfield=[field],
         interp="linear",
         applymode='calonly',
         calwt=False)

contimagename = field+'.spw0thru19.20000.robust0.0.thr0.075mJy.mfs.I.selfcal1'
threshold = '0.075mJy'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists(cleanmask)

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000,
           threshold=threshold, interactive = False, gridder = 'standard',
           scales=[0,6,18],
           pbcor = True, savemodel='none', mask=cleanmask)
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           scales=[0,6,18],
           savemodel='modelcolumn', mask=cleanmask)

    gaincal(vis=contvis, field=field, caltable=field+'_b3_selfcal_phase2_T.cal',
            gaintable=[field+'_b3_selfcal_phase1_T.cal'], calmode="p",
            gaintype="T", solint="inf", solnorm=True)

