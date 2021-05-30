import os
# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'

print("field", field)
print("cleanmask", cleanmask)


# Set the ms and continuum image name.
contvis = 'calibrated_final_cont_'+field+'selfcal.ms'
if not os.path.exists(contvis):
    split(vis="calibrated_final_cont.ms", 
          outputvis=contvis,
          field=field,
          datacolumn='data')


preselfcal_imagename = contimagename = field+'.spw0thru19.20000.robust0.0.thr0.1mJy.mfs.I'
threshold='0.1mJy' # noise ~0.015 mJy

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

caltable = field+'_b3_selfcal_phase1_T.cal'
if not os.path.exists(caltable):
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           scales=[0,6,18],
           savemodel='modelcolumn')

    gaincal(vis=contvis, field=field, caltable=caltable, calmode="p",
            gaintype="T", solint="inf", solnorm=True)

cals = [caltable]

applycal(vis=contvis,
         gaintable=cals,
         gainfield=[field],
         interp="linear",
         applymode='calonly',
         calwt=False)



# first iteration of self-calibrated imaging

prevcal_imagename = contimagename = field+'.spw0thru19.20000.robust0.0.thr0.075mJy.mfs.I.selfcal1'
threshold = '0.075mJy'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists(cleanmask)

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           startmodel=[preselfcal_imagename + ".model.tt0",
                       preselfcal_imagename + ".model.tt1"],
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000,
           threshold=threshold, interactive = False, gridder = 'standard',
           scales=[0,6,18],
           pbcor = True, savemodel='none', mask=cleanmask)

caltable = field+'_b3_selfcal_phase2_T.cal'
if not os.path.exists(caltable):
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           scales=[0,6,18],
           savemodel='modelcolumn')

    gaincal(vis=contvis, field=field, caltable=caltable,
            gaintable=cals, calmode="p",
            gaintype="T", solint="inf", solnorm=True)


cals += [caltable]


applycal(vis=contvis,
         gaintable=cals,
         gainfield=[field]*len(cals),
         interp="linear",
         applymode='calonly',
         calwt=False)


prevcal_imagename = contimagename = field+'.spw0thru19.20000.robust0.0.thr0.075mJy.mfs.I.selfcal2'
threshold = '0.075mJy'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists(cleanmask)

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           startmodel=[prevcal_imagename + ".model.tt0",
                       prevcal_imagename + ".model.tt1"],
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000,
           threshold=threshold, interactive = False, gridder = 'standard',
           scales=[0,6,18],
           pbcor = True, savemodel='none', mask=cleanmask)

caltable = field+'_b3_selfcal_phase3_T.cal'
if not os.path.exists(caltable):
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           scales=[0,6,18],
           savemodel='modelcolumn')

    gaincal(vis=contvis, field=field, caltable=caltable,
            gaintable=cals, calmode="p",
            gaintype="T", solint="inf", solnorm=True)


cals += [caltable]




applycal(vis=contvis,
         gaintable=cals,
         gainfield=[field]*len(cals),
         interp="linear",
         applymode='calonly',
         calwt=False)


prevcal_imagename = contimagename = field+'.spw0thru19.20000.robust0.0.thr0.075mJy.mfs.I.selfcal3'
threshold = '0.075mJy'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists(cleanmask)

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           startmodel=[prevcal_imagename + ".model.tt0",
                       prevcal_imagename + ".model.tt1"],
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000,
           threshold=threshold, interactive = False, gridder = 'standard',
           scales=[0,6,18],
           pbcor = True, savemodel='none', mask=cleanmask)

caltable = field+'_b3_selfcal_phase4_T_int.cal'
if not os.path.exists(caltable):
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           scales=[0,6,18],
           savemodel='modelcolumn')

    gaincal(vis=contvis, field=field, caltable=caltable,
            gaintable=cals, calmode="p",
            gaintype="T", solint="int", solnorm=True)


cals += [caltable]






applycal(vis=contvis,
         gaintable=cals,
         gainfield=[field]*len(cals),
         interp="linear",
         applymode='calonly',
         calwt=False)


prevcal_imagename = contimagename = field+'.spw0thru19.20000.robust0.0.thr0.075mJy.mfs.I.3term.selfcal4'
threshold = '0.075mJy'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists(cleanmask)

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           startmodel=[prevcal_imagename + ".model.tt0",
                       prevcal_imagename + ".model.tt1"],
           deconvolver='mtmfs', nterms=3, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000,
           threshold=threshold, interactive = False, gridder = 'standard',
           scales=[0,6,18],
           pbcor = True, savemodel='none', mask=cleanmask)

caltable = field+'_b3_selfcal_phase5_T_int.cal'
if not os.path.exists(caltable):
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', nterms=3, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           scales=[0,6,18],
           savemodel='modelcolumn')

    gaincal(vis=contvis, field=field, caltable=caltable,
            gaintable=cals, calmode="p",
            gaintype="T", solint="int", solnorm=True)


cals += [caltable]




