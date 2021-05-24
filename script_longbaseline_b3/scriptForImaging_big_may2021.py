import shutil

# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'

contimagename = 'w51north.spw0thru19.20000.robust0.0.thr0.5mJy.mfs.I.may21.2021'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists('cleanmask_north.crtf')

    tclean(vis=contvis,
           imagename=contimagename,
           field='w51n',
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
           imsize = 20000,
           cell= '0.005arcsec',
           weighting = 'briggs',
           robust = 0.0,
           niter = 200000,
           threshold = '0.5mJy',
           interactive = False,
           gridder = 'standard',
           pbcor = True,
           savemodel='none',
           mask='cleanmask_north.crtf')

contimagename = 'w51e2.spw0thru19.20000.robust0.0.thr0.5mJy.mfs.I.may21.2021'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists('cleanmask_e2.crtf')

    tclean(vis=contvis,
           imagename=contimagename,
           field='w51e2',
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
           imsize = 20000,
           cell= '0.005arcsec',
           weighting = 'briggs',
           robust = 0.0,
           niter = 200000,
           threshold = '0.5mJy',
           interactive = False,
           gridder = 'standard',
           pbcor = True,
           savemodel='none',
           mask='cleanmask_e2.crtf')




if not os.path.exists('calibrated_final_cont_selfcal.ms'):
    shutil.copytree('calibrated_final_cont.ms', 'calibrated_final_cont_selfcal.ms')

# Set the ms and continuum image name.
contvis = 'calibrated_final_cont_selfcal.ms'

contimagename = 'w51north.spw0thru19.20000.robust0.0.thr0.2mJy.mfs.I.may23.2021'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists('cleanmask_north.crtf')

    tclean(vis=contvis, imagename=contimagename, field='w51n', specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000, threshold =
           '0.2mJy', interactive = False, gridder = 'standard', pbcor = True,
           savemodel='none', mask='cleanmask_north.crtf')
    tclean(vis=contvis, imagename=contimagename, field='w51n', specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold = '0.2mJy',
           interactive = False, gridder = 'standard', pbcor = True,
           savemodel='modelcolumn', mask='cleanmask_north.crtf')

contimagename = 'w51e2.spw0thru19.20000.robust0.0.thr0.2mJy.mfs.I.may23.2021'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists('cleanmask_e2.crtf')

    tclean(vis=contvis, imagename=contimagename, field='w51e2', specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000, threshold =
           '0.2mJy', interactive = False, gridder = 'standard', pbcor = True,
           savemodel='none', mask='cleanmask_e2.crtf')
    tclean(vis=contvis, imagename=contimagename, field='w51e2', specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold = '0.2mJy',
           interactive = False, gridder = 'standard', pbcor = True,
           savemodel='modelcolumn', mask='cleanmask_e2.crtf')



gaincal(vis=contvis, caltable='w51_b3_selfcal_phase1_T.cal', calmode="p",
        gaintype="T", solint="inf", solnorm=True)

cals = ['w51_b3_selfcal_phase1_T.cal']

applycal(vis=contvis,
         gaintable=cals,
         interp="linear",
         applymode='calonly',
         calwt=False)

contimagename = 'w51north.spw0thru19.20000.robust0.0.thr0.1mJy.mfs.I.may23.2021.selfcal1'
threshold = '0.1mJy'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists('cleanmask_north.crtf')

    tclean(vis=contvis, imagename=contimagename, field='w51n', specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000,
           threshold=threshold, interactive = False, gridder = 'standard',
           pbcor = True, savemodel='none', mask='cleanmask_north.crtf')
    tclean(vis=contvis, imagename=contimagename, field='w51n', specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           savemodel='modelcolumn', mask='cleanmask_north.crtf')

contimagename = 'w51e2.spw0thru19.20000.robust0.0.thr0.2mJy.mfs.I.may23.2021'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists('cleanmask_e2.crtf')

    tclean(vis=contvis, imagename=contimagename, field='w51e2', specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 200000,
           threshold=threshold, interactive = False, gridder = 'standard',
           pbcor = True, savemodel='none', mask='cleanmask_e2.crtf')
    tclean(vis=contvis, imagename=contimagename, field='w51e2', specmode='mfs',
           deconvolver='mtmfs', nterms=2, imsize = 20000, cell= '0.005arcsec',
           weighting = 'briggs', robust = 0.0, niter = 1, threshold=threshold,
           interactive = False, gridder = 'standard', pbcor = True,
           savemodel='modelcolumn', mask='cleanmask_e2.crtf')

gaincal(vis=contvis, caltable='w51_b3_selfcal_phase2_T.cal', calmode="p",
        gaintype="T", solint="inf", solnorm=True)

