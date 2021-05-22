
# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'

contimagename = 'w51north_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.20000.robust0.0.thr0.5mJy.mfs.I.may21.2021'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

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

contimagename = 'w51north_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.20000.robust0.0.thr0.5mJy.mfs.I.may21.2021'

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

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
