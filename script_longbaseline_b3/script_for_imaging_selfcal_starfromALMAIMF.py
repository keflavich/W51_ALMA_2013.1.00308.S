import os
# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'

print("field", field)
print("cleanmask", cleanmask)
print("startmodel: ", startmodel)


# Set the ms and continuum image name.
contvis = 'calibrated_final_cont_'+field+'selfcal_startmod.ms'
if not os.path.exists(contvis):
    split(vis="calibrated_final_cont.ms", 
          outputvis=contvis,
          field=field,
          datacolumn='data')


# beamsize: Beam : 0.0669461 arcsec, 0.0421181 arcsec, -44.3619 deg
# >>> 0.0421181 / 0.007 / 2**0.5
# 4.254570588670446
selfcalpars = {0: {'imaging': {'threshold': '0.1mJy', 'nterms': 2, 'robust': 0, 'weighting':'briggs',
                            'cell':'0.007arcsec', 'imsize':14500, 'scales':[0,6,18], 'niter':200000},
                'calibration': {'calmode': 'p', 'gaintype': 'T', 'solint': 'inf', 'solnorm': True},},

               1: {'imaging': {'threshold': '0.075mJy', 'nterms': 2, 'robust': 0,
                               'weighting':'briggs', 'cell':'0.007arcsec',
                               'imsize':14500, 'scales':[0,6,18],
                               'niter':200000},
                   'calibration': {'calmode': 'p', 'gaintype': 'T', 'solint': 'inf', 'solnorm': True},}
               2: {'imaging': {'threshold': '0.075mJy', 'nterms': 2, 'robust': 0,
                               'weighting':'briggs', 'cell':'0.007arcsec',
                               'imsize':14500, 'scales':[0,6,18],
                               'niter':200000},
                   'calibration': {'calmode': 'p', 'gaintype': 'T', 'solint': 'inf', 'solnorm': True},}
               3: {'imaging': {'threshold': '0.075mJy', 'nterms': 2, 'robust': 0,
                               'weighting':'briggs', 'cell':'0.007arcsec',
                               'imsize':14500, 'scales':[0,6,18],
                               'niter':200000},
                   'calibration': {'calmode': 'p', 'gaintype': 'G', 'solint': 'inf', 'solnorm': True},}
               4: {'imaging': {'threshold': '0.075mJy', 'nterms': 2, 'robust': 0,
                               'weighting':'briggs', 'cell':'0.007arcsec',
                               'imsize':14500, 'scales':[0,6,18],
                               'niter':200000},
                   'calibration': {'calmode': 'p', 'gaintype': 'T',
                                   'solint': 'int', 'solnorm': True},}
               5: {'imaging': {'threshold': '0.075mJy', 'nterms': 3, 'robust': 0,
                               'weighting':'briggs', 'cell':'0.007arcsec',
                               'imsize':14500, 'scales':[0,6,18],
                               'niter':200000},
                   'calibration': {'calmode': 'p', 'gaintype': 'T',
                                   'solint': 'int', 'solnorm': True},}
              }

selfcaliter=0
impars = selfcalpars[selfcaliter]['imaging']
calpars = selfcalpars[selfcaliter]['calibration']

dirtyimage = field+'.spw0thru19.{imsize}.robust{robust}.thr{threshold}.mfs.I.dirty'.format(**impars)

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists(cleanmask)

    impars = impars.copy()
    impars['niter'] = 0

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs',
           interactive = False, gridder = 'standard', pbcor = True,
           savemodel='none', mask=cleanmask,
           **impars
          )

for sm in startmodel:
    imregrid(imagename=sm, template=dirtyimage+".image.tt0", output=sm+".regrid")

startmodel = [sm+".regrid" for sm in startmodel]


impars = selfcalpars[selfcaliter]['imaging']
preselfcal_imagename = contimagename = field+'.spw0thru19.{imsize}.robust{robust}.thr{threshold}.mfs.I'.format(
    **impars)

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        if os.path.exists(contimagename+ext):
            rmtables(contimagename+ext)

    assert os.path.exists(cleanmask)

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs',
           interactive = False, gridder = 'standard', pbcor = True,
           savemodel='none', mask=cleanmask,
           startmodel=startmodel,
           **impars
          )

caltable = field+'_b3_selfcal_phase{selfcaliter}_T.cal'.format(selfcaliter=selfcaliter)
if not os.path.exists(caltable):
    impars['niter'] = 1
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', interactive=False, gridder='standard',
           pbcor=True, savemodel='modelcolumn', **impars)

    gaincal(vis=contvis, field=field, caltable=caltable, **calpars)

cals = [caltable]

applycal(vis=contvis,
         gaintable=cals,
         gainfield=[field],
         interp="linear",
         applymode='calonly',
         calwt=False)



for selfcaliter in selfcalpars:

    impars = selfcalpars[selfcaliter]['imaging']
    calpars = selfcalpars[selfcaliter]['calibration']
               
    prevcal_imagename = contimagename = field+'.spw0thru19.{imsize}.robust{robust}.thr{threshold}.mfs.I.selfcal{selfcaliter}'.format(
        selfcaliter=selfcaliter,
        **impars)

    if not os.path.exists(contimagename+'.image.tt0.pbcor'):
        for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
            if os.path.exists(contimagename+ext):
                rmtables(contimagename+ext)

        assert os.path.exists(cleanmask)

        tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
               startmodel=[preselfcal_imagename + ".model.tt0",
                           preselfcal_imagename + ".model.tt1"],
               deconvolver='mtmfs', interactive=False, gridder='standard',
               pbcor=True, savemodel='none', mask=cleanmask,
               **impars
              )

    caltable = field+'_b3_selfcal_phase{selfcaliter}_T.cal'.format(selfcaliter=selfcaliter)
    if not os.path.exists(caltable):
        impars['niter'] = 1
        tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
               deconvolver='mtmfs', interactive=False, gridder='standard',
               pbcor=True, savemodel='modelcolumn', mask=cleanmask,
               **impars
              )

        gaincal(vis=contvis, field=field, caltable=caltable,
                gaintable=cals, **calpars)


    cals += [caltable]

    applycal(vis=contvis,
             gaintable=cals,
             gainfield=[field]*len(cals),
             interp="linear",
             applymode='calonly',
             calwt=False)

