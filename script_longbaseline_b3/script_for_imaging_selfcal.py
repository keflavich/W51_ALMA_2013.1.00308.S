import os



def test_tclean_success():
    # An EXTREMELY HACKY way to test whether tclean succeeded on the previous iteration
    with open(casalog.logfile(), "r") as fh:
        lines = fh.readlines()

    for line in lines[-5:]:
        if 'SEVERE  tclean::::      An error occurred running task tclean.' in line:
            raise ValueError("tclean failed.  See log for detailed error report.\n{0}".format(line))
        if 'SEVERE' in line:
            raise ValueError("SEVERE error message encountered: {0}".format(line))

def check_model_is_populated(msfile):
    casalog.post("Checking for populated model column in "+msfile)
    ms.open(msfile)
    modelphase = ms.getdata(items=['model_phase'])
    if 'model_phase' not in modelphase:
        raise ValueError("model_phase not acquired")
    if modelphase['model_phase'].shape == (0,):
        raise ValueError("Model phase column was not populated")
    if np.all(modelphase['model_phase'] == 0):
        raise ValueError("Phase is zero.")
    ms.close()






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
                   'calibration': {'calmode': 'p', 'gaintype': 'T', 'solint': 'inf', 'solnorm': True},},
               2: {'imaging': {'threshold': '0.075mJy', 'nterms': 2, 'robust': 0,
                               'weighting':'briggs', 'cell':'0.007arcsec',
                               'imsize':14500, 'scales':[0,6,18],
                               'niter':200000},
                   'calibration': {'calmode': 'p', 'gaintype': 'T', 'solint': 'inf', 'solnorm': True},},
               3: {'imaging': {'threshold': '0.075mJy', 'nterms': 2, 'robust': 0,
                               'weighting':'briggs', 'cell':'0.007arcsec',
                               'imsize':14500, 'scales':[0,6,18],
                               'niter':200000},
                   'calibration': {'calmode': 'p', 'gaintype': 'G', 'solint': 'inf', 'solnorm': True},},
               4: {'imaging': {'threshold': '0.075mJy', 'nterms': 2, 'robust': 0,
                               'weighting':'briggs', 'cell':'0.007arcsec',
                               'imsize':14500, 'scales':[0,6,18],
                               'niter':200000},
                   'calibration': {'calmode': 'p', 'gaintype': 'T',
                                   'solint': 'int', 'solnorm': True},},
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
preselfcal_imagename = contimagename = field+'.spw0thru19.{imsize}.robust{robust}.thr{threshold}.mfs.I'.format(
    **impars)

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt',]:
        for extext in ('','.tt0','.tt1','.tt2'):
            if os.path.exists(contimagename+ext+extext):
                rmtables(contimagename+ext+extext)

    assert os.path.exists(cleanmask)

    niter = impars.pop('niter')

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs',
           interactive = False, gridder = 'standard', pbcor = True,
           savemodel='none', mask=cleanmask,
           niter=10,
           **impars
          )
    test_tclean_success()

    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs',
           interactive = False, gridder = 'standard', pbcor = True,
           savemodel='none',
           # mask=cleanmask, # mask is already made
           niter=niter,
           **impars
          )
    test_tclean_success()

caltable = field+'_b3_selfcal_phase{selfcaliter}_T.cal'.format(selfcaliter=selfcaliter)
if not os.path.exists(caltable):
    impars['niter'] = 0
    tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
           deconvolver='mtmfs', interactive=False, gridder='standard',
           calcres=False, calcpsf=False,
           pbcor=True, savemodel='modelcolumn', **impars)
    test_tclean_success()
    check_model_is_populated(contvis)

    gaincal(vis=contvis, field=field, caltable=caltable, **calpars)

cals = [caltable]

applycal(vis=contvis,
         gaintable=cals,
         gainfield=[field],
         interp="linear",
         applymode='calonly',
         calwt=False)



for selfcaliter in selfcalpars:

    if selfcaliter == 0:
        continue

    impars = selfcalpars[selfcaliter]['imaging']
    calpars = selfcalpars[selfcaliter]['calibration']

    prevcal_imagename = contimagename = field+'.spw0thru19.{imsize}.robust{robust}.thr{threshold}.mfs.I.selfcal{selfcaliter}'.format(
        selfcaliter=selfcaliter,
        **impars)

    if not os.path.exists(contimagename+'.image.tt0.pbcor'):
        for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
            for extext in ('','.tt0','.tt1','.tt2'):
                if os.path.exists(contimagename+ext+extext):
                    rmtables(contimagename+ext+extext)

        assert os.path.exists(cleanmask)

        tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
               startmodel=[preselfcal_imagename + ".model.tt0",
                           preselfcal_imagename + ".model.tt1"],
               deconvolver='mtmfs', interactive=False, gridder='standard',
               pbcor=True, savemodel='none', mask=cleanmask,
               **impars
              )
        test_tclean_success()

    caltable = field+'_b3_selfcal_phase{selfcaliter}_T.cal'.format(selfcaliter=selfcaliter)
    if not os.path.exists(caltable):
        impars['niter'] = 0
        tclean(vis=contvis, imagename=contimagename, field=field, specmode='mfs',
               deconvolver='mtmfs', interactive=False, gridder='standard',
               pbcor=True, savemodel='modelcolumn',
               calcres=False, calcpsf=False,
               **impars
              )
        test_tclean_success()
        check_model_is_populated(contvis)

        gaincal(vis=contvis, field=field, caltable=caltable,
                gaintable=cals, **calpars)


    cals += [caltable]

    # double-check that there are no repeated entries
    assert len(cals) == len(set(cals))

    applycal(vis=contvis,
             gaintable=cals,
             gainfield=[field]*len(cals),
             interp="linear",
             applymode='calonly',
             calwt=False)

