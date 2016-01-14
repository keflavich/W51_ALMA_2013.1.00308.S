
if not os.path.exists('w51_merge_test_small.ms'):
    assert split(vis='continuum_7m12m.ms',
                 outputvis='w51_merge_test_small.ms',
          #field=','.join([str(x-4) for x in (31,32,33,39,40,24,25)]),
          #field='28', # 32-4
          #field='27,35,36,28,55,52,56', # 32 I know, but 56...? must be 7m!
                 field='31,32,33,39,40,24,25,51,52,55,56',
                 spw='', # I guess we're using all spws?
                 datacolumn='data',
                )

for cleanmode in ('clark', 'hogbom', 'mtmfs', 'mem', 'multiscale'):
    imagename = 'tclean_test_mfs_uniform_{0}'.format(cleanmode)
    clearcal('w51_merge_test_small.ms')
    tclean(vis='w51_merge_test_small.ms',
          imagename=imagename,
          field='w51',
          spw='',
          specmode='mfs',
          deconvolver=cleanmode,
          imsize = [512,512],
          cell= '0.06arcsec',
          weighting = 'uniform',
          phasecenter='J2000 19h23m43.905 +14d30m28.08',
          scales=[0,3,9,27,81],
          robust = -2.0,
          niter = 50000,
          threshold = '1.0mJy',
          interactive = False,
          gridder = 'mosaic',
          savemodel='none',
          )
    if 'mtmfs' in cleanmode:
        imagetypes = ('image.tt0', 'model.tt0', 'residual.tt0')
    else:
        imagetypes = ('image', 'model', 'residual')
    for imagetype in imagetypes:
        exportfits(imagename+".{0}".format(imagetype),
                   imagename+".{0}.fits".format(imagetype), dropdeg=True,
                   overwrite=True)

for cleanmode in ('csclean', 'mosaic', 'multiscale', 'mutiscale_csclean'):
    imagename = 'clean_test_mfs_uniform_{0}'.format(cleanmode)
    if 'multiscale' in cleanmode:
        multiscale = [0,3,9,27,81]
    else:
        multiscale = []
    if 'csclean' in cleanmode:
        imagermode = 'csclean'
    else:
        imagermode = 'mosaic'

    clearcal('w51_merge_test_small.ms')
    clean(vis='w51_merge_test_small.ms',
          imagename=imagename,
          field='w51',
          spw='',
          mode='mfs',
          imsize = [512,512],
          cell= '0.06arcsec',
          weighting = 'uniform',
          phasecenter='J2000 19h23m43.905 +14d30m28.08',
          multiscale=multiscale,
          robust = -2.0,
          niter = 50000,
          threshold = '1.0mJy',
          interactive = False,
          imagermode = imagermode,
          savemodel='none',
          )
    for imagetype in ('image', 'model', 'residual'):
        exportfits(imagename+".{0}".format(imagetype),
                   imagename+".{0}.fits".format(imagetype), dropdeg=True,
                   overwrite=True)
