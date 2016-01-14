def clean(vis, imagename, **kwargs):
    tclean(vis = vis,
           imagename = imagename,
           mask='auto-pb', # masks at minpb=0.2.  0.4 or 0.5 are desired, but very difficult to configure
           **kwargs)
 #         field = '',
 #         spw = '', # there should be only one
 #         specmode = 'cube',
 #         width = width,
 #         start = startfreq,
 #         nchan = nchans_per_cube,
 #         veltype = 'radio',
 #         outframe = 'LSRK',
 #          gridder='mosaic',
 #          deconvolver='clark',
 #         interactive = F,
 #         niter = 25000,
 #         imsize = imsize,
 #         cell = cell,
 #         weighting = weighting,
 #         phasecenter = phasecenter,
 #         robust = robust,
 #         threshold = threshold,
 #         savemodel='none',
 #         overwrite=True)

"""
Attempt to image the continuum with NO flagging
"""

mergevis = 'continuum_7m12m_noflag.ms'
if not os.path.exists(mergevis):
    finalvis12m='calibrated_12m.ms'
    contvis12m='w51_spw3_continuum_12m.split'
    split(vis=finalvis12m,
          spw='3,7',
          outputvis=contvis12m,
          width=[192,192],
          datacolumn='data')

    finalvis7m='calibrated_7m.ms'
    contvis7m='w51_spw3_continuum_7m.split'
    split(vis=finalvis7m,
          spw='3',
          outputvis=contvis7m,
          width=[192],
          datacolumn='data')


    concat(vis=[contvis7m,contvis12m], concatvis=mergevis)

extensions = ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf',
              '.residual', '.flux.pbcoverage', '.sumwt', '.weight', '.pb',
              '.pbcoverage']

contimagename = 'w51_spw3_continuum_7m12m_noflag_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      specmode='mfs',
      deconvolver='clark',
      imsize = [1280,1280],
      cell= '0.15arcsec',
      weighting = 'natural',
      robust = 2.0,
      niter = 10000,
      threshold = '1.0mJy',
      interactive = False,
      gridder = 'mosaic',
      savemodel='none',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)


contimagename = 'w51_spw3_continuum_7m12m_noflag_r0_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      specmode='mfs',
      deconvolver='clark',
      imsize = [3072,3072],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 10000,
      threshold = '1.0mJy',
      interactive = False,
      gridder = 'mosaic',
      savemodel='none',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

contimagename = 'w51_spw3_continuum_7m12m_noflag_r0_dirty_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      specmode='mfs',
      deconvolver='clark',
      imsize = [3072,3072],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 0,
      threshold = '1.0mJy',
      interactive = False,
      gridder = 'mosaic',
      savemodel='none',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

contimagename = 'w51_spw3_continuum_7m12m_noflag_r0_multiscale_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      scales=[0,5,15,45],
      phasecenter='',
      specmode='mfs',
      deconvolver='multiscale',
      imsize = [3072,3072],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 10000,
      threshold = '10.0mJy',
      interactive = False,
      gridder = 'mosaic',
      savemodel='none',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

#contimagename = 'w51_spw3_continuum_7m12m_noflag_r0_MEM_tclean'
#
#for ext in extensions:
#    rmtables(contimagename+ext)
#
#clean(vis=mergevis,
#      imagename=contimagename,
#      field='w51',
#      scales=[0,3,6,9,12,15,18],
#      phasecenter='',
#      specmode='mfs',
#      deconvolver='mem',
#      imsize = [3072,3072],
#      cell= '0.052arcsec',
#      weighting = 'briggs',
#      robust = 0.0,
#      niter = 10000,
#      threshold = '10.0mJy',
#      interactive = False,
#      gridder = 'mosaic',
#      savemodel='none',
#      )
#exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)


contimagename = 'w51_spw3_continuum_7m12m_noflag_uniform_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      specmode='mfs',
      deconvolver='clark',
      imsize = [3072,3072],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = -2.0,
      niter = 50000,
      threshold = '20.0mJy',
      interactive = False,
      gridder = 'mosaic',
      savemodel='none',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

