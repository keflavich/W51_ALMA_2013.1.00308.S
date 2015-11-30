def clean(vis, imagename, **kwargs):
    tclean(vis = vis,
           imagename = imagename,
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


mergevis = 'continuum_7m12m_noflag.ms'
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
      imsize = [960,960],
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
      imsize = [2560,2560],
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
      imsize = [2560,2560],
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

contimagename = 'w51_spw3_continuum_7m12m_noflag_r0_mulstiscale_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      scales=[0,3,6,9,12,15,18],
      phasecenter='',
      specmode='mfs',
      deconvolver='multiscale',
      imsize = [2560,2560],
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
