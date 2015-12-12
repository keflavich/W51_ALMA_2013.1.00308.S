phasecenter = "J2000 19:23:41.585 +14:30:41.00"
phasecenter=''
def clean(vis, imagename, **kwargs):
    tclean(vis = vis,
           imagename = imagename,
           phasecenter = phasecenter,
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
 #         robust = robust,
 #         threshold = threshold,
 #         savemodel='none',
 #         overwrite=True)

"""
Attempt to image the continuum with NO flagging
"""

mergevis='w51_concat.ms.split.cal'

extensions = ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf',
              '.residual', '.flux.pbcoverage', '.sumwt', '.weight', '.pb',
              '.pbcoverage']

contimagename = 'w51_spw3_continuum_noflag_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      spw='3,7',
      specmode='mfs',
      deconvolver='clark',
      imsize = [1280,1280],
      cell= '0.15arcsec',
      weighting = 'natural',
      robust = 2.0,
      niter = 50000,
      threshold = '1.0mJy',
      interactive = False,
      gridder = 'mosaic',
      savemodel='none',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)


contimagename = 'w51_spw3_continuum_noflag_r0_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      spw='3,7',
      specmode='mfs',
      deconvolver='clark',
      imsize = [3072,3072],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 50000,
      threshold = '1.0mJy',
      interactive = False,
      gridder = 'mosaic',
      savemodel='none',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

contimagename = 'w51_spw3_continuum_noflag_r0_dirty_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      spw='3,7',
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

contimagename = 'w51_spw3_continuum_noflag_r0_multiscale_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      scales=[0,5,15,45],
      spw='3,7',
      specmode='mfs',
      deconvolver='multiscale',
      imsize = [3072,3072],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 50000,
      threshold = '10.0mJy',
      interactive = False,
      gridder = 'mosaic',
      savemodel='none',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

#contimagename = 'w51_spw3_continuum_noflag_r0_MEM_tclean'
#
#for ext in extensions:
#    rmtables(contimagename+ext)
#
#clean(vis=mergevis,
#      imagename=contimagename,
#      field='w51',
#      spw='3,7',
#      scales=[0,3,6,9,12,15,18],
#      specmode='mfs',
#      deconvolver='mem',
#      imsize = [3072,3072],
#      cell= '0.052arcsec',
#      weighting = 'briggs',
#      robust = 0.0,
#      niter = 50000,
#      threshold = '10.0mJy',
#      interactive = False,
#      gridder = 'mosaic',
#      savemodel='none',
#      )
#exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)


contimagename = 'w51_spw3_continuum_noflag_uniform_tclean'

for ext in extensions:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      spw='3,7',
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
