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


imsize_lo = [1280,1280]
cell_lo = '0.15arcsec'
multiscale_lo = [0, 3, 9, 27, 81]

imsize_hi = [3072,3072]
cell_hi = '0.052arcsec'
multiscale_hi = [0, 5, 15, 45, 135]

contimagename = 'w51_spw3_continuum_7m12m_noflag'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      mode='mfs',
      psfmode='clark',
      imsize = imsize_lo,
      cell= cell_lo,
      weighting = 'natural',
      robust = 2.0,
      niter = 50000,
      threshold = '20.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)


contimagename = 'w51_spw3_continuum_7m12m_noflag_r0'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      mode='mfs',
      psfmode='clark',
      imsize = imsize_hi,
      cell= cell_hi,
      weighting = 'briggs',
      robust = 0.0,
      niter = 50000,
      threshold = '2.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

contimagename = 'w51_spw3_continuum_7m12m_noflag_r0_dirty'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      mode='mfs',
      psfmode='clark',
      imsize = imsize_hi,
      cell= cell_hi,
      weighting = 'briggs',
      robust = 0.0,
      niter = 0,
      threshold = '1.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

contimagename = 'w51_spw3_continuum_7m12m_noflag_r0_multiscale'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      multiscale=multiscale_hi,
      phasecenter='',
      mode='mfs',
      psfmode='clark',
      imsize = imsize_hi,
      cell= cell_hi,
      weighting = 'briggs',
      robust = 0.0,
      niter = 50000,
      threshold = '20.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

