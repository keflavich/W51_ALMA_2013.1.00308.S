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


contimagename = 'w51_spw3_continuum_7m12m_noflag'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      mode='mfs',
      psfmode='clark',
      imsize = [960,960],
      cell= '0.15arcsec',
      weighting = 'natural',
      robust = 2.0,
      niter = 10000,
      threshold = '1.0mJy',
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
      imsize = [2560,2560],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 10000,
      threshold = '1.0mJy',
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
      imsize = [2560,2560],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 0,
      threshold = '1.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

contimagename = 'w51_spw3_continuum_7m12m_noflag_r0_mulstiscale'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      multiscale=[0,5,15,45],
      phasecenter='',
      mode='mfs',
      psfmode='clark',
      imsize = [2560,2560],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 10000,
      threshold = '10.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

