"""
Attempt to image the continuum by flagging out lines, splitting, then doing imaging things...
"""

with open('linechannels12m_spw3','r') as f:
    linechannels = f.read()


finalvis='w51_concat.ms.split.cal'
contvis='w51_spw3_continuum_flagged.split'
flagmanager(vis=finalvis,mode='save',
            versionname='before_cont_flags')

flagdata(vis=finalvis,mode='manual',
         spw=linechannels,flagbackup=False)

split(vis=finalvis,
      spw='3,7',
      field='w51',
      outputvis=contvis,
      width=[192,192],
      datacolumn='data')


flagmanager(vis=finalvis,mode='restore',
            versionname='before_cont_flags')

split(vis=finalvis,
      spw='3,7',
      field='w51',
      outputvis='w51_spw3_continuum_noflag.split',
      width=[192,192],
      datacolumn='data')


contimagename = 'w51_spw3_continuum'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=contvis,
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


contimagename = 'w51_spw3_continuum_r0'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=contvis,
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

contimagename = 'w51_spw3_continuum_r0_dirty'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=contvis,
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

contimagename = 'w51_spw3_continuum_r0_mulstiscale'

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
    rmtables(contimagename+ext)

clean(vis=contvis,
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
