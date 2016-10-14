"""
Attempt to image the continuum by flagging out lines, splitting, then doing
imaging things...
"""

mergevis = 'continuum_7m12m.ms'
if not os.path.exists(mergevis):
    with open('linechannels12m_spw3','r') as f:
        linechannels12m = f.read()


    finalvis12m='calibrated_12m.ms'
    contvis12m='w51_spw3_continuum_flagged_12m.split'
    flagmanager(vis=finalvis12m,mode='save',
                versionname='before_cont_flags')


    flagdata(vis=finalvis12m,mode='manual',
             spw=linechannels12m,flagbackup=False)

    split(vis=finalvis12m,
          spw='3,7',
          outputvis=contvis12m,
          width=[192,192],
          datacolumn='data')


    flagmanager(vis=finalvis12m,mode='restore',
                versionname='before_cont_flags')

    with open('linechannels7m_spw3','r') as f:
        linechannels7m = f.read()
    finalvis7m='calibrated_7m.ms'
    contvis7m='w51_spw3_continuum_flagged_7m.split'
    flagmanager(vis=finalvis7m,mode='save',
                versionname='before_cont_flags')

    flagdata(vis=finalvis7m,mode='manual',
             spw=linechannels7m,flagbackup=False)

    split(vis=finalvis7m,
          spw='3',
          outputvis=contvis7m,
          width=[192],
          datacolumn='data')


    flagmanager(vis=finalvis7m,mode='restore',
                versionname='before_cont_flags')


    concat(vis=[contvis7m,contvis12m], concatvis=mergevis)


extnames = ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf',
            '.residual', '.flux.pbcoverage']

contimagename = 'w51_spw3_continuum_7m12m'

for ext in extnames:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='',
      mode='mfs',
      psfmode='clark',
      imsize = [1280,1280],
      cell= '0.15arcsec',
      weighting = 'natural',
      robust = 2.0,
      niter = 10000,
      threshold = '2.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)


contimagename = 'w51_spw3_continuum_7m12m_r0'

for ext in extnames:
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
      threshold = '2.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

contimagename = 'w51_spw3_continuum_7m12m_r0_dirty'

for ext in extnames:
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


dirtyimage = 'w51_spw3_continuum_7m12m_r0_dirty.image'
ia.open(dirtyimage)
ia.calcmask(mask=dirtyimage+" > 0.05", name='dirty_mask_50mJy')
ia.close()
makemask(mode='copy', inpimage=dirtyimage,
         inpmask=dirtyimage+":dirty_mask_50mJy", output='dirty_50mJy.mask')


# looks pretty good.  This is the best one I've seen yet.
contimagename = 'w51_spw3_continuum_7m12m_r0_multiscale_minpb'

for ext in extnames:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      multiscale=[0,5,15,45,135],
      phasecenter='4',
      mode='mfs',
      psfmode='clark',
      imsize = [2560,2560],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 50000,
      threshold = '20.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      pbcor=True,
      #mask='dirty_50mJy.mask',
      mask=0.5, # minpb
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

contimagename = 'w51_spw3_continuum_7m12m_r0_multiscale_mask50mJy'

for ext in extnames:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      multiscale=[0,5,15,45,135],
      phasecenter='4',
      mode='mfs',
      psfmode='clark',
      imsize = [2560,2560],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 50000,
      threshold = '20.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      pbcor=True,
      mask='dirty_50mJy.mask',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)


# results in gridded artifacts
contimagename = 'w51_spw3_continuum_7m12m_r0_monoscale_mask50mJy'

for ext in extnames:
    rmtables(contimagename+ext)

clean(vis=mergevis,
      imagename=contimagename,
      field='w51',
      phasecenter='4',
      mode='mfs',
      psfmode='clark',
      imsize = [2560,2560],
      cell= '0.052arcsec',
      weighting = 'briggs',
      robust = 0.0,
      niter = 50000,
      threshold = '1.0mJy',
      interactive = False,
      imagermode = 'mosaic',
      usescratch=False,
      pbcor=True,
      mask='dirty_50mJy.mask',
      )
exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)
