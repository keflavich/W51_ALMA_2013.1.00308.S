import os

weighting = 'briggs'
robust=0.0
niter=50000
threshold = '20.0mJy'
phasecenter='4'

imsize_lo = [1280,1280]
cell_lo = '0.15arcsec'
multiscale_lo = [0, 3, 9, 27, 81]

imsize_hi = [3072,3072]
cell_hi = '0.052arcsec'
multiscale_hi = [0, 5, 15, 45, 135]

extnames = ['.flux', '.image', '.mask', '.model', '.pbcor', '.psf',
            '.residual', '.flux.pbcoverage']

spws = {0: '0,4',
        1: '1,5',
        2: '2,6',
        3: '3,7',
       }
finalvis='w51_concat.ms.split.cal'
# don't know how to contsub yet
mergevis = linevis = finalvis#+'.contsub'

for spwnum in '3210':
    print "# running clean on all lines in spw{0}".format(spwnum)
    spw = spws[int(spwnum)]

    contimagename = 'w51_spw{0}_continuum_noflag'.format(spwnum)

    for ext in extnames:
        rmtables(contimagename+ext)

    clean(vis=mergevis,
          imagename=contimagename,
          field='w51',
          phasecenter=phasecenter,
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
          spw=spw,
          )
    exportfits(contimagename+".image", contimagename+".image.fits",
               dropdeg=True, overwrite=True)


    contimagename = 'w51_spw{0}_continuum_noflag_r0'.format(spwnum)

    for ext in extnames:
        rmtables(contimagename+ext)

    clean(vis=mergevis,
          imagename=contimagename,
          field='w51',
          phasecenter=phasecenter,
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
          spw=spw,
          )
    exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

    contimagename = 'w51_spw{0}_continuum_noflag_r0_dirty'.format(spwnum)

    for ext in extnames:
        rmtables(contimagename+ext)

    clean(vis=mergevis,
          imagename=contimagename,
          field='w51',
          phasecenter=phasecenter,
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
          spw=spw,
          )
    exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

    contimagename = 'w51_spw{0}_continuum_noflag_r0_multiscale'.format(spwnum)

    for ext in extnames:
        rmtables(contimagename+ext)

    clean(vis=mergevis,
          imagename=contimagename,
          field='w51',
          multiscale=multiscale_hi,
          phasecenter=phasecenter,
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
          spw=spw,
          )
    exportfits(contimagename+".image", contimagename+".image.fits", dropdeg=True, overwrite=True)

