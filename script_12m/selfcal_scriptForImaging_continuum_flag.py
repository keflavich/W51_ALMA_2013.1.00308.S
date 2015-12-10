import time
t0 = time.time()

phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"

contvis='w51_spw3_continuum_flagged.split'
vis0 = 'w51_contvis_selfcal_0.ms'

os.system('rm -rf {0}'.format(vis0))
os.system('rm -rf {0}.flagversions'.format(vis0))
assert split(vis=contvis,
      outputvis=vis0,
      #field=','.join([str(x-4) for x in (31,32,33,39,40,24,25)]),
      field='w51', # 32-4
      spw='',
      datacolumn='data',
     )

print("Done splitting")

imsize = [3072,3072]
cell = '0.05arcsec'
solint = 'int'
threshold = '50.0mJy'
multiscale = [0,5,15,45]
#multiscale = []

clearcal(vis=vis0)
#flagmanager(vis=vis0, versionname='flagdata_1', mode='restore')
myimagebase = "selfcal_spw3_dirty"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis=vis0, imagename=myimagebase, field="", spw='',
      mode='mfs', outframe='LSRK', interpolation='linear', imagermode='mosaic',
      interactive=False, niter=0, threshold=threshold, imsize=imsize,
      cell=cell, phasecenter=phasecenter,
      minpb=0.4,
      weighting='briggs', usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = "selfcal_spw3_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis=vis0, imagename=myimagebase, field="", spw='',
      mode='mfs', outframe='LSRK', interpolation='linear', imagermode='mosaic',
      multiscale=multiscale,
      interactive=False, niter=10000, threshold=threshold, imsize=imsize,
      minpb=0.4,
      cell=cell, phasecenter=phasecenter,
      weighting='briggs', usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables('phase.cal')
gaincal(vis=vis0, caltable="phase.cal", field="", solint='inf',
        calmode="p", refant="", gaintype="G", minsnr=5)

#plotcal(caltable="phase.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis=vis0, mode='save', versionname='backup')
applycal(vis=vis0, field="", gaintable=["phase.cal"],
         interp="linear", applymode='calonly')
flagmanager(vis=vis0, mode='restore', versionname='backup')
vis1 = 'w51_contvis_selfcal_1.ms'
os.system('rm -rf {0}'.format(vis1))
os.system('rm -rf {0}.flagversions'.format(vis1))
split(vis=vis0, outputvis=vis1,
      datacolumn="corrected")

myimagebase="selfcal_spw3_selfcal_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis=vis1, imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=imsize, cell=cell,
      phasecenter=phasecenter, weighting='briggs',
      minpb=0.4,
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables('phase_2.cal')
gaincal(vis=vis1, caltable="phase_2.cal", field="",
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5)
#plotcal(caltable="phase_2.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)


flagmanager(vis=vis1, mode='save', versionname='backup')
applycal(vis=vis1, field="", gaintable=["phase_2.cal"],
         interp="linear", applymode='calonly')
flagmanager(vis=vis1, mode='restore', versionname='backup')
vis2 = 'w51_contvis_selfcal_2.ms'
os.system('rm -rf {0}'.format(vis2))
os.system('rm -rf {0}.flagversions'.format(vis2))
split(vis=vis1, outputvis=vis2,
      datacolumn="corrected")

myimagebase = "selfcal_spw3_selfcal_2_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis=vis2, imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=imsize, cell=cell,
      minpb=0.4,
      phasecenter=phasecenter, weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables("phase_3.cal")
gaincal(vis=vis2, caltable="phase_3.cal", field="",
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5)
#plotcal(caltable="phase_3.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis=vis2, mode='save', versionname='backup')
applycal(vis=vis2, field="", gaintable=["phase_3.cal"],
         interp="linear", applymode='calonly')
flagmanager(vis=vis2, mode='restore', versionname='backup')
vis3 = 'w51_contvis_selfcal_3.ms'
os.system('rm -rf {0}'.format(vis3))
os.system('rm -rf {0}.flagversions'.format(vis3))
split(vis=vis2, outputvis=vis3,
      datacolumn="corrected")

myimagebase = "selfcal_spw3_selfcal_3_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis=vis3, imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      interpolation='linear', imagermode='mosaic', interactive=False,
      minpb=0.4,
      niter=10000, threshold=threshold, imsize=imsize, cell=cell,
      phasecenter=phasecenter, weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables("phase_4.cal")
gaincal(vis=vis3, caltable="phase_4.cal", field="",
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5)
#plotcal(caltable="phase_4.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)


rmtables("ampphase.cal")
gaincal(vis=vis3, caltable="ampphase.cal", field="",
        solint=solint, solnorm=True, calmode="ap", refant="", gaintype="G", minsnr=5)

flagmanager(vis=vis3, mode='save', versionname='backup')
applycal(vis=vis3, field="", gaintable=["phase_4.cal", 'ampphase.cal'],
         interp="linear", applymode='calonly')
flagmanager(vis=vis3, mode='restore', versionname='backup')
vis4 = 'w51_contvis_selfcal_4.ms'
os.system('rm -rf {0}'.format(vis4))
os.system('rm -rf {0}.flagversions'.format(vis4))
split(vis=vis3, outputvis=vis4,
      datacolumn="corrected")

myimagebase = "selfcal_spw3_selfcal_4ampphase_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
clean(vis=vis4, imagename=myimagebase,
      field="", spw='', mode='mfs', outframe='LSRK',
      multiscale=multiscale,
      minpb=0.4,
      interpolation='linear', imagermode='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=imsize, cell=cell,
      phasecenter=phasecenter, weighting='briggs',
      usescratch=True, pbcor=True, robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)



import numpy as np
from astropy.io import fits
print("Stats (mfs):")
slc = slice(470,612), slice(740,855)
sigma, peak = (fits.getdata('selfcal_spw3_mfs_dirty.image.fits')[slc].std(),     np.nanmax(fits.getdata('selfcal_spw3_mfs_dirty.image.fits')))
print("dirty:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_spw3_mfs.image.pbcor.fits')[slc].std(),           np.nanmax(fits.getdata('selfcal_spw3_mfs.image.pbcor.fits')))
print("clean:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_spw3_selfcal_mfs.image.pbcor.fits')[slc].std(),   np.nanmax(fits.getdata('selfcal_spw3_selfcal_mfs.image.pbcor.fits')))
print("selfcal:           peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_spw3_selfcal_2_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('selfcal_spw3_selfcal_2_mfs.image.pbcor.fits')))
print("selfcal2:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_spw3_selfcal_3_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('selfcal_spw3_selfcal_3_mfs.image.pbcor.fits')))
print("selfcal3:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_spw3_selfcal_4ampphase_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('selfcal_spw3_selfcal_4ampphase_mfs.image.pbcor.fits')))
print("selfcal4 ampphase: peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
print("Completed in {0}s".format(time.time()-t0))
