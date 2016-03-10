import time
import selfcal_heuristics
import numpy as np
t0 = time.time()

phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"

fields_to_selfcal = (31,32,33,39,40,24,25,20,13,21,27)
fields_after_split = [f-4 for f in fields_to_selfcal]
field = ",".join([str(x) for x in fields_after_split])

contvis='continuum_7m12m.ms'
vis0 = 'w51_continuum_7m12m_contvis_selfcal_0.ms'

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
summary_init = flagdata(vis=vis0, mode='summary')
print("{flagged}/{total} of flagged points in vis0".format(**summary_init))

imsize = [3072,3072]
cell = '0.05arcsec'
solint = 'int'
threshold = '20.0mJy'
multiscale = [0,5,15,45,135]
# multiscale clean... does not work without a mask (but it's not bad with a
# mask)
#multiscale = []

clearcal(vis=vis0)
#flagmanager(vis=vis0, versionname='flagdata_1', mode='restore')
myimagebase = "merge_selfcal_spw3_dirty"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='',
       specmode='mfs', outframe='LSRK', interpolation='linear', gridder='mosaic',
       interactive=False, niter=0, threshold=threshold, imsize=imsize,
       cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

dirtyimage = 'merge_selfcal_spw3_dirty.residual'
ia.open(dirtyimage)
ia.calcmask(mask=dirtyimage+" > 0.1", name='dirty_mask_100mJy')
ia.close()
makemask(mode='copy', inpimage=dirtyimage,
         inpmask=dirtyimage+":dirty_mask_100mJy", output='dirty_100mJy.mask',
         overwrite=True)
exportfits('dirty_100mJy.mask', 'dirty_100mJy.mask.fits', dropdeg=True, overwrite=True)


# rms ~3.7 mJy/beam
myimagebase = "merge_selfcal_spw3_mfs"
# threshold = 50mJy with no other restrictions -> infinite divergence
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       scales=multiscale, interactive=False, niter=10000,
       threshold='100mJy', imsize=imsize, specmode='mfs',
       mask='dirty_100mJy.mask', savemodel='modelcolumn',
       cell=cell, phasecenter=phasecenter, weighting='briggs', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables('selfcal_spw3_phase.cal')
gaincal(vis=vis0, caltable="selfcal_spw3_phase.cal", field=field, solint='inf',
        calmode="p", refant="", gaintype="G", minsnr=5, uvrange='100~5000m')

image1 = 'merge_selfcal_spw3_mfs.image'
ia.open(image1)
ia.calcmask(mask=image1+" > 0.05", name='clean_mask_50mJy')
ia.close()
cleanimage = myimagebase+".image"
makemask(mode='copy', inpimage=image1,
         inpmask=cleanimage+":clean_mask_50mJy", output='clean_50mJy.mask',
         overwrite=True)

flagmanager(vis=vis0, mode='save', versionname='backup')
applycal(vis=vis0, field="", gaintable=["selfcal_spw3_phase.cal"],
         interp="linear", applymode='calonly', calwt=False)
flagmanager(vis=vis0, mode='restore', versionname='backup')
summary0 = flagdata(vis=vis0, mode='summary')
print("{flagged}/{total} flagged points in vis0".format(**summary0))
vis1 = 'w51_continuum_7m12m_contvis_selfcal_1.ms'
os.system('rm -rf {0}'.format(vis1))
os.system('rm -rf {0}.flagversions'.format(vis1))
split(vis=vis0, outputvis=vis1,
      datacolumn="corrected")

myimagebase="selfcal_spw3_selfcal_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis1, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       scales=multiscale, interactive=False, niter=10000,
       mask='clean_50mJy.mask',
       threshold='50mJy', imsize=imsize, specmode='mfs',
       cell=cell, phasecenter=phasecenter, weighting='briggs', robust=-2.0,
       savemodel='modelcolumn')
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables('selfcal_spw3_phase_2.cal')
gaincal(vis=vis1, caltable="selfcal_spw3_phase_2.cal", field=field,
        solint=solint, calmode="p", refant="DV07", gaintype="G", minsnr=5,
        uvrange='100~5000m')
#plotcal(caltable="selfcal_spw3_phase_2.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)
assert len(selfcal_heuristics.goodenough_field_solutions("selfcal_spw3_phase_2.cal", minsnr=5, maxphasenoise=np.pi/4., pols=[0])) > 0
print("Goodenough field solns: ",selfcal_heuristics.goodenough_field_solutions("selfcal_spw3_phase_2.cal", minsnr=5, maxphasenoise=np.pi/4., pols=[0]))


flagmanager(vis=vis1, mode='save', versionname='backup')
applycal(vis=vis1, field="", gaintable=["selfcal_spw3_phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)
flagmanager(vis=vis1, mode='restore', versionname='backup')
summary1 = flagdata(vis=vis1, mode='summary')
print("{flagged}/{total} flagged points in vis1".format(**summary1))
vis2 = 'w51_continuum_7m12m_contvis_selfcal_2.ms'
os.system('rm -rf {0}'.format(vis2))
os.system('rm -rf {0}.flagversions'.format(vis2))
split(vis=vis1, outputvis=vis2,
      datacolumn="corrected")

myimagebase = "merge_selfcal_spw3_selfcal_2_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis2, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       scales=multiscale, interactive=False, niter=10000,
       threshold='50mJy', imsize=imsize, specmode='mfs',
       mask='clean_50mJy.mask', savemodel='modelcolumn',
       cell=cell, phasecenter=phasecenter, weighting='briggs', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables("selfcal_spw3_phase_3.cal")
gaincal(vis=vis2, caltable="selfcal_spw3_phase_3.cal", field=field,
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5,
        uvrange='100~5000m')
#plotcal(caltable="phase_3.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis=vis2, mode='save', versionname='backup')
applycal(vis=vis2, field="", gaintable=["selfcal_spw3_phase_3.cal"],
         interp="linear", applymode='calonly', calwt=False)
flagmanager(vis=vis2, mode='restore', versionname='backup')
summary2 = flagdata(vis=vis2, mode='summary')
print("{flagged}/{total} flagged points in vis2".format(**summary2))
vis3 = 'w51_continuum_7m12m_contvis_selfcal_3.ms'
os.system('rm -rf {0}'.format(vis3))
os.system('rm -rf {0}.flagversions'.format(vis3))
split(vis=vis2, outputvis=vis3,
      datacolumn="corrected")

myimagebase = "merge_selfcal_spw3_selfcal_3_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis3, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       scales=multiscale, interactive=False, niter=10000,
       threshold='50mJy', imsize=imsize, specmode='mfs',
       mask='clean_50mJy.mask', savemodel='modelcolumn',
       cell=cell, phasecenter=phasecenter, weighting='briggs', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables("selfcal_spw3_phase_4.cal")
gaincal(vis=vis3, caltable="selfcal_spw3_phase_4.cal", field=field,
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5,
        uvrange='100~5000m')

rmtables("selfcal_spw3_ampphase.cal")
gaincal(vis=vis3, caltable="selfcal_spw3_ampphase.cal", field=field,
        solint=solint, solnorm=True, calmode="ap", refant="", gaintype="G",
        minsnr=5, uvrange='100~5000m')

flagmanager(vis=vis3, mode='save', versionname='backup')
applycal(vis=vis3, field="", gaintable=["selfcal_spw3_phase_4.cal", 'selfcal_spw3_ampphase.cal'],
         interp="linear", applymode='calonly', calwt=False)
summary3 = flagdata(vis=vis3, mode='summary')
print("{flagged}/{total} flagged points in vis3 afer applycal".format(**summary3))
flagmanager(vis=vis3, mode='restore', versionname='backup')
summary3 = flagdata(vis=vis3, mode='summary')
print("{flagged}/{total} flagged points in vis3 after restoration".format(**summary3))
vis4 = 'w51_continuum_7m12m_contvis_selfcal_4.ms'
os.system('rm -rf {0}'.format(vis4))
os.system('rm -rf {0}.flagversions'.format(vis4))
split(vis=vis3, outputvis=vis4,
      datacolumn="corrected")
summary4 = flagdata(vis=vis4, mode='summary')
print("{flagged}/{total} flagged points in vis4 after restoration".format(**summary4))

myimagebase = "merge_selfcal_spw3_selfcal_4ampphase_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw='',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       scales=multiscale, interactive=False, niter=10000,
       threshold='50mJy', imsize=imsize, specmode='mfs',
       mask='clean_50mJy.mask', savemodel='modelcolumn',
       cell=cell, phasecenter=phasecenter, weighting='briggs', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


myimagebase = "merge_selfcal_spw3_selfcal_4ampphase_mfs_tclean"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='multiscale', gridder='mosaic', outframe='LSRK',
       scales=multiscale, pblimit=0.4, interpolation='linear',
       interactive=False, niter=10000,
       threshold='20mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       mask='clean_50mJy.mask',
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = "merge_selfcal_spw3_selfcal_4ampphase_mfs_tclean_deeper"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='15mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       mask='clean_50mJy.mask',
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)




# using vis2.corrected = vis3.data
myimagebase = "merge_selfcal_spw3_selfcal_3_mfs_deeper"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis2, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=50000,
       threshold='5mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       mask='clean_50mJy.mask',
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


import numpy as np
from astropy.io import fits
print("Stats (mfs):")
slc = slice(1007,1434), slice(1644,1900)
sigma, peak = (fits.getdata('merge_selfcal_spw3_dirty.image.fits')[slc].std(),     np.nanmax(fits.getdata('merge_selfcal_spw3_dirty.image.fits')))
print("dirty:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('merge_selfcal_spw3_mfs.image.pbcor.fits')[slc].std(),           np.nanmax(fits.getdata('merge_selfcal_spw3_mfs.image.pbcor.fits')))
print("clean:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('merge_selfcal_spw3_mfs.image.pbcor.fits')[slc].std(),   np.nanmax(fits.getdata('merge_selfcal_spw3_mfs.image.pbcor.fits')))
print("selfcal:           peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('merge_selfcal_spw3_selfcal_2_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('merge_selfcal_spw3_selfcal_2_mfs.image.pbcor.fits')))
print("selfcal2:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('merge_selfcal_spw3_selfcal_3_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('merge_selfcal_spw3_selfcal_3_mfs.image.pbcor.fits')))
print("selfcal3:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('merge_selfcal_spw3_selfcal_4ampphase_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('merge_selfcal_spw3_selfcal_4ampphase_mfs.image.pbcor.fits')))
print("selfcal4 ampphase: peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
print("Completed in {0}s".format(time.time()-t0))
