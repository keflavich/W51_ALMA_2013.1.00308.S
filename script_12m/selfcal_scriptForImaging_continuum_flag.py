"""
run this then iterative_clean to get the "best" image
"""
import time
t0 = time.time()

phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"

fields_to_selfcal = (31,32,33,39,40,24,25,20,13,21,27)
fields_after_split = [f-4 for f in fields_to_selfcal]
field = ",".join([str(x) for x in fields_after_split])

contvis=['w51_spw{0}_continuum_flagged.split'.format(ii) for ii in range(3)]
vis0 = 'w51_contvis_selfcal_0.ms'

os.system('rm -rf {0}'.format(vis0))
os.system('rm -rf {0}.flagversions'.format(vis0))
assert concat(vis=contvis, concatvis=vis0,)

print("Done splitting")
summary_init = flagdata(vis=vis0, mode='summary')
print("{flagged}/{total} of flagged points in vis0".format(**summary_init))

imsize = [3072,3072]
cell = '0.05arcsec'
solint = 'int'
threshold = '50.0mJy'

clearcal(vis=vis0)
#flagmanager(vis=vis0, versionname='flagdata_1', mode='restore')
myimagebase = "selfcal_allspw_dirty"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='',
       specmode='mfs', outframe='LSRK', interpolation='linear', gridder='mosaic',
       interactive=False, niter=0, threshold=threshold, imsize=imsize,
       cell=cell, phasecenter=phasecenter,
       pblimit=0.4,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = "selfcal_allspw_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='',
      specmode='mfs', outframe='LSRK', interpolation='linear', gridder='mosaic',
      interactive=False, niter=10000, threshold=threshold, imsize=imsize,
      pblimit=0.4,
      cell=cell, phasecenter=phasecenter,
      weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables('selfcal_allspw_phase.cal')
gaincal(vis=vis0, caltable="selfcal_allspw_phase.cal", field=field, solint='inf',
        calmode="p", refant="", gaintype="G", minsnr=5)

#plotcal(caltable="phase.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis=vis0, mode='save', versionname='backup')
applycal(vis=vis0, field="", gaintable=["selfcal_allspw_phase.cal"],
         interp="linear", applymode='calonly', calwt=False)
flagmanager(vis=vis0, mode='restore', versionname='backup')
summary0 = flagdata(vis=vis0, mode='summary')
print("{flagged}/{total} flagged points in vis0".format(**summary0))
vis1 = 'w51_contvis_selfcal_1.ms'
os.system('rm -rf {0}'.format(vis1))
os.system('rm -rf {0}.flagversions'.format(vis1))
split(vis=vis0, outputvis=vis1,
      datacolumn="corrected")

myimagebase="selfcal_allspw_selfcal_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis1, imagename=myimagebase,
      field="", spw='', specmode='mfs', outframe='LSRK',
      
      interpolation='linear', gridder='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=imsize, cell=cell,
      phasecenter=phasecenter, weighting='briggs',
      pblimit=0.4,
      savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables('selfcal_allspw_phase_2.cal')
gaincal(vis=vis1, caltable="selfcal_allspw_phase_2.cal", field=field,
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5)
#plotcal(caltable="phase_2.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)


flagmanager(vis=vis1, mode='save', versionname='backup')
applycal(vis=vis1, field="", gaintable=["selfcal_allspw_phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)
flagmanager(vis=vis1, mode='restore', versionname='backup')
summary1 = flagdata(vis=vis1, mode='summary')
print("{flagged}/{total} flagged points in vis1".format(**summary1))
vis2 = 'w51_contvis_selfcal_2.ms'
os.system('rm -rf {0}'.format(vis2))
os.system('rm -rf {0}.flagversions'.format(vis2))
split(vis=vis1, outputvis=vis2,
      datacolumn="corrected")

myimagebase = "selfcal_allspw_selfcal_2_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis2, imagename=myimagebase,
      field="", spw='', specmode='mfs', outframe='LSRK',
      
      interpolation='linear', gridder='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=imsize, cell=cell,
      pblimit=0.4,
      phasecenter=phasecenter, weighting='briggs',
      savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables("selfcal_allspw_phase_3.cal")
gaincal(vis=vis2, caltable="selfcal_allspw_phase_3.cal", field=field,
        solint=solint, calmode="p", refant="", gaintype="G", minsnr=5)
#plotcal(caltable="phase_3.cal", xaxis="time", yaxis="phase", subplot=331,
#        iteration="antenna", plotrange=[0,0,-30,30], markersize=5,
#        fontsize=10.0,)

flagmanager(vis=vis2, mode='save', versionname='backup')
applycal(vis=vis2, field="", gaintable=["selfcal_allspw_phase_3.cal"],
         interp="linear", applymode='calonly', calwt=False)
flagmanager(vis=vis2, mode='restore', versionname='backup')
summary2 = flagdata(vis=vis2, mode='summary')
print("{flagged}/{total} flagged points in vis2".format(**summary2))
vis3 = 'w51_contvis_selfcal_3.ms'
os.system('rm -rf {0}'.format(vis3))
os.system('rm -rf {0}.flagversions'.format(vis3))
split(vis=vis2, outputvis=vis3,
      datacolumn="corrected")

myimagebase = "selfcal_allspw_selfcal_3_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis3, imagename=myimagebase,
      field="", spw='', specmode='mfs', outframe='LSRK',
      
      interpolation='linear', gridder='mosaic', interactive=False,
      pblimit=0.4,
      niter=10000, threshold=threshold, imsize=imsize, cell=cell,
      phasecenter=phasecenter, weighting='briggs',
      savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

rmtables("selfcal_allspw_phase_4.cal")
gaincal(vis=vis3, caltable="selfcal_allspw_phase_4.cal", field=field, solint=solint, calmode="p",
        refant="", gaintype="G", minsnr=5)

rmtables("selfcal_allspw_ampphase.cal")
gaincal(vis=vis3, caltable="selfcal_allspw_ampphase.cal", field=field, solint=solint,
        solnorm=True, calmode="ap", refant="", gaintype="G", minsnr=5)

flagmanager(vis=vis3, mode='save', versionname='backup')
applycal(vis=vis3, field="", gaintable=["selfcal_allspw_phase_4.cal", 'selfcal_allspw_ampphase.cal'],
         interp="linear", applymode='calonly', calwt=False)
summary3 = flagdata(vis=vis3, mode='summary')
print("{flagged}/{total} flagged points in vis3 afer applycal".format(**summary3))
flagmanager(vis=vis3, mode='restore', versionname='backup')
summary3 = flagdata(vis=vis3, mode='summary')
print("{flagged}/{total} flagged points in vis3 after restoration".format(**summary3))
vis4 = 'w51_contvis_selfcal_4.ms'
os.system('rm -rf {0}'.format(vis4))
os.system('rm -rf {0}.flagversions'.format(vis4))
split(vis=vis3, outputvis=vis4,
      datacolumn="corrected")
summary4 = flagdata(vis=vis4, mode='summary')
print("{flagged}/{total} flagged points in vis4 after restoration".format(**summary4))

# WHY IS THIS EMPTY?!
myimagebase = "selfcal_allspw_selfcal_4ampphase_mfs"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase,
      field="", spw='', specmode='mfs', outframe='LSRK',
      
      pblimit=0.4,
      interpolation='linear', gridder='mosaic', interactive=False,
      niter=10000, threshold=threshold, imsize=imsize, cell=cell,
      phasecenter=phasecenter, weighting='briggs',
      savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


myimagebase = "selfcal_allspw_selfcal_4ampphase_mfs_tclean"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=10000,
       threshold=threshold, imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

clearcal(vis=vis0)
myimagebase = "iter0_selfcal_allspw_mfs_deeper_5mjy"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis0, imagename=myimagebase, field="", spw='', specmode='mfs',
       outframe='LSRK', interpolation='linear', gridder='mosaic',
       interactive=False, niter=100000,
       threshold='5mJy', imsize=imsize, pblimit=0.4, cell=cell,
       phasecenter=phasecenter, weighting='briggs', savemodel='modelcolumn',
       robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = "selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='15mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = "selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_5mJy"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='5mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


myimagebase = "selfcal_allspw_selfcal_4ampphase_mfs_tclean_deeper_4mJy"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase, field="", spw="", specmode='mfs',
       deconvolver='clark', gridder='mosaic', outframe='LSRK',
       pblimit=0.4, interpolation='linear',
       interactive=False, niter=100000,
       threshold='4mJy', imsize=imsize, cell=cell, phasecenter=phasecenter,
       weighting='briggs', savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


# since ampphase fails by flagging out good data, try deeper here...
# No, something else has caused problems.  Screw it, try ampphase.
# oooooh, vis3.corrected = vis4.data, and clean automatically selected corrected... sigh
myimagebase = "selfcal_allspw_selfcal_4ampphase_mfs_deeper"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis4, imagename=myimagebase,
      field="", spw='', specmode='mfs', outframe='LSRK',
      
      interpolation='linear', gridder='mosaic', interactive=False,
      pblimit=0.4,
      niter=50000, threshold='5mJy', imsize=imsize, cell=cell,
      phasecenter=phasecenter, weighting='briggs',
      savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

# using vis2.corrected = vis3.data
myimagebase = "selfcal_allspw_selfcal_3_mfs_deeper"
os.system('rm -rf {0}.*'.format(myimagebase))
tclean(vis=vis2, imagename=myimagebase,
      field="", spw='', specmode='mfs', outframe='LSRK',
      
      interpolation='linear', gridder='mosaic', interactive=False,
      pblimit=0.4,
      niter=50000, threshold='5mJy', imsize=imsize, cell=cell,
      phasecenter=phasecenter, weighting='briggs',
      savemodel='modelcolumn', robust=-2.0)
exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


import numpy as np
from astropy.io import fits
print("Stats (mfs):")
slc = slice(1007,1434), slice(1644,1900)
sigma, peak = (fits.getdata('selfcal_allspw_dirty.image.fits')[slc].std(),     np.nanmax(fits.getdata('selfcal_allspw_dirty.image.fits')))
print("dirty:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_allspw_mfs.image.pbcor.fits')[slc].std(),           np.nanmax(fits.getdata('selfcal_allspw_mfs.image.pbcor.fits')))
print("clean:             peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_allspw_selfcal_mfs.image.pbcor.fits')[slc].std(),   np.nanmax(fits.getdata('selfcal_allspw_selfcal_mfs.image.pbcor.fits')))
print("selfcal:           peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_allspw_selfcal_2_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('selfcal_allspw_selfcal_2_mfs.image.pbcor.fits')))
print("selfcal2:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_allspw_selfcal_3_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('selfcal_allspw_selfcal_3_mfs.image.pbcor.fits')))
print("selfcal3:          peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
sigma, peak = (fits.getdata('selfcal_allspw_selfcal_4ampphase_mfs.image.pbcor.fits')[slc].std(), np.nanmax(fits.getdata('selfcal_allspw_selfcal_4ampphase_mfs.image.pbcor.fits')))
print("selfcal4 ampphase: peak={1:0.5f} sigma={0:0.5f} s/n={2:0.5f}".format(sigma, peak, peak/sigma))
print("Completed in {0}s".format(time.time()-t0))

"""
Stats (mfs):
dirty:             peak=0.57990 sigma=0.00056 s/n=1031.64363
tclean:             peak=0.62937 sigma=0.00059 s/n=1060.74211
selfcal:           peak=0.64963 sigma=0.00059 s/n=1104.06681
selfcal2:          peak=0.65766 sigma=0.00059 s/n=1115.24844
selfcal3:          peak=0.66613 sigma=0.00059 s/n=1129.29259
selfcal4 ampphase: peak=0.70497 sigma=0.00063 s/n=1126.32878
Completed in 2088.65619898s
"""


# # Extras to try to reproduce bug:
# phasecenter = "J2000 19:23:41.629000 +14.30.42.38000"
# imsize = [3072,3072]
# cell = '0.05arcsec'
# solint = 'int'
# threshold = '50.0mJy'
# multiscale = [0,5,15,45]
# vis4 = 'w51_contvis_selfcal_4.ms'
# myimagebase = "selfcal_allspw_selfcal_4ampphase_mfs_dirty"
# os.system('rm -rf {0}.*'.format(myimagebase))
# tclean(vis=vis4, imagename=myimagebase,
#       field="", spw='', specmode='mfs', outframe='LSRK',
#       
#       interpolation='linear', gridder='mosaic', interactive=False,
#       pblimit=0.4,
#       niter=0, threshold='5mJy', imsize=imsize, cell=cell,
#       phasecenter=phasecenter, weighting='briggs',
#       savemodel='modelcolumn', robust=-2.0)
# exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
# impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
#         outfile=myimagebase+'.image.pbcor', overwrite=True)
# exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
# exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
# exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)
# 
# 
# myimagebase = "selfcal_allspw_selfcal_4ampphase_mfs_dirty_field32"
# os.system('rm -rf {0}.*'.format(myimagebase))
# tclean(vis=vis4, imagename=myimagebase,
#       field="28", spw='', specmode='mfs', outframe='LSRK',
#       
#       interpolation='linear', gridder='mosaic', interactive=False,
#       pblimit=0.4,
#       niter=0, threshold='5mJy', imsize=imsize, cell=cell,
#       phasecenter="", weighting='briggs',
#       savemodel='modelcolumn', robust=-2.0)
# exportfits(myimagebase+'.image', myimagebase+'.image.fits', dropdeg=True, overwrite=True)
# impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
#         outfile=myimagebase+'.image.pbcor', overwrite=True)
# exportfits(myimagebase+'.image.pbcor', myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
# exportfits(myimagebase+'.model', myimagebase+'.model.fits', dropdeg=True, overwrite=True)
# exportfits(myimagebase+'.residual', myimagebase+'.residual.fits', dropdeg=True, overwrite=True)
