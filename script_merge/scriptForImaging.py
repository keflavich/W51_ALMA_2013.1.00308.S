# these are made in scriptForImaging_fullcube
# w51_concat_7m12m.spw0.merge/  w51_concat_7m12m.spw1.merge/  w51_concat_7m12m.spw2.merge/  w51_concat_7m12m.spw3.merge/


#Wavelength ~ 0.137 mm
#primary beam ~ 28.3 arcsec, however mode=mosaic with roughly 140 by 140 arcsec dimensions
#max uvdistance is 1500 m (Xb4b EB), however the majority of the data has shorter uvdistances
#cell size 0.15 arcsec is enough to sample the beam of 0.5 arcsec, thus imagesize is [1280,1280]
#No clean mask is used since the emission is complex

###################################
#w51 LSR velocity= 65.4 km/s, hel
##################################

#H2CO(3(0,3)-2(0,2)- Rest 218.222GHz - SPW 0

#H2CO(3(0,3)-2(0,2) -  Xb4b only, high res data set, cont. subtracted

fitspw='0:0~100;450~700;1250~1550;2500~2700;3000~3350'
uvcontsub(vis='w51_concat_7m12m.spw0.merge',
         field='',
         spw='', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=True)

linesub='w51_concat_7m12m.spw0.merge.contsub' #result of uvcontsub contains only the selected uvdata by the task

os.system("rm -rf w51_merge7m12m_H2CO_303_202_contsub_uniform.*")
clean(vis= linesub,
 imagename = "w51_merge7m12m_H2CO_303_202_contsub_uniform",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 70,
 start = '20km/s',
 width = '1.0km/s',
 restfreq = '218.222GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 15000,
 # go to >5-sigma, especially for uniform...
 threshold = '10mJy', #req rms 5.85 mJy, 38 antennas, 33.8min tos, pwv 2 (for EB Xb4b),0.2 arcsec res, 0.728 MHz BW gives 1.8mJy!
 imsize = [3840,3840], # make it much bigger to avoid edge effects
 cell = '0.052arcsec',  #synth beam expected to be 0.2 arcsec, so 0.2/3= 0.06 arcsec cell
 weighting = 'uniform',
 minpb=0.4,
 pbcor=False,
 robust=-2.0)


myimagebase= 'w51_merge7m12m_H2CO_303_202_contsub_uniform'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


myimagebase = "w51_merge7m12m_HC3N_24_23_contsub"
os.system("rm -rf {0}.*".format(myimagebase))
tclean(vis = linesub,
       imagename = myimagebase,
       field = "w51",
       spw = '0,4',
       specmode = 'cube',
       nchan = 140,
       start = '20km/s',
       width = '0.5km/s',
       restfreq = '218.32471GHz',
       outframe = 'LSRK',
       interpolation = 'linear',
       gridder='mosaic',
       pblimit=0.4,
       interactive = False,
       niter = 5000,
       threshold = '10mJy',
       imsize = [2880,2880],
       cell = '0.05arcsec',
       weighting = 'robust',
       robust = -0.5)
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb',
        outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',
           fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits') # export the PB image



os.system("rm -rf w51_merge7m12m_H2CO_303_202_contsub_briggs0.*")
clean(vis= linesub,
 imagename = "w51_merge7m12m_H2CO_303_202_contsub_briggs0",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 70,
 start = '20km/s',
 width = '1.0km/s',
 restfreq = '218.222GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 15000,
 # go to >5-sigma, especially for uniform...
 threshold = '10mJy', #req rms 5.85 mJy, 38 antennas, 33.8min tos, pwv 2 (for EB Xb4b),0.2 arcsec res, 0.728 MHz BW gives 1.8mJy!
 imsize = [3840,3840], # make it much bigger to avoid edge effects
 cell = '0.052arcsec',  #synth beam expected to be 0.2 arcsec, so 0.2/3= 0.06 arcsec cell
 weighting = 'briggs',
 minpb=0.4,
 pbcor=False,
 robust=0.0)


myimagebase= 'w51_merge7m12m_H2CO_303_202_contsub_briggs0'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


#H2CO(3(0,3)-2(0,2) - concatenated datasets

os.system("rm -rf w51_H2CO_303_202_merge7m12m_nocontsub.*")
clean(vis= 'w51_concat_7m12m.spw0.merge',
 imagename = 'w51_H2CO_303_202_merge7m12m_nocontsub',
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 70,
 start = '20km/s',
 width = '1.0km/s',
 restfreq = '218.222GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.5mJy', #req rms 5.85 mJy, 38 antennas, 66min tos, pwv auto (combination of two EB, which both had different PWV),1arcsec res, 0.728 MHz BW gives 1.5mJy!
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)


myimagebase= 'w51_H2CO_303_202_merge7m12m_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

os.system("rm -rf w51_H2CO_303_202_merge7m12m_contsub.*")
clean(vis= linesub,
 imagename = "w51_H2CO_303_202_merge7m12m_contsub",
 multiscale = [0,3,9,27,81,243],
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 70,
 start = '20km/s',
 width = '1.0km/s',
 restfreq = '218.222GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.5mJy', #req threshold 5.85 mJy, 38 antennas, 66min tos, pwv auto,1arcsec res, 0.728 MHz BW gives 1.5mJy!
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)

myimagebase= 'w51_H2CO_303_202_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

os.system("rm -rf w51_H2CO_303_202_merge7m12m_contsub_vresolved.*")
clean(vis= linesub,
 imagename = "w51_H2CO_303_202_merge7m12m_contsub_vresolved",
 multiscale = [0,3,9,27,81,243],
 field = "w51",
 spw = '0,1',
 mode = 'velocity',
 nchan = 350,
 start = '20km/s',
 width = '0.2km/s',
 restfreq = '218.222GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 10000,
 threshold = '10mJy',
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)


myimagebase= 'w51_H2CO_303_202_merge7m12m_contsub_vresolved'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


myimagebase= 'w51_H2CO_322_221_merge7m12m_contsub'
os.system("rm -rf {0}.*".format(myimagebase))
clean(vis= linesub,
 imagename = myimagebase,
 multiscale = [0,3,9,27,81,243],
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 70,
 start = '20km/s',
 width = '1.0km/s',
 restfreq = '218.47563GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.5mJy', #req threshold 5.85 mJy, 38 antennas, 66min tos, pwv auto,1arcsec res, 0.728 MHz BW gives 1.5mJy!
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)

impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image



myimagebase= 'w51_H2CO_322_221_merge7m12m_contsub_briggs0'
os.system("rm -rf {0}.*".format(myimagebase))
clean(vis= linesub,
 imagename = myimagebase,
 multiscale = [0,3,9,27,81,243], # multiscale has not yet been tuned to appropriate beam size
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 70,
 start = '20km/s',
 width = '1.0km/s',
 restfreq = '218.47563GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 15000,
 threshold = '10mJy', # higher noise with uniform... but I want 5-sig here
 imsize = [3840,3840],
 cell = '0.052arcsec',
 weighting = 'briggs',
 minpb=0.4,
 pbcor=False,
 robust=0.0)

impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image




#C18O(2-1)- Rest 219.560 GHz - SPW 1

os.system("rm -rf w51_C18O_21_merge7m12m_nocontsub.*")
clean(vis= 'w51_concat_7m12m.spw1.merge',
 imagename = 'w51_C18O_21_merge7m12m_nocontsub',
 field = "w51",
 spw = '',
 mode = 'velocity',
 multiscale = [0,3,9,27,81,243],
 nchan = 75,
 start = '20km/s',
 width = '1.335km/s',
 restfreq = '219.560GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 10000,
 threshold = '10mJy', #req threshold 5.85 mJy, 38 antennas, 66min tos, pwv auto,1arcsec res, 0.997 MHz BW gives 1.3mJy!
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)

myimagebase= 'w51_C18O_21_merge7m12m_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


fitspw='0:100~550;650~750;1100~1400;1700~1800;2750~2900;3300~3400;3520~3691'
uvcontsub(vis='w51_concat_7m12m.spw1.merge',
         field='',
         spw='', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=True)

linesub='w51_concat_7m12m.spw1.merge.contsub'

os.system("rm -rf w51_H2CO_321_220_merge7m12m_contsub.*")
clean(vis= linesub,
 imagename = "w51_H2CO_321_220_merge7m12m_contsub",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 70,
 start = '20km/s',
 width = '1.0km/s',
 restfreq = '218.76007GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.5mJy', #req threshold 5.85 mJy, 38 antennas, 66min tos, pwv auto,1arcsec res, 0.728 MHz BW gives 1.5mJy!
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)

myimagebase= 'w51_H2CO_321_220_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


linesub='w51_concat_7m12m.spw1.merge.contsub'
os.system("rm -rf w51_C18O_21_merge7m12m_contsub.*")
clean(vis= linesub,
 imagename = "w51_C18O_21_merge7m12m_contsub",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 75,
 start = '20km/s',
 width = '1.335km/s',
 restfreq = '219.560GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 10000,
 threshold = '50mJy',
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 multiscale = [0,3,9,27,81,243],
 minpb=0.4,
 pbcor=False,
 robust=2.0)

myimagebase= 'w51_C18O_21_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

linesub='w51_concat_7m12m.spw1.merge.contsub'
os.system("rm -rf w51_SO_65-54_merge7m12m_contsub_hires.*")
clean(vis= linesub,
 imagename = "w51_SO_65-54_merge7m12m_contsub_hires",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 105,
 start = '0km/s',
 width = '1.335km/s',
 restfreq = '219.94944GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 10000,
 threshold = '35mJy',
 imsize = [3840,3840],
 cell = '0.052arcsec',
 weighting = 'briggs',
 minpb=0.4,
 pbcor=False,
 robust=0.0)

myimagebase= 'w51_SO_65-54_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


##SPW2
#CO(2-1) spectral line

os.system("rm -rf w51_12CO_21_merge7m12m_nocontsub.*")
clean(vis= 'w51_concat_7m12m.spw2.merge',
 imagename = 'w51_12CO_21_merge7m12m_nocontsub',
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 95,
 start = '20km/s',
 width = '1.265km/s',
 restfreq = '230.538GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.3mJy', #req threshold 5.85 mJy, 38 antennas, 66min tos, pwv auto,1arcsec res, 1.33km/s MHz BW gives 1.3mJy!
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)

myimagebase= 'w51_12CO_21_merge7m12m_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


os.system("rm -rf w51_13CS_54_merge7m12m_nocontsub.*")
clean(vis= 'w51_concat_7m12m.spw2.merge',
 imagename = 'w51_13CS_54_merge7m12m_nocontsub',
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 55,
 start = '30km/s',
 width = '1.265km/s',
 restfreq = '231.221GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.3mJy', #req threshold 5.85 mJy, 38 antennas, 66min tos, pwv auto,1arcsec res, 1.33km/s MHz BW gives 1.3mJy!
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)

myimagebase= 'w51_13CS_54_merge7m12m_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


fitspw='0:0~100;260~410;500~650;725~800;1300~1450;1900~2000;2350~2420;2510~2570;3200~3230;2580~3620'
uvcontsub(vis='w51_concat_7m12m.spw2.merge',
         field='',
         spw='', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=True)

linesub='w51_concat_7m12m.spw2.merge.contsub'

os.system("rm -rf w51_12CO_21_merge7m12m_contsub.*")
clean(vis= linesub,
 imagename = "w51_12CO_21_merge7m12m_contsub",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 95,
 start = '20km/s',
 width = '1.265km/s',
 restfreq = '230.538GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.3mJy',
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)

myimagebase= 'w51_12CO_21_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

linesub='w51_concat_7m12m.spw2.merge.contsub'
os.system("rm -rf w51_12CO_21_merge7m12m_contsub_hires.*")
clean(vis= linesub,
 imagename = "w51_12CO_21_merge7m12m_contsub_hires",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 250,
 start = '-100km/s',
 width = '1.265km/s',
 restfreq = '230.538GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 10000,
 threshold = '35mJy', # rms ~7mJy
 imsize = [3840,3840],
 cell = '0.052',
 weighting = 'briggs',
 minpb=0.4,
 pbcor=False,
 robust = 0.0)

myimagebase= 'w51_12CO_21_merge7m12m_contsub_hires'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

os.system("rm -rf w51_13CS_54_merge7m12m_contsub.*")
clean(vis= linesub,
 imagename = "w51_13CS_54_merge7m12m_contsub",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 55,
 start = '30km/s',
 width = '1.265km/s',
 restfreq = '231.221GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.3mJy',
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)


myimagebase= 'w51_13CS_54_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

os.system("rm -rf w51_h30alpha_merge7m12m_contsub.*")
clean(vis= linesub,
 imagename = "w51_h30alpha_merge7m12m_contsub",
 field = "w51",
 spw = '',
 mode = 'velocity',
 nchan = 250,
 start = '-150km/s',
 width = '1.265km/s',
 restfreq = '231.90093GHz',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.3mJy',
 imsize = [1280,1280],
 cell = '0.15arcsec',
 weighting = 'natural',
 minpb=0.4,
 pbcor=False,
 robust=2.0)



myimagebase= 'w51_h30alpha_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image
