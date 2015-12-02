# these are made in scriptForImaging_fullcube
# w51_concat_7m12m.spw0.merge/  w51_concat_7m12m.spw1.merge/  w51_concat_7m12m.spw2.merge/  w51_concat_7m12m.spw3.merge/


#Wavelength ~ 0.137 mm
#primary beam ~ 28.3 arcsec, however mode=mosaic with roughly 140 by 140 arcsec dimensions
#max uvdistance is 1500 m (Xb4b EB), however the majority of the data has shorter uvdistances
#cell size 0.15 arcsec is enough to sample the beam of 0.5 arcsec, thus imagesize is [960,960]
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
         want_cont=False) 

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
 niter = 5000,
 threshold = '1.8mJy', #req rms 5.85 mJy, 38 antennas, 33.8min tos, pwv 2 (for EB Xb4b),0.2 arcsec res, 0.728 MHz BW gives 1.8mJy!   
 imsize = [2560,2560], #144 arcsec by 144 arcsec mosaic
 cell = '0.052arcsec',  #synth beam expected to be 0.2 arcsec, so 0.2/3= 0.06 arcsec cell
 weighting = 'uniform',
 pbcor=False,
 robust = 0.0)


#----------testhigh resolution

myimagebase= 'w51_merge7m12m_H2CO_303_202_contsub_uniform'
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
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)


myimagebase= 'w51_H2CO_303_202_merge7m12m_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

os.system("rm -rf w51_H2CO_303_202_merge7m12m_contsub.*")
clean(vis= linesub,
 imagename = "w51_H2CO_303_202_merge7m12m_contsub",
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
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)

myimagebase= 'w51_H2CO_303_202_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

os.system("rm -rf w51_H2CO_303_202_merge7m12m_contsub_vresolved.*")
clean(vis= linesub,
 imagename = "w51_H2CO_303_202_merge7m12m_contsub_vresolved",
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
 niter = 5000,
 threshold = '3.4mJy',
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)


myimagebase= 'w51_H2CO_303_202_merge7m12m_contsub_vresolved'
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
 nchan = 75,
 start = '20km/s', 
 width = '1.335km/s',
 restfreq = '219.560GHz',
 outframe = 'LSRK',
 interpolation = 'linear', 
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.3mJy', #req threshold 5.85 mJy, 38 antennas, 66min tos, pwv auto,1arcsec res, 0.997 MHz BW gives 1.3mJy!   
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)

myimagebase= 'w51_C18O_21_merge7m12m_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


fitspw='0:174~203;875~982;1371~1402;2414~2590;2992~3127'
uvcontsub(vis='w51_concat_7m12m.spw1.merge',
         field='',
         spw='', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=False) 

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
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)

myimagebase= 'w51_H2CO_321_220_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


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
 niter = 5000,
 threshold = '1.3mJy',    
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)

myimagebase= 'w51_C18O_21_merge7m12m_contsub'
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
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)

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
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)

myimagebase= 'w51_13CS_54_merge7m12m_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


fitspw='0:255~409;528~676;2313~2584;3356~3482'
uvcontsub(vis='w51_concat_7m12m.spw2.merge',
         field='',
         spw='', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=False)

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
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)

myimagebase= 'w51_12CO_21_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

os.system("rm -rf w51_12CO_21_merge7m12m_contsub_hires.*")
clean(vis= linesub,
 imagename = "w51_12CO_21_merge7m12m_contsub_hires",
 field = "w51",
 spw = '', 
 mode = 'velocity',
 nchan = 250,
 start = '-150km/s', 
 width = '1.265km/s',
 restfreq = '230.538GHz',
 outframe = 'LSRK',
 interpolation = 'linear', 
 imagermode='mosaic',
 interactive = False,
 niter = 5000,
 threshold = '1.3mJy',    
 imsize = [2560,2560],
 cell = '0.052',
 weighting = 'briggs',
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
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)


myimagebase= 'w51_13CS_54_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

os.system("rm -rf w51_h41alpha_merge7m12m_contsub.*")
clean(vis= linesub,
 imagename = "w51_h41alpha_merge7m12m_contsub",
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
 imsize = [960,960],
 cell = '0.15arcsec',
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)



myimagebase= 'w51_h41alpha_merge7m12m_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image
