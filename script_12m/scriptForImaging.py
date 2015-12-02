#Imaging script

########################################
#Check CASA version

import re

if(re.search('^4.4', casadef.casa_version))  == None:
   sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.2 or 4.3')


########################################
#Definitions

msnames=['uid___A002_X9ee74a_X26f0.ms.split.cal','uid___A002_Xa8df68_Xb4b.ms.split.cal']

#EBX26f0  
#
#Fields:40
# ID   Code Name                RA               Decl           Epoch   SrcId      nRows
# 0    none J1751+0939          17:51:32.818570 +09.39.00.72840 J2000   0         148200
# 2    none Titan               16:07:42.286521 -18.43.07.40736 J2000   2          74100
# 3    none J1922+1530          19:22:34.699390 +15.30.10.03180 J2000   3         103740
# 4    none w51                 19:23:41.629000 +14.30.42.38000 J2000   4          29640
# 5    none w51                 19:23:38.489426 +14.30.54.59611 J2000   4          29640
# ... till 40 all w51 mosaic
#SpectralWindows:  (4 unique spectral windows and 1 unique polarization setups)
# SpwID  Name                           #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs  
# 0      ALMA_RB_06#BB_1#SW-01#FULL_RES   3840   TOPO  218604.028      -122.070    468750.0 218369.7142        1  XX  YY
# 1      ALMA_RB_06#BB_2#SW-01#FULL_RES   3840   TOPO  220296.833      -488.281   1875000.0 219359.5769        2  XX  YY
# 2      ALMA_RB_06#BB_3#SW-01#FULL_RES   3840   TOPO  230435.532       488.281   1875000.0 231372.7880        3  XX  YY
# 3      ALMA_RB_06#BB_4#SW-01#FULL_RES   3840   TOPO  233040.032       488.281   1875000.0 233977.2880        4  XX  YY 
#
#EBXb4b
#Fields:40
# ID   Code Name                RA               Decl           Epoch   SrcId      nRows
# 0    none J1733-1304          17:33:02.705790 -13.04.49.54820 J2000   0         180600
# 2    none Titan               15:47:13.128564 -18.00.17.55326 J2000   2          90300
# 3    none J1922+1530          19:22:34.699390 +15.30.10.03180 J2000   3          72240
# 4    none w51                 19:23:41.629000 +14.30.42.38000 J2000   4          18060
# 5    none w51                 19:23:38.489426 +14.30.54.59611 J2000   4          36120
#...till 40 all w51 mosaic
#SpectralWindows:  (4 unique spectral windows and 1 unique polarization setups)
# SpwID  Name                           #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs  
# 0      ALMA_RB_06#BB_1#SW-01#FULL_RES   3840   TOPO  218575.868      -122.070    468750.0 218341.5544        1  XX  YY
# 1      ALMA_RB_06#BB_2#SW-01#FULL_RES   3840   TOPO  220268.545      -488.281   1875000.0 219331.2895        2  XX  YY
# 2      ALMA_RB_06#BB_3#SW-01#FULL_RES   3840   TOPO  230406.218       488.281   1875000.0 231343.4741        3  XX  YY
# 3      ALMA_RB_06#BB_4#SW-01#FULL_RES   3840   TOPO  233010.718       488.281   1875000.0 233947.9741        4  XX  YY

#Concatenateddata
#2015-10-0112:19:05 INFO listobs	  0      ALMA_RB_06#BB_1#SW-01#FULL_RES   3840   TOPO  218604.028      -122.070    468750.0 218369.7142        1  XX  YY
#2015-10-0112:19:05 INFO listobs	  1      ALMA_RB_06#BB_2#SW-01#FULL_RES   3840   TOPO  220296.833      -488.281   1875000.0 219359.5769        2  XX  YY
#2015-10-0112:19:05 INFO listobs	  2      ALMA_RB_06#BB_3#SW-01#FULL_RES   3840   TOPO  230435.532       488.281   1875000.0 231372.7880        3  XX  YY
#2015-10-0112:19:05 INFO listobs	  3      ALMA_RB_06#BB_4#SW-01#FULL_RES   3840   TOPO  233040.032       488.281   1875000.0 233977.2880        4  XX  YY
#2015-10-0112:19:05 INFO listobs	  4      ALMA_RB_06#BB_1#SW-01#FULL_RES   3840   TOPO  218575.868      -122.070    468750.0 218341.5544        1  XX  YY
#2015-10-0112:19:05 INFO listobs	  5      ALMA_RB_06#BB_2#SW-01#FULL_RES   3840   TOPO  220268.545      -488.281   1875000.0 219331.2895        2  XX  YY
#2015-10-0112:19:05 INFO listobs	  6      ALMA_RB_06#BB_3#SW-01#FULL_RES   3840   TOPO  230406.218       488.281   1875000.0 231343.4741        3  XX  YY
#2015-10-0112:19:05 INFO listobs	  7      ALMA_RB_06#BB_4#SW-01#FULL_RES   3840   TOPO  233010.718       488.281   1875000.0 233947.9741        4  XX  YY


#Science Goal requirements: 
#RMS 5.85 mJy in 0.728 MHz bandwidth at rep. frequency 218.352 GHz
#Angular resolution 1.0 arcsec, later modified to 0.2 arcsec 

#Summaryof rms and beamsize of the produced images, based on concatenated data of the two Execution blocks. RMS in spectral line cubes is the average of four line free channels.  

#Target:w51
#Spw0, line  H2CO 3(0,3)-2(0,2), detected
#nocontinuum subtraction: 4.2 mJy/beam, beamsize 0.59x0.51 arcsec
#continuumsubtraction: 3.8 mJy/beam, beamsize 0.59x0.51 arcsec
#Spw0, line  H2CO 3(0,3)-2(0,2), high resolution dataset from EB Xb4b 
#continuumsubtraction: 4.8 mJy/beam, beamsize 0.31x0.27 arcsec

#Spw1, line C18O(2-1), detected
#nocontinuum subtraction: 3.7 mJy/beam, beamsize 0.68x0.53 arcsec
#continuumsubtraction: 3.1 mJy/beam, beamsize 0.68x0.53 arcsec

#Spw2, line 12CO(2-1), detected
#nocontinuum subtraction: 5.0 mJy/beam, beamsize 0.59x0.51 arcsec
#continuumsubtraction:  4.5 mJy/beam, beamsize 0.59x0.51 arcsec
#Spw2, 13CS(5-4), detected, 
#nocontinuum subtraction: 4.5 mJy/beam, beamsize 0.59x0.51 arcsec
#continuumsubtraction:  3.8 mJy/beam, beamsize 0.59x0.51 arcsec

#Spw3 was designated for continuum and PN(5-4) and neither could be imaged:
#1)there are locations in the mosaic where there are too many lines to obtain confidently the continuum emission, 
#2)the central frequency of spw 3 was changed due to technical reasons (JIRA 1708) causing the PN(5-4) line to be outside of the frequency range of spw 3.



##################################
#Concantenate
os.system('rm -rf w51_concat.ms.split.cal')
concat(vis=msnames,concatvis='w51_concat.ms.split.cal', freqtol='',copypointing=False)

#try continuum in a full spw, ignoring line contamination:

os.system("rm -rf w51_cont_spw3_hires.*")
clean(vis= 'w51_concat.ms.split.cal',
 imagename = "w51_cont_spw3_hires",
 field = "w51",
 spw = '3,7',
 mode = 'mfs',
 outframe = 'LSRK',
 interpolation = 'linear',
 imagermode='mosaic',
 interactive = False,
 niter = 10000,
 threshold = '20.0mJy', # got to 11 mJy, but it showed really bad instabilities.  15 was still too deep
 imsize = [2560,2560],
 cell = '0.052arcsec',
 weighting = 'briggs',
 pbcor=False,
 robust = 0.0)

#Imaging parameters:

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

fitspw='0:403~454;2893~3056;3821~2839'
uvcontsub(vis='uid___A002_Xa8df68_Xb4b.ms.split.cal',
         field='4~40',
         spw='0', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=False) 

os.system("rm -rf uid___A002_Xa8df68_Xb4b.ms.split.cal.spw0.contsub")
os.system("mv uid___A002_Xa8df68_Xb4b.ms.split.cal.contsub uid___A002_Xa8df68_Xb4b.ms.split.cal.spw0.contsub")
linesub='uid___A002_Xa8df68_Xb4b.ms.split.cal.spw0.contsub' #result of uvcontsub contains only the selected uvdata by the task

os.system("rm -rf w51_Xb4b_H2CO_303_202_contsub.*")
clean(vis= linesub,
 imagename = "w51_Xb4b_H2CO_303_202_contsub",
 field = "w51",
 spw = '0', 
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
 imsize = [2400,2400], #144 arcsec by 144 arcsec mosaic
 cell = '0.06arcsec',  #synth beam expected to be 0.2 arcsec, so 0.2/3= 0.06 arcsec cell
 weighting = 'natural',
 pbcor=False,
 robust = 0.0)




#H2CO(3(0,3)-2(0,2) - concatenated datasets

os.system("rm -rf w51_H2CO_303_202_nocontsub.*")
clean(vis= 'w51_concat.ms.split.cal',
 imagename = 'w51_H2CO_303_202_nocontsub',
 field = "w51",
 spw = '0,4',
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

fitspw='0:403~454;2893~3056;3821~2839,4:403~454;2893~3056;3821~2839'
uvcontsub(vis='w51_concat.ms.split.cal',
         field='4~40',
         spw='0,4', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=False) 

os.system("rm -rf w51_concat.ms.split.cal.spw0.contsub")
os.system("mv w51_concat.ms.split.cal.contsub w51_concat.ms.split.cal.spw0.contsub")
linesub='w51_concat.ms.split.cal.spw0.contsub' #result of uvcontsub contains only the selected uvdata by the task


os.system("rm -rf w51_H2CO_303_202_contsub.*")
clean(vis= linesub,
 imagename = "w51_H2CO_303_202_contsub",
 field = "w51",
 spw = '0,1', 
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

linesub='w51_concat.ms.split.cal.spw0.contsub' #result of uvcontsub contains only the selected uvdata by the task
os.system("rm -rf w51_H2CO_303_202_contsub_vresolved.*")
clean(vis= linesub,
 imagename = "w51_H2CO_303_202_contsub_vresolved",
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


#C18O(2-1)- Rest 219.560 GHz - SPW 1

os.system("rm -rf w51_C18O_21_nocontsub.*")
clean(vis= 'w51_concat.ms.split.cal',
 imagename = 'w51_C18O_21_nocontsub',
 field = "w51",
 spw = '1,5',
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

fitspw='1:174~203;875~982;1371~1402;2414~2590;2992~3127,5:174~203;875~982;1371~1402;2414~2590;2992~3127'
uvcontsub(vis='w51_concat.ms.split.cal',
         field='4~40',
         spw='1,5', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=False) 

os.system("rm -rf w51_concat.ms.split.cal.spw1.contsub")
os.system("mv w51_concat.ms.split.cal.contsub w51_concat.ms.split.cal.spw1.contsub")
linesub='w51_concat.ms.split.cal.spw1.contsub'

os.system("rm -rf w51_H2CO_321_220_contsub.*")
clean(vis= linesub,
 imagename = "w51_H2CO_321_220_contsub",
 field = "w51",
 spw = '0,1', 
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


os.system("rm -rf w51_C18O_21_contsub.*")
clean(vis= linesub,
 imagename = "w51_C18O_21_contsub",
 field = "w51",
 spw = '0,1', 
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

##SPW2
#CO(2-1) spectral line

os.system("rm -rf w51_12CO_21_nocontsub.*")
clean(vis= 'w51_concat.ms.split.cal',
 imagename = 'w51_12CO_21_nocontsub',
 field = "w51",
 spw = '2,6',
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

os.system("rm -rf w51_13CS_54_nocontsub.*")
clean(vis= 'w51_concat.ms.split.cal',
 imagename = 'w51_13CS_54_nocontsub',
 field = "w51",
 spw = '2,6',
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

fitspw='2:255~409;528~676;2313~2584;3356~3482,6:255~409;528~676;2313~2584;3356~3482'
uvcontsub(vis='w51_concat.ms.split.cal',
         field='4~40',
         spw='2,6', # spw to do continuum subtraction on
         fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
         solint='int',
         fitorder=1,
         want_cont=False)

os.system("rm -rf w51_concat.ms.split.cal.spw2.contsub")
os.system("mv w51_concat.ms.split.cal.contsub w51_concat.ms.split.cal.spw2.contsub")
linesub='w51_concat.ms.split.cal.spw2.contsub'

os.system("rm -rf w51_12CO_21_contsub.*")
clean(vis= linesub,
 imagename = "w51_12CO_21_contsub",
 field = "w51",
 spw = '0,1', 
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

os.system("rm -rf w51_12CO_21_contsub_hires.*")
clean(vis= linesub,
 imagename = "w51_12CO_21_contsub_hires",
 field = "w51",
 spw = '0,1', 
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

os.system("rm -rf w51_13CS_54_contsub.*")
clean(vis= linesub,
 imagename = "w51_13CS_54_contsub",
 field = "w51",
 spw = '0,1', 
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

os.system("rm -rf w51_h41alpha_contsub.*")
clean(vis= linesub,
 imagename = "w51_h41alpha_contsub",
 field = "w51",
 spw = '0,1', 
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


##SPW 3
#Continuumspectral window - there are too many lines to define the continuum.
#ThePN(5-4) line is not in the frequency range of the spectral window



#write out pbcorrected fits files and convert cleaned images to fits
#
#----------spw 0
myimagebase= 'w51_H2CO_303_202_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


myimagebase= 'w51_H2CO_303_202_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

myimagebase= 'w51_H2CO_303_202_contsub_vresolved'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

myimagebase= 'w51_H2CO_321_220_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

#----------spw1

myimagebase= 'w51_C18O_21_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


myimagebase= 'w51_C18O_21_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

#----------spw2

myimagebase= 'w51_12CO_21_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


myimagebase= 'w51_12CO_21_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

myimagebase= 'w51_12CO_21_contsub_hires'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


myimagebase= 'w51_13CS_54_nocontsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image


myimagebase= 'w51_13CS_54_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

myimagebase= 'w51_h41alpha_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image

#----------testhigh resolution

myimagebase= 'w51_Xb4b_H2CO_303_202_contsub'
impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor',fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
exportfits(imagename=myimagebase+'.flux',fitsimage=myimagebase+'.flux.fits') # export the PB image





