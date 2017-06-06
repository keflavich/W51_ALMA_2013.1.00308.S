# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
import re
import glob


## this file is to be run last after other data has been looked at

thesteps = []

step_title = {0: 'Image cont natural',
              1: 'Image cont briggs R=0.5',
              2: 'Image cont uniform',
              3: 'Image cont super-uniform',
              4: 'Image cont tapered',
             }
#              20: ' glob all files and make fits'}


try:
    print 'List of steps to be executed ...', mysteps
    thesteps = mysteps
except:
    print 'global variable mysteps not set.'
if (thesteps==[]):
    thesteps = range(0,len(step_title))
    print 'Executing all steps: ', thesteps





visname='w51e2_w51n_cax_new.ms'
# Fields are:  0 -> w51e2  1 -> w51n

# spw maps due to 3 EB now - nov 2016 - was june but lost files
# 0  1  2  3  4  5  6  7  8  9
# 10 11 12 13 14 15 16 17 18 19
# 20 21 22 23 24 25 26 27 28 29

souname1 = 'W51e2_cont_'
souname2 = 'W51n_cont_'

spwcont1='0:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,1:20~30;188~200;380~420,3:20~30;80~87;300~320,7:52~64;226~238;300~320,9:105~107;150~170;492~497;550~580;600~610;1257~1263,10:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,11:20~30;188~200;380~420,13:20~30;80~87;300~320,17:52~64;226~238;300~320,19:105~107;150~170;492~497;550~580;600~610;1257~1263,20:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,21:20~30;188~200;380~420,23:20~30;80~87;300~320,27:52~64;226~238;300~320,29:105~107;150~170;492~497;550~580;600~610;1257~1263'
spwcont2='0:815~820;880~890;900~910;1465~1500;1640~1650,1:18~24;395~415,2:100~105;239~242,3:461~468,6:191~195;403~409,7:305~315,9:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732,10:815~820;880~890;900~910;1465~1500;1640~1650,11:18~24;395~415,12:100~105;239~242,13:461~468,16:191~195;403~409,17:305~315,19:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732,20:815~820;880~890;900~910;1465~1500;1640~1650,21:18~24;395~415,22:100~105;239~242,23:461~468,26:191~195;403~409,27:305~315,29:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732'

# Image continuum of the target source 1 and source 2

mystep = 0  ## NATURAL
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print 'Step ', mystep, step_title[mystep]
# CHECK NOISE LEVELLLLLLLLLLAND RESIDUALS
    cell='0.005arcsec'
    imagesize=5120
    thre='0.15mJy'## measured noise ~0.13mJy
    weighting_scheme = 'natural'
    os.system('rm -rf '+souname1+weighting_scheme+'*')
    os.system('rm -rf '+souname2+weighting_scheme+'*')

    tclean(vis=visname,
           spw = spwcont1,
           imagename = souname1+weighting_scheme,
           field='0',
           cell=cell,
           imsize=imagesize,
           outframe='LSRK',
           niter=20000,
           interactive=False,
           threshold=thre,
           weighting=weighting_scheme,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
           uvrange='>300m',
           mask=['box[[2000pix,2200pix],[2954pix,3102pix]]',
                 'box[[2421pix,933pix],[3099pix,1590pix]]'])  ## 2 boxes as a list


    tclean(vis=visname,
           spw = spwcont2,
           imagename = souname2+weighting_scheme,
           field='1',
           cell=cell,
           imsize=imagesize,
           outframe='LSRK',
           niter=20000,
           interactive=False,
           threshold=thre,
           weighting=weighting_scheme,
           #robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
           uvrange='>300m',
           mask = 'box[[2191pix,2254pix],[3947pix,3154pix]]')


mystep = 1  ## BRIGGS R=0.5
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print 'Step ', mystep, step_title[mystep]

    cell='0.005arcsec'
    imagesize=5120
    thre='0.1mJy'## measured noise ~0.13mJy
    weighting_scheme = 'briggs'
    os.system('rm -rf '+souname1+weighting_scheme+'*')
    os.system('rm -rf '+souname2+weighting_scheme+'*')

    tclean(vis=visname,
      spw = spwcont1,
      imagename = souname1+weighting_scheme,
      field='0',
      cell=cell,
      imsize=imagesize,
      outframe='LSRK',
      niter=20000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
      robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
      uvrange='>300m',
      mask= ['box[[2000pix,2200pix],[2954pix,3102pix]]','box[[2421pix,933pix],[3099pix,1590pix]]'])  ## 2 boxes as a list


    tclean(vis=visname,
      spw = spwcont2,
      imagename = souname2+weighting_scheme,
      field='1',
      cell=cell,
      imsize=imagesize,
      outframe='LSRK',
      niter=20000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
      robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
      uvrange='>300m',
      mask = 'box[[2191pix,2254pix],[3947pix,3154pix]]')

mystep = 2  ## uniform
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print 'Step ', mystep, step_title[mystep]

    cell='0.005arcsec'
    imagesize=5120
    thre='0.25mJy'## measured noise ~0.3mJy
    weighting_scheme = 'uniform'
    os.system('rm -rf '+souname1+weighting_scheme+'*')
    os.system('rm -rf '+souname2+weighting_scheme+'*')

    tclean(vis=visname,
      spw = spwcont1,
      imagename = souname1+weighting_scheme,
      field='0',
      cell=cell,
      imsize=imagesize,
      outframe='LSRK',
      niter=20000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
  #    robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
      uvrange='>300m',
      mask='W51e2cax.cont_super.mask')  ## 2 boxes as a list


    tclean(vis=visname,
      spw = spwcont2,
      imagename = souname2+weighting_scheme,
      field='1',
      cell=cell,
      imsize=imagesize,
      outframe='LSRK',
      niter=20000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
  #    robust=0.5,
      uvrange='>300m',
      mask = 'W51ncax.cont_super100.mask')




mystep = 3  ##superuniform
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print 'Step ', mystep, step_title[mystep]

    weighting_scheme = 'superuniform'
    os.system('rm -rf '+souname1+weighting_scheme+'*')
    os.system('rm -rf '+souname2+weighting_scheme+'*')

    #IMPLEMENT PARAMETERS FROM BOTTOM LIST
    tclean(vis=visname,
      spw = spwcont1,
      imagename = souname1+weighting_scheme,
      field='0',
      cell=cell,
      imsize=imagesize,
      outframe='LSRK',
      niter=20000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
  #    robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
      uvrange='>300m',
      mask= 'W51e2cax.cont_super.mask')  ## 2 boxes as a list


    tclean(vis=visname,
      spw = spwcont2,
      imagename = souname2+weighting_scheme,
      field='1',
      cell=cell,
      imsize=imagesize,
      outframe='LSRK',
      niter=20000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
  #    robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
      uvrange='>300m',
      mask = 'W51ncax.cont_super100.mask')


#os.system('rm -rf W51ncax.cont_super100*')


#clean(vis="w51e2_w51n_cax_new.ms",imagename="W51ncax.cont_super100",outlierfile="",field="1",spw="0:815~820;880~890;900~910;1465~1500;1640~1650,1:18~24;395~415,2:100~105;239~242,3:461~468,6:191~195;403~409,7:305~315,9:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732,10:815~820;880~890;900~910;1465~1500;1640~1650,11:18~24;395~415,12:100~105;239~242,13:461~468,16:191~195;403~409,17:305~315,19:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732,20:815~820;880~890;900~910;1465~1500;1640~1650,21:18~24;395~415,22:100~105;239~242,23:461~468,26:191~195;403~409,27:305~315,29:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732",selectdata=True,timerange="",uvrange=">10m",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=1000,gain=0.1,threshold="0.5mJy",psfmode="clark",imagermode="csclean",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0],negcomponent=-1,smallscalebias=0.6,interactive=True,mask="",nchan=-1,start=0,width=1,outframe="LSRK",veltype="radio",imsize=5120,cell="0.002arcsec",phasecenter="",restfreq="",stokes="I",weighting="superuniform",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,flatnoise=True,allowchunk=False)

mystep = 4  ## tapered
if(mystep in thesteps):
    casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
    print 'Step ', mystep, step_title[mystep]
# CHECK NOISE LEVELLLLLLLLLLAND RESIDUALS
    cell='0.005arcsec'
    imagesize=5120
    thre='0.5mJy'
    weighting_scheme = 'natural'
    os.system('rm -rf '+souname1+weighting_scheme+'tapered*')
    os.system('rm -rf '+souname2+weighting_scheme+'tapered*')

    tclean(vis=visname,
           spw = spwcont1,
           imagename = souname1+weighting_scheme+"tapered",
           field='0',
           cell=cell,
           imsize=imagesize,
           outframe='LSRK',
           niter=20000,
           interactive=False,
           threshold=thre,
           weighting=weighting_scheme,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
           scales=[0,3,9,15],
           uvrange='<3000m',
           )


    tclean(vis=visname,
           spw = spwcont2,
           imagename = souname2+weighting_scheme+"tapered",
           field='1',
           cell=cell,
           imsize=imagesize,
           outframe='LSRK',
           niter=20000,
           interactive=False,
           threshold=thre,
           weighting=weighting_scheme,
           #robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
           scales=[0,3,9,15],
           uvrange='<3000m',
          )

