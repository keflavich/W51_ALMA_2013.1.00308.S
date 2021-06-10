invisname='w51e2_w51n_cax_new.ms'
visname='w51e2_w51n_cax_selfcal.ms'
if not os.path.exists(visname):
    import shutil
    shutil.copytree(invisname, visname)
    # no corrected col exsits
    #split(vis=invisname, outputvis=visname, datacolumn='corrected')
# Fields are:  0 -> w51e2  1 -> w51n

# spw maps due to 3 EB now - nov 2016 - was june but lost files
# 0  1  2  3  4  5  6  7  8  9
# 10 11 12 13 14 15 16 17 18 19
# 20 21 22 23 24 25 26 27 28 29

souname1 = 'W51e2_cont_big'
souname2 = 'W51n_cont_big'

spwcont1='0:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,1:20~30;188~200;380~420,3:20~30;80~87;300~320,7:52~64;226~238;300~320,9:105~107;150~170;492~497;550~580;600~610;1257~1263,10:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,11:20~30;188~200;380~420,13:20~30;80~87;300~320,17:52~64;226~238;300~320,19:105~107;150~170;492~497;550~580;600~610;1257~1263,20:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,21:20~30;188~200;380~420,23:20~30;80~87;300~320,27:52~64;226~238;300~320,29:105~107;150~170;492~497;550~580;600~610;1257~1263'
spwcont2='0:815~820;880~890;900~910;1465~1500;1640~1650,1:18~24;395~415,2:100~105;239~242,3:461~468,6:191~195;403~409,7:305~315,9:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732,10:815~820;880~890;900~910;1465~1500;1640~1650,11:18~24;395~415,12:100~105;239~242,13:461~468,16:191~195;403~409,17:305~315,19:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732,20:815~820;880~890;900~910;1465~1500;1640~1650,21:18~24;395~415,22:100~105;239~242,23:461~468,26:191~195;403~409,27:305~315,29:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732'

# Image continuum of the target source 1 and source 2

cell='0.005arcsec'
imagesize=5120
imagesize=12000
thre='0.3mJy'## measured noise ~0.13mJy

tclean(vis=visname,
       spw = spwcont1,
       imagename = 'w51e2_cont_big_robust0',
       field='w51e2',
       cell='0.005arcsec',
       imsize=12000,
       outframe='LSRK',
       niter=200000,
       interactive=False,
       threshold=thre,
       weighting='briggs',
       robust=0,
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       mask='w51e_b6_longbaseline_cleanregions.crtf',
       pblimit=0.01,
       savemodel='modelcolumn',
      )


tclean(vis=visname,
       spw = spwcont2,
       imagename = 'w51n_cont_big_robust0',
       field='w51n',
       cell='0.005arcsec',
       imsize=12000,
       outframe='LSRK',
       niter=200000,
       interactive=False,
       threshold=thre,
       weighting='briggs',
       robust=0,
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       mask='w51n_b6_longbaseline_cleanregions.crtf',
       pblimit=0.01,
       savemodel='modelcolumn',
       )

gaincal(vis=visname, caltable='w51_b6_selfcal_phase1_T.cal', calmode="p",
        gaintype="T", solint="inf", solnorm=True)

thesteps=(0,1,2,3)
mystep = 1
if(mystep in thesteps):

    cell='0.005arcsec'
    imagesize=5120
    imagesize=12000
    thre='0.5mJy'## measured noise ~0.13mjy
    weighting_scheme = 'natural'
    os.system('rm -rf '+souname1+weighting_scheme+'*')
    os.system('rm -rf '+souname2+weighting_scheme+'*')

    tclean(vis=visname,
      spw = spwcont1,
      imagename = souname1+weighting_scheme,
      field='0',
      cell=cell,
      imsize=imagesize,
      outframe='lsrk',
      niter=200000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       mask='w51e_b6_longbaseline_cleanregions.crtf',
       pblimit=0.01,
      #uvrange='>300m',
     # mask= ['box[[2000pix,2200pix],[2954pix,3102pix]]','box[[2421pix,933pix],[3099pix,1590pix]]'])  ## 2 boxes as a list
    )


    tclean(vis=visname,
      spw = spwcont2,
      imagename = souname2+weighting_scheme,
      field='1',
      cell=cell,
      imsize=imagesize,
      outframe='lsrk',
      niter=200000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       mask='w51n_b6_longbaseline_cleanregions.crtf',
       pblimit=0.01,
      #uvrange='>300m',
      #mask = 'box[[2191pix,2254pix],[3947pix,3154pix]]')
      )

mystep = 2  ## uniform
if(mystep in thesteps):

    cell='0.005arcsec'
    imagesize=5120
    imagesize=12000
    thre='0.5mJy'## measured noise ~0.3mjy
    weighting_scheme = 'uniform'
    os.system('rm -rf '+souname1+weighting_scheme+'*')
    os.system('rm -rf '+souname2+weighting_scheme+'*')

    tclean(vis=visname,
      spw = spwcont1,
      imagename = souname1+weighting_scheme,
      field='0',
      cell=cell,
      imsize=imagesize,
      outframe='lsrk',
      niter=200000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
       mask='w51e_b6_longbaseline_cleanregions.crtf',
  #    robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
       pblimit=0.01,
      #uvrange='>300m',
      #mask='w51e2cax.cont_super.mask'  ## 2 boxes as a list
          )


    tclean(vis=visname,
      spw = spwcont2,
      imagename = souname2+weighting_scheme,
      field='1',
      cell=cell,
      imsize=imagesize,
      outframe='lsrk',
      niter=200000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
       mask='w51n_b6_longbaseline_cleanregions.crtf',
       pblimit=0.01,
  #    robust=0.5,
      #uvrange='>300m',
      #mask = 'w51ncax.cont_super100.mask'
          )

mystep = 3  ##superuniform
if(mystep in thesteps):

    weighting_scheme = 'superuniform'
    os.system('rm -rf '+souname1+weighting_scheme+'*')
    os.system('rm -rf '+souname2+weighting_scheme+'*')

    #implement parameters from bottom list
    tclean(vis=visname,
      spw = spwcont1,
      imagename = souname1+weighting_scheme,
      field='0',
       mask='w51e_b6_longbaseline_cleanregions.crtf',
      cell=cell,
      imsize=imagesize,
      outframe='lsrk',
      niter=200000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
  #    robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
       pblimit=0.01,
      #uvrange='>300m',
      #mask= 'w51e2cax.cont_super.mask'
          )  ## 2 boxes as a list


    tclean(vis=visname,
      spw = spwcont2,
      imagename = souname2+weighting_scheme,
      field='1',
      cell=cell,
      imsize=imagesize,
      outframe='lsrk',
      niter=200000,
      interactive=False,
      threshold=thre,
      weighting=weighting_scheme,
  #    robust=0.5,
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
       mask='w51n_b6_longbaseline_cleanregions.crtf',
       pblimit=0.01,
      #mask = 'w51ncax.cont_super100.mask'
          )


