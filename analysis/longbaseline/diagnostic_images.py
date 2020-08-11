def makefits(myimagebase):
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.alpha', fitsimage=myimagebase+'.alpha.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.alpha.error', fitsimage=myimagebase+'.alpha.error.fits', dropdeg=True, overwrite=True)




cell='0.01arcsec'
imagesize=2560
thre='50.0mJy'
weighting_scheme = 'natural'
visname='w51e2_w51n_cax_new.ms'
souname1 = 'W51e2_cont_'
souname2 = 'W51n_cont_'

spwcont1='0:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,1:20~30;188~200;380~420,3:20~30;80~87;300~320,7:52~64;226~238;300~320,9:105~107;150~170;492~497;550~580;600~610;1257~1263,10:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,11:20~30;188~200;380~420,13:20~30;80~87;300~320,17:52~64;226~238;300~320,19:105~107;150~170;492~497;550~580;600~610;1257~1263,20:88~92;167~170;285~295;525~530;814~817;872~875;943~947;1085~1090;1470~1478;1845~1850,21:20~30;188~200;380~420,23:20~30;80~87;300~320,27:52~64;226~238;300~320,29:105~107;150~170;492~497;550~580;600~610;1257~1263'
spwcont2='0:815~820;880~890;900~910;1465~1500;1640~1650,1:18~24;395~415,2:100~105;239~242,3:461~468,6:191~195;403~409,7:305~315,9:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732,10:815~820;880~890;900~910;1465~1500;1640~1650,11:18~24;395~415,12:100~105;239~242,13:461~468,16:191~195;403~409,17:305~315,19:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732,20:815~820;880~890;900~910;1465~1500;1640~1650,21:18~24;395~415,22:100~105;239~242,23:461~468,26:191~195;403~409,27:305~315,29:315~320;395~400;700~708;790~796;1380~1385;1560~1565;1725~1732'


for timerange,eb in [('2015/10/27/23:06:40~2015/10/27/25:20:00', '1'),
                     ('2015/10/30/20:00:00~2015/10/30/22:46:40', '2'),
                     ('2015/10/30/22:46:40~2015/10/31/22:46:40', '3'),
                    ]:

    imagename = souname1+weighting_scheme+'tapered300m_EB{0}'.format(eb)
    os.system('rm -rf '+imagename)

    tclean(vis=visname,
           spw = spwcont1,
           imagename = imagename,
           field='0',
           cell=cell,
           timerange=timerange,
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
           savemodel='none',
           )
    makefits(imagename)
