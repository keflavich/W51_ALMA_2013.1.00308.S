

########################################
# Check CASA version

import re
import casadef

if casadef.casa_version < '4.4.0' :
    sys.exit("Please use CASA version greater than or equal to 4.4.0 with this script")


##################################################
# Create an Averaged Continuum MS


finalvis='calibrated_final.ms' # This is your output ms from the data
                               # preparation script.
contvis='calibrated_final_cont.ms'

contspws_lo = '0,1,4,5,8,9,12,13,16,17'
contspws_hi = '2,3,6,7,10,11,14,15,18,19'

if not os.path.exists(contvis):
    # Use plotms to identify line and continuum spectral windows.
    plotms(vis=finalvis, xaxis='channel', yaxis='amplitude',
           ydatacolumn='data',
           avgtime='1e8', avgscan=True, avgchannel='1', 
           iteraxis='spw' )






    # Set spws to be used to form continuum
    contspws = '0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19'

    # If you have complex line emission and no dedicated continuum
    # windows, you will need to flag the line channels prior to averaging.
    flagmanager(vis=finalvis,mode='save',
                versionname='before_cont_flags')

    initweights(vis=finalvis,wtmode='weight',dowtsp=True)

    # Flag the "line channels"
    flagchannels='0:62~76;180~192;219~247;317~334;366~390;465~480;512~519;548~555;572~580;605~617;649~669;795~806;812~826;925~931;948~960;1125~1132;1209~1239;1250~1259;1345~1380;1398~1420;1453~1467;1479~1497;1558~1573;1596~1612;1716~1739;1796~1816;1847~1885;2063~2078;2442~2606;2793~2813;2928~2936;3105~3122;3203~3223;3261~3285;3663~3698;3770~3787,1:85~96;146~180;196~211;380~413;560~577;887~897;911~922;946~954;961~970;1213~1232;1414~1435;1490~1521;1665~1675;1686~1693;1740~1749,2:268~283;369~382;393~401;675~695;728~740;880~906;971~999;1135~1159;1233~1254;1311~1343;1384~1443;1455~1479;1504~1524;1695~1708;1722~1807,3:0~8;138~183;229~280;312~347;430~444;466~485;520~536;544~560;672~692;719~738;819~831;984~995;1029~1050;1184~1213;1251~1292;1426~1472;1581~1592;1608~1621;1651~1666;1678~1692,4:62~77;135~153;176~192;238~246;366~386;607~618;650~668;793~811;1211~1237;1355~1385;1401~1417;1454~1467;1481~1500;1563~1574;1597~1612;1718~1736;1795~1815;1848~1885;2040~2083;2459~2605;2795~2808;3204~3232;3257~3287;3582~3595;3665~3697;3765~3787,5:81~93;147~174;197~206;379~411;561~576;848~863;882~893;909~921;943~975;1216~1233;1419~1432;1496~1519;1665~1697,6:269~292;347~417;675~695;728~741;874~905;976~1004;1130~1156;1237~1247;1315~1342;1376~1413;1424~1440;1462~1475;1508~1521;1697~1711;1720~1804,7:0~10;149~181;234~251;262~274;429~439;519~530;676~694;716~740;822~835;982~992;1029~1051;1180~1208;1254~1292;1428~1472;1611~1620;1651~1667;1677~1692,8:61~84;177~195;218~250;366~392;543~561;601~616;644~665;794~807;946~963;1209~1239;1351~1381;1399~1415;1450~1465;1483~1496;1564~1578;1593~1611;1714~1734;1800~1812;1843~1881;2061~2078;2441~2605;2789~2814;3187~3218;3249~3292;3657~3698;3754~3787,9:85~95;147~179;196~213;377~416;560~574;880~979;1220~1235;1413~1432;1492~1521;1669~1696,10:272~283;367~381;396~417;667~703;727~739;875~912;971~1000;1130~1161;1227~1254;1274~1291;1314~1337;1379~1480;1511~1522;1724~1808,11:0~12;137~182;231~282;314~350;520~534;666~739;982~998;1030~1050;1253~1292;1425~1477;1605~1629;1653~1665;1681~1693,12:55~81;172~197;236~246;360~392;608~618;650~664;1213~1246;1357~1415;1449~1497;1719~1735;1801~1816;1848~1888;2057~2083;2470~2600;2790~2816;3204~3294;3658~3704;3763~3792,13:79~94;144~180;197~213;374~413;562~574;885~953;1220~1230;1414~1432;1492~1523;1663~1697,14:269~292;371~409;669~695;875~891;982~996;1127~1162;1305~1350;1371~1484;1720~1807,15:0~9;139~192;230~277;519~528;670~693;718~730;1008~1052;1256~1294;1420~1468;1603~1625;1649~1666,16:62~80;176~195;366~392;607~621;644~663;799~807;1214~1235;1357~1383;1400~1415;1498~1502;1593~1619;1714~1736;1793~1816;1846~1891;2064~2083;2451~2605;2789~2816;3202~3222;3261~3286;3666~3697;3760~3800,17:79~94;155~178;197~208;348~421;562~577;879~973;1217~1233;1416~1433;1490~1521;1668~1696,18:268~283;367~404;674~695;878~908;978~999;1131~1160;1236~1247;1307~1336;1372~1443;1461~1479;1503~1521;1697~1809,19:0~13;130~187;231~283;672~689;722~735;980~999;1031~1058;1254~1297;1425~1471;1605~1624;1647~1666;1681~1693'


    flagdata(vis=finalvis,mode='manual',
              spw=flagchannels,flagbackup=False)

    # check that flags are as expected, NOTE must check reload on plotms
    # gui if its still open.
    plotms(vis=finalvis,yaxis='amp',xaxis='channel',
           avgchannel='1',avgtime='1e8',avgscan=True,iteraxis='spw') 

    # Average the channels within spws
    rmtables(contvis)
    os.system('rm -rf ' + contvis + '.flagversions')


    split(vis=finalvis,
         spw=contspws,      
         outputvis=contvis,
          width=[256,128,128,128,256,128,128,128,256,128,128,128,256,128,128,128,256,128,128,128], # number of channels to average together. The final channel width should be less than 125MHz in Bands 3, 4, and 6 and 250MHz in Band 7.
         datacolumn='data')


    # Check the weights. You will need to change antenna and field to
    # appropriate values
    plotms(vis=contvis, yaxis='wtsp',xaxis='freq',spw='',antenna='DA61',field='w51n')

    # If you flagged any line channels, restore the previous flags
    flagmanager(vis=finalvis,mode='restore',
                versionname='before_cont_flags')

    # Inspect continuum for any problems
    plotms(vis=contvis,xaxis='uvdist',yaxis='amp',coloraxis='spw')

# #############################################
# Image Parameters


# source parameters
# ------------------

field='w51n' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL POTENTIALLY TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
gridder='standard' # uncomment if single field 

# image parameters.
# ----------------




cell='0.005arcsec' # cell size for imaging.
imsize = [14700,14700] # size of image in pixels.

# velocity parameters
# -------------------

outframe='lsrk' # velocity reference frame. See science goals.
veltype='radio' # velocity type. 


# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean. 

weighting = 'briggs'
robust=0
niter=1000
threshold = '0.0mJy'

#############################################
# Imaging the Continuuum

# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'         

contimagename = 'w51north_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual' 

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=contvis)
#delmod(vis=contvis)

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
        rmtables(contimagename+ext)



    tclean(vis=contvis,
           imagename=contimagename,
           field=field,
           #  phasecenter=phasecenter, # uncomment if mosaic.      
           specmode='mfs',
           #deconvolver='hogbom', 
           # Uncomment the below to image with nterms>1 when the fractional bandwidth is greater than 10%.
           deconvolver='mtmfs',
           nterms=2,
           imsize = imsize, 
           cell= cell, 
           weighting = weighting,
           robust = robust,
           niter = niter, 
           threshold = threshold,
           interactive = False,
           gridder = gridder,
           pbcor = True,
           mask='cont.mask')
           

    # 2018-07-31 01:58:18     INFO    SIImageStore::restore   Common Beam for chan : 0 : 0.0666256 arcsec, 0.0414155 arcsec, -44.5511 deg
    # 2018-07-31 01:57:53     INFO    task_tclean::SIImageStoreMultiTerm::printImageStats     [w51north_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual] Peak residual (max,min) within mask : (0.000433739,-0.00042061) over full image : (0.00133055,-0.000830991)
    # 2018-07-31 01:57:57     INFO    task_tclean::SIImageStoreMultiTerm::printImageStats     [w51north_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual] Total Model Flux : 0.156686(tt0)0.0910296(tt1)
    # 2018-07-31 01:57:57     INFO    tclean::::      Reached global stopping criterion : iteration limit


    # If you'd like to redo your clean, but don't want to make a new mask
    # use the following commands to save your original mask. This is an optional step.
    #contmaskname = 'cont.mask'
    ##rmtables(contmaskname) # if you want to delete the old mask
    #os.system('cp -ir ' + contimagename + '.mask ' + contmaskname)


contimagename = 'w51north_sci.spw{0}.mfs.I.manual'.format(contspws_hi.replace(",","_"))

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    tclean(vis=contvis,
           imagename=contimagename,
           field=field,
           spw=contspws_hi,
           specmode='mfs',
           #deconvolver='hogbom', 
           # Uncomment the below to image with nterms>1 when the fractional bandwidth is greater than 10%.
           deconvolver='mtmfs',
           nterms=2,
           imsize = imsize, 
           cell= cell, 
           weighting = weighting,
           robust = robust,
           niter = niter, 
           threshold = threshold,
           interactive = True,
           gridder = gridder,
           pbcor = True,
           mask='cont.mask')

contimagename = 'w51north_sci.spw{0}.mfs.I.manual'.format(contspws_lo.replace(",","_"))

if not os.path.exists(contimagename+'.image.tt0.pbcor'):
    tclean(vis=contvis,
           imagename=contimagename,
           field=field,
           spw=contspws_lo,
           specmode='mfs',
           #deconvolver='hogbom', 
           # Uncomment the below to image with nterms>1 when the fractional bandwidth is greater than 10%.
           deconvolver='mtmfs',
           nterms=2,
           imsize = imsize, 
           cell= cell, 
           weighting = weighting,
           robust = robust,
           niter = niter, 
           threshold = threshold,
           interactive = True,
           gridder = gridder,
           pbcor = True,
           mask='cont.mask')


    ########################################
    # Continuum Subtraction for Line Imaging




    #fitspw = '0:0~61;77~179;193~218;248~316;335~365;391~464;481~511;520~547;556~571;581~604;618~648;670~794;807~811;827~924;932~947;961~1124;1133~1208;1240~1249;1260~1344;1381~1397;1421~1452;1468~1478;1498~1557;1574~1595;1613~1715;1740~1795;1817~1846;1886~2062;2079~2441;2607~2792;2814~2927;2937~3104;3123~3202;3224~3260;3286~3662;3699~3769;3788~3839,1:0~84;97~145;181~195;212~379;414~559;578~886;898~910;923~945;955~960;971~1212;1233~1413;1436~1489;1522~1664;1676~1685;1694~1739;1750~1919,2:0~267;284~368;383~392;402~674;696~727;741~879;907~970;1000~1134;1160~1232;1255~1310;1344~1383;1444~1454;1480~1503;1525~1694;1709~1721;1808~1919,3:9~137;184~228;281~311;348~429;445~465;486~519;537~543;561~671;693~718;739~818;832~983;996~1028;1051~1183;1214~1250;1293~1425;1473~1580;1593~1607;1622~1650;1667~1677;1693~1919,4:0~61;78~134;154~175;193~237;247~365;387~606;619~649;669~792;812~1210;1238~1354;1386~1400;1418~1453;1468~1480;1501~1562;1575~1596;1613~1717;1737~1794;1816~1847;1886~2039;2084~2458;2606~2794;2809~3203;3233~3256;3288~3581;3596~3664;3698~3764;3788~3839,5:0~80;94~146;175~196;207~378;412~560;577~847;864~881;894~908;922~942;976~1215;1234~1418;1433~1495;1520~1664;1698~1919,6:0~268;293~346;418~674;696~727;742~873;906~975;1005~1129;1157~1236;1248~1314;1343~1375;1414~1423;1441~1461;1476~1507;1522~1696;1712~1719;1805~1919,7:11~148;182~233;252~261;275~428;440~518;531~675;695~715;741~821;836~981;993~1028;1052~1179;1209~1253;1293~1427;1473~1610;1621~1650;1668~1676;1693~1919,8:0~60;85~176;196~217;251~365;393~542;562~600;617~643;666~793;808~945;964~1208;1240~1350;1382~1398;1416~1449;1466~1482;1497~1563;1579~1592;1612~1713;1735~1799;1813~1842;1882~2060;2079~2440;2606~2788;2815~3186;3219~3248;3293~3656;3699~3753;3788~3839,9:0~84;96~146;180~195;214~376;417~559;575~879;980~1219;1236~1412;1433~1491;1522~1668;1697~1919,10:0~271;284~366;382~395;418~666;704~726;740~874;913~970;1001~1129;1162~1226;1255~1273;1292~1313;1338~1378;1481~1510;1523~1723;1809~1919,11:13~136;183~230;283~313;351~519;535~665;740~981;999~1029;1051~1252;1293~1424;1478~1604;1630~1652;1666~1680;1694~1919,12:0~54;82~171;198~235;247~359;393~607;619~649;665~1212;1247~1356;1416~1448;1498~1718;1736~1800;1817~1847;1889~2056;2084~2469;2601~2789;2817~3203;3295~3657;3705~3762;3793~3839,13:0~78;95~143;181~196;214~373;414~561;575~884;954~1219;1231~1413;1433~1491;1524~1662;1698~1919,14:0~268;293~370;410~668;696~874;892~981;997~1126;1163~1304;1351~1370;1485~1719;1808~1919,15:10~138;193~229;278~518;529~669;694~717;731~1007;1053~1255;1295~1419;1469~1602;1626~1648;1667~1919,16:0~61;81~175;196~365;393~606;622~643;664~798;808~1213;1236~1356;1384~1399;1416~1497;1503~1592;1620~1713;1737~1792;1817~1845;1892~2063;2084~2450;2606~2788;2817~3201;3223~3260;3287~3665;3698~3759;3801~3839,17:0~78;95~154;179~196;209~347;422~561;578~878;974~1216;1234~1415;1434~1489;1522~1667;1697~1919,18:0~267;284~366;405~673;696~877;909~977;1000~1130;1161~1235;1248~1306;1337~1371;1444~1460;1480~1502;1522~1696;1810~1919,19:14~129;188~230;284~671;690~721;736~979;1000~1030;1059~1253;1298~1424;1472~1604;1625~1646;1667~1680;1694~1919' # *line-free* channels for fitting continuum
    #linespw = '0~19' # line spectral windows. You can subtract the continuum from multiple spectral line windows at once.

    #finalvis='calibrated_final.ms'

    #uvcontsub(vis=finalvis,
    #          spw=linespw, # spw to do continuum subtraction on
    #          fitspw=fitspw, # regions without lines.
    #          excludechans=False, # fit the regions in fitspw
    #          combine='spw', 
    #          solint='int',
    #          fitorder=1,
    #          want_cont=False) # This value should not be changed.


    # NOTE: Imaging the continuum produced by uvcontsub with
    # want_cont=True will lead to extremely poor continuum images because
    # of bandwidth smearing effects. For imaging the continuum, you should
    # always create a line-free continuum data set using the process
    # outlined above.

    ##############################################
    # Image line emission [REPEAT AS NECESSARY]


finalvis = 'calibrated_final.ms'
linevis = finalvis # uncomment if you neither continuum subtracted nor self-calibrated your data.
#linevis = finalvis + '.contsub' # uncomment if continuum subtracted
# linevis = finalvis + '.contsub.selfcal' # uncommment if both continuum subtracted and self-calibrated
# linevis = finalvis + '.selfcal' # uncomment if just self-calibrated (no continuum subtraction)

robust=-1.5
#restfreq='85.98GHz' # Typically the rest frequency of the line of
#                        # interest. If the source has a significant
#                        # redshift (z>0.2), use the observed sky
#                        # frequency (nu_rest/(1+z)) instead of the
#                        # rest frequency of the
#                        # line.

# Skip GIGANTIC cube making...
if False:
    for spw in ('0,4,8,12,16',
                '1,5,9,13,17',
                '2,6,10,14,18',
                '3,7,11,15,19',
               ):

        lineimagename =  'w51north_sci.spw{0}.robustm1.5.cube.I.manual'.format(spw.replace(",","_"))

        #start='-250km/s' # start velocity. See science goals for appropriate value.
        #width='2km/s' # velocity width. See science goals.
        #nchan = 175  # number of channels. See science goals for appropriate value.


        # If necessary, run the following commands to get rid of older clean
        # data.

        #clearcal(vis=linevis)
        #delmod(vis=linevis)

        if not os.path.exists(lineimagename+".image"):
            for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
                rmtables(lineimagename + ext)

            tclean(vis=linevis,
                   imagename=lineimagename, 
                   field=field,
                   spw=spw,
                   # phasecenter=phasecenter, # uncomment if mosaic.      
                   specmode='cube',
                   outframe=outframe,
                   veltype=veltype, 
                   niter=niter,  
                   threshold=threshold, 
                   interactive=False,
                   cell=cell,
                   imsize=imsize, 
                   weighting=weighting,
                   robust=robust,
                   gridder=gridder,
                   pbcor=True,
                   restoringbeam='common',
                   chanchunks=-1) # break up large cubes automatically so that you don't run out of memory.

        # If you'd like to redo your clean, but don't want to make a new mask
        # use the following commands to save your original mask. This is an
        # optional step.
        # linemaskname = 'line.mask'
        ## rmtables(linemaskname) # uncomment if you want to overwrite the mask.
        # os.system('cp -ir ' + lineimagename + '.mask ' + linemaskname)

    ##############################################
    # Export the images

import glob

myimages = glob.glob("w51north*.pbcor")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("w51north*.pb")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True) 

##############################################
# Create Diagnostic PNGs

os.system("rm -rf *.png")
mycontimages = glob.glob("*mfs*manual.image")
for cimage in mycontimages:
    mymax=imstat(cimage)['max'][0]
    mymin=-0.1*mymax
    outimage = cimage+'.png'
    os.system('rm -rf '+outimage)
    imview(raster={'file':cimage,'range':[mymin,mymax]},out=outimage)

mylineimages = glob.glob("*cube*manual.image")
for limage in mylineimages:
    mom8=limage+'.mom8'
    os.system("rm -rf "+mom8)
    immoments(limage,moments=[8],outfile=mom8)
    mymax=imstat(mom8)['max'][0]
    mymin=-0.1*mymax
    os.system("rm -rf "+mom8+".png")
    imview(raster={'file':mom8,'range':[mymin,mymax]},out=mom8+'.png')


##############################################
# Analysis

# For examples of how to get started analyzing your data, see
#     https://casaguides.nrao.edu/index.php/TWHydraBand7_Imaging_4.3
#     
