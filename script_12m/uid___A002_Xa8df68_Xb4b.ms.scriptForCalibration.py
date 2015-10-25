# ALMA Data Reduction Script

# Calibration

thesteps = []
step_title = {0: 'Import of the ASDM',
              1: 'Fix of SYSCAL table times',
              2: 'listobs',
              3: 'A priori flagging',
              4: 'Generation and time averaging of the WVR cal table',
              5: 'Generation of the Tsys cal table',
              6: 'Generation of the antenna position cal table',
              7: 'Application of the WVR, Tsys and antpos cal tables',
              8: 'Split out science SPWs and time average',
              9: 'Listobs, clear pointing table, and save original flags',
              10: 'Initial flagging',
              11: 'Putting a model for the flux calibrator(s)',
              12: 'Save flags before bandpass cal',
              13: 'Bandpass calibration',
              14: 'Save flags before gain cal',
              15: 'Gain calibration',
              16: 'Save flags before applycal',
              17: 'Application of the bandpass and gain cal tables',
              18: 'Split out corrected column',
              19: 'Save flags after applycal'}

if 'applyonly' not in globals(): applyonly = False
try:
  print 'List of steps to be executed ...', mysteps
  thesteps = mysteps
except:
  print 'global variable mysteps not set.'
if (thesteps==[]):
  thesteps = range(0,len(step_title))
  print 'Executing all steps: ', thesteps

# The Python variable 'mysteps' will control which steps
# are executed when you start the script using
#   execfile('scriptForCalibration.py')
# e.g. setting
#   mysteps = [2,3,4]# before starting the script will make the script execute
# only steps 2, 3, and 4
# Setting mysteps = [] will make it execute all steps.

import re

import os

if applyonly != True: es = aU.stuffForScienceDataReduction() 


if re.search('^4.4.0', casadef.casa_version) == None:
 sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.4.0')


# CALIBRATE_AMPLI: Titan
# CALIBRATE_ATMOSPHERE: J1733-1304,Titan,w51
# CALIBRATE_BANDPASS: J1733-1304
# CALIBRATE_FLUX: Titan
# CALIBRATE_FOCUS: 
# CALIBRATE_PHASE: J1922+1530
# CALIBRATE_POINTING: J1517-2422,J1733-1304,J1922+1530
# OBSERVE_TARGET: w51

# Using reference antenna = DA55

# Import of the ASDM
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  if os.path.exists('uid___A002_Xa8df68_Xb4b.ms') == False:
    importasdm('uid___A002_Xa8df68_Xb4b', asis='Antenna Station Receiver Source CalAtmosphere CalWVR CorrelatorMode SBSummary', bdfflags=True, lazy=False, process_caldevice=False)
  if applyonly != True: es.fixForCSV2555('uid___A002_Xa8df68_Xb4b.ms')

# Fix of SYSCAL table times
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  from recipes.almahelpers import fixsyscaltimes
  fixsyscaltimes(vis = 'uid___A002_Xa8df68_Xb4b.ms')

print "# A priori calibration"

# listobs
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.listobs')
  listobs(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    listfile = 'uid___A002_Xa8df68_Xb4b.ms.listobs')
  
  

# A priori flagging
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    mode = 'manual',
    spw = '1~24',
    autocorr = T,
    flagbackup = F)
  
  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    mode = 'manual',
    intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*',
    flagbackup = F)
  
  flagcmd(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'plot',
    plotfile = 'uid___A002_Xa8df68_Xb4b.ms.flagcmd.png')
  
  flagcmd(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'apply')

 #Antenna has a low gain
  flagdata(vis= 'uid___A002_Xa8df68_Xb4b.ms', mode= 'manual', antenna='DV12',flagbackup= F)


  flagmanager(vis = 'uid___A002_Xa8df68_Xb4b.ms',mode = 'save', versionname = 'Apriori')

 

# Generation and time averaging of the WVR cal table
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.wvr') 
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.wvrgcal') 
  
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_Xa8df68_Xb4b.ms.wvrgcal')
  
  wvrgcal(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.wvr',
    spw = [17, 19, 21, 23],
    smooth = '6.048s',
    toffset = 0,
    tie = ['w51,J1922+1530'],
    statsource = 'w51')
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: aU.plotWVRSolutions(caltable='uid___A002_Xa8df68_Xb4b.ms.wvr', spw='17', antenna='DA55',
    yrange=[-199,199],subplot=22, interactive=False,
    figfile='uid___A002_Xa8df68_Xb4b.ms.wvr.plots/uid___A002_Xa8df68_Xb4b.ms.wvr') 
  
  #Note: If you see wraps in these plots, try changing yrange or unwrap=True 
  #Note: If all plots look strange, it may be a bad WVR on the reference antenna.
  #      To check, you can set antenna='' to show all baselines.
  

# Generation of the Tsys cal table
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.tsys') 
  gencal(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    caltype = 'tsys')
  
  # Flagging edge channels
  
  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    mode = 'manual',
    spw = '9:0~3;124~127,11:0~3;124~127,13:0~3;124~127,15:0~3;124~127',
    flagbackup = F)
  
  if applyonly != True: aU.plotbandpass(caltable='uid___A002_Xa8df68_Xb4b.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='92.1875%',showfdm=True, showBasebandNumber=True, showimage=False, 
    field='', figfile='uid___A002_Xa8df68_Xb4b.ms.tsys.plots.overlayTimeBeforeTsysFlag/uid___A002_Xa8df68_Xb4b.ms.tsys') 
  
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.tsys', msName='uid___A002_Xa8df68_Xb4b.ms', interactive=False) 
  

# Flagging ms with bad tsys calibration 


  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    mode = 'manual',
    antenna = 'DA41',
    spw = '15:94~120',
    flagbackup = F)

  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    mode = 'manual',
    timerange='01:34:0~01:36:0',
    antenna = '',
    spw = '13:4~12',
    flagbackup = F)

  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    mode = 'manual',
    timerange='01:51:0~01:53:0',
    antenna = '',
    spw = '13:4~12',
    flagbackup = F)

  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    mode = 'manual',
    antenna = 'DV14',
    spw = '11',
    flagbackup = F)

  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    mode = 'manual',
    antenna = 'DV15',
    spw = '11:40~60',
    flagbackup = F)

  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    mode = 'manual',
    antenna = 'DV18',
    spw = '11:30~42',
    flagbackup = F)

  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.tsys',
    mode = 'manual',
    antenna = 'DV18',
    spw = '13:42~53',
    flagbackup = F)


  if applyonly != True: aU.plotbandpass(caltable='uid___A002_Xa8df68_Xb4b.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='92.1875%',showfdm=True, showBasebandNumber=True, showimage=False, 
    field='', figfile='uid___A002_Xa8df68_Xb4b.ms.tsys.plots.overlayTime/uid___A002_Xa8df68_Xb4b.ms.tsys') 
  


# Generation of the antenna position cal table
mystep = 6
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Position for antenna DV12 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV18 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA64 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA63 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA62 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV11 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA44 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA47 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA41 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA42 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV23 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV15 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV25 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV08 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV09 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV10 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV22 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV04 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA52 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA53 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA50 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA51 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV24 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA57 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA54 is derived from baseline run made on 2015-07-12 17:20:43.
  
  # Position for antenna DA55 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV06 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV07 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DA58 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV02 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV03 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV01 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV20 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV16 is derived from baseline run made on 2015-08-04 05:15:55.
  
  # Position for antenna DV13 is derived from baseline run made on 2015-08-04 05:15:55.
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.antpos') 
  gencal(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.antpos',
    caltype = 'antpos',
    antenna = 'DV20,DV18,DA64,DA63,DA62,DV11,DA44,DA47,DV12,DA41,DA50,DA42,DA51,DA57,DA54,DA55,DV07,DV10,DA58,DA52,DA53,DV22,DV23,DV24,DV25,DV08,DV09,DV06,DV15,DV04,DV02,DV03,DV01,DV16,DV13',
    parameter = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  #  parameter = [-0.000221839660157,0.000346910921173,0.000264363344571,-0.000128212173825,-0.000122380614286,0.000122464785674,-0.000207690850126,0.000103665931555,1.07650085352e-05,-0.000246277335521,0.000261395039717,0.000111126536973,-0.000289658688066,0.000105918890098,0.000191924550667,-0.000234128898797,0.000229045253644,-0.000364768324463,8.14211870312e-05,-0.000458046852059,8.47860343533e-05,0.000110554889203,-0.000307305869212,-0.000286443221636,-0.000374110897902,2.71041790439e-05,0.000314648056643,-1.0021229211e-05,6.0101771253e-06,5.49677426861e-05,-0.00025095781748,0.000156579918409,3.07555360296e-05,-0.000127083060314,0.000112575746756,0.000194570288131,-8.72450817554e-07,0.000229284510995,0.000162461847359,3.96474638441e-05,-0.000310785704733,0.000200015810236,-0.000211115402606,1.05159365076e-06,3.76575801814e-06,-0.000324716975517,0.0002867266026,3.15395176378e-05,-0.000345673877515,0.000118807595051,4.41511694324e-05,-6.65575907608e-05,-0.00015163124275,2.30942034569e-07,-0.000100381419225,-0.000116473821503,-8.91763397926e-05,-0.000126692168962,-0.000107596628631,0.000271293309595,-0.000238965118247,0.000119352480187,0.000177109168973,-0.000367546660154,0.000111422764918,0.000179411315328,-0.000101462197537,-0.000100607334197,3.44323507292e-05,-5.85658339886e-05,0.000267330256411,0.000199575993287,-0.00025508289008,0.000193107674288,0.000330066880071,-0.00014279138017,0.000206464449675,-4.84567102646e-05,-0.000110757963056,-4.61881493877e-05,9.1920228872e-05,7.16563922017e-05,9.71289968133e-06,-0.000296396729126,-0.000281386521826,0.000326110698093,7.14663412003e-05,-9.74708950745e-05,-7.60311018594e-05,-9.55442105173e-05,-0.000258917375863,-0.000267744559472,-2.31007218331e-05,-5.64765679256e-05,0.000169737680127,0.000176299509983,0.000102295532779,6.15939281393e-05,-9.21456021188e-05,-0.000135090947074,0.000182340344215,0.00021442582608,-0.000187695036972,5.21095375799e-05,0.000182800900407])
  

# Application of the WVR, Tsys and antpos cal tables
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  
  from recipes.almahelpers import tsysspwmap
  tsysmap = tsysspwmap(vis = 'uid___A002_Xa8df68_Xb4b.ms', tsystable = 'uid___A002_Xa8df68_Xb4b.ms.tsys', tsysChanTol = 1)
  
  
  
  applycal(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    field = '0',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_Xa8df68_Xb4b.ms.tsys', 'uid___A002_Xa8df68_Xb4b.ms.wvr', 'uid___A002_Xa8df68_Xb4b.ms.antpos'],
    gainfield = ['0', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1517-2422 didn't have any Tsys measurement, and I couldn't find any close measurement. But this is not a science target, so this is probably Ok.
  
  applycal(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    field = '2',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_Xa8df68_Xb4b.ms.tsys', 'uid___A002_Xa8df68_Xb4b.ms.wvr', 'uid___A002_Xa8df68_Xb4b.ms.antpos'],
    gainfield = ['2', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1922+1530 didn't have any Tsys measurement, so I used the one made on w51. This is probably Ok.
  
  applycal(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    field = '3',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_Xa8df68_Xb4b.ms.tsys', 'uid___A002_Xa8df68_Xb4b.ms.wvr', 'uid___A002_Xa8df68_Xb4b.ms.antpos'],
    gainfield = ['4', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    field = '4~40',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_Xa8df68_Xb4b.ms.tsys', 'uid___A002_Xa8df68_Xb4b.ms.wvr', 'uid___A002_Xa8df68_Xb4b.ms.antpos'],
    gainfield = ['4', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  if applyonly != True: es.getCalWeightStats('uid___A002_Xa8df68_Xb4b.ms') 
  

# Split out science SPWs and time average
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split') 
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.flagversions') 
  split(vis = 'uid___A002_Xa8df68_Xb4b.ms',
    outputvis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    datacolumn = 'corrected',
    spw = '17,19,21,23',
    keepflags = T)
  
  

print "# Calibration"

# Listobs, clear pointing table, and save original flags
mystep = 9
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.listobs')
  listobs(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    listfile = 'uid___A002_Xa8df68_Xb4b.ms.split.listobs')
  
  tb.open('uid___A002_Xa8df68_Xb4b.ms.split/POINTING', nomodify = False)
  a = tb.rownumbers()
  tb.removerows(a)
  tb.close()
  
  if not os.path.exists('uid___A002_Xa8df68_Xb4b.ms.split.flagversions/Original.flags'):
    flagmanager(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
      mode = 'save',
      versionname = 'Original')
  
  

# Initial flagging
mystep = 10
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Flagging shadowed data
  
  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'shadow',
    flagbackup = F)
  
  # Flagging atmospheric line(s)
  
  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'manual',
    spw = '2:0~2671',
    field = '2',
    flagbackup = F)
  
  # Flagging atmospheric line(s)
  
  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'manual',
    spw = '2:0~2671',
    field = '2',
    flagbackup = F)
  
  # Flags, leftover Tsys artefacts and atmospheric lines in calibrators
  # chans spw 2: 180~278 cannot be flagged since they overlap with the science targets CO(2-1) lines
  # The lines near 231.28 GHz cannot be flagged as they are on top of 13CS(5-4), affetcted channels: 2:1684~1926

  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'manual',
    spw = '3:1730~2096',
    flagbackup = F)


  flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'manual',
    spw = '2:1434~1492',
    antenna='DV18',
    field='0',
    flagbackup = F)

 # Flags on flux calibrator, lines which are not in the model
 
flagdata(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'manual',
    spw = '0:2188~2209',
    field = '2',
    flagbackup = F)


# Putting a model for the flux calibrator(s)
mystep = 11
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  setjy(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    field = '2', # Titan
    spw = '0,1,2,3',
    standard = 'Butler-JPL-Horizons 2012')
  
  if applyonly != True:
    os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.setjy.field*.png') 
    for i in ['2']:
      plotms(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
        xaxis = 'uvdist',
        yaxis = 'amp',
        ydatacolumn = 'model',
        field = str(i),
        spw = '0,1,2,3',
        avgchannel = '9999',
        coloraxis = 'spw',
        plotfile = 'uid___A002_Xa8df68_Xb4b.ms.split.setjy.field'+i+'.png')
  

# Save flags before bandpass cal
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]




# Save flags
  
  flagmanager(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'save',
    versionname = 'BeforeBandpassCalibration')
  
  

# Bandpass calibration
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.ap_pre_bandpass') 
  
  gaincal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.ap_pre_bandpass',
    field = '0', # J1733-1304
    spw = '0:1536~2304,1:1536~2304,2:1500~2304,3:1536~2304',
    scan = '1,2,4',
    solint = 'int',
    refant = 'DA55',
    calmode = 'p')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.split.ap_pre_bandpass', msName='uid___A002_Xa8df68_Xb4b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.bandpass') 
  bandpass(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.bandpass',
    field = '0', # J1733-1304
    scan = '1,2,4',
    solint = 'inf',
    combine = 'scan',
    refant = 'DA55',
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_Xa8df68_Xb4b.ms.split.ap_pre_bandpass')
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch') 
  
  bandpass(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch',
    field = '0', # J1733-1304
    scan = '1,2,4',
    solint = 'inf,20ch',
    combine = 'scan',
    refant = 'DA55',
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_Xa8df68_Xb4b.ms.split.ap_pre_bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch', msName='uid___A002_Xa8df68_Xb4b.ms.split', interactive=False) 
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.split.bandpass', msName='uid___A002_Xa8df68_Xb4b.ms.split', interactive=False) 
  

# Save flags before gain cal
mystep = 14
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'save',
    versionname = 'BeforeGainCalibration')
  
  

# Gain calibration
mystep = 15
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Note: the Solar system object used for flux calibration is highly resolved on some baselines.
  # Note: we will first determine the flux of the phase calibrator(s) on a subset of antennas.
  
  delmod('uid___A002_Xa8df68_Xb4b.ms.split',field='3')
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.phase_short_int') 
  gaincal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.phase_short_int',
    field = '2', # Titan
    selectdata = T,
    #antenna = 'DA41,DA46,DA50,DA53,DA55,DA59,DA61,DA63,DV01,DV04,DV06,DV08,DV09,DV16,DV18,DV24&',
    solint = 'int',
    refant = 'DA55',
    gaintype = 'G',
    calmode = 'p',
    uvrange = '0~450m', #Titan resolved (check amp model vs uvdist plot)
    gaintable = 'uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch')
  
  gaincal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.phase_short_int',
    field = '0,3', # J1733-1304,J1922+1530
    selectdata = T,
    solint = 'int',
    refant = 'DA55',
    gaintype = 'G',
    calmode = 'p',
    append = T,
    gaintable = 'uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.split.phase_short_int', msName='uid___A002_Xa8df68_Xb4b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.ampli_short_inf') 
  gaincal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.ampli_short_inf',
    field = '0,2,3', # J1733-1304,Titan,J1922+1530
    selectdata = T,
    solint = 'inf',
    refant = 'DA55',
    gaintype = 'G', #was T
    calmode = 'ap', # was a
    gaintable = ['uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch', 'uid___A002_Xa8df68_Xb4b.ms.split.phase_short_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.split.ampli_short_inf', msName='uid___A002_Xa8df68_Xb4b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.flux_short_inf') 
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.fluxscale') 
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_Xa8df68_Xb4b.ms.split.fluxscale')
  
  fluxscaleDict = fluxscale(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.ampli_short_inf',
    fluxtable = 'uid___A002_Xa8df68_Xb4b.ms.split.flux_short_inf',
    reference = '2') # Titan
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: es.fluxscale2(caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.ampli_short_inf', removeOutliers=True, msName='uid___A002_Xa8df68_Xb4b.ms', writeToFile=True, preavg=10000)
  
  f = open('uid___A002_Xa8df68_Xb4b.ms.split.fluxscale')
  fc = f.readlines()
  f.close()
  
  for phaseCalName in ['J1922+1530']:
    for i in range(len(fc)):
      if fc[i].find('Flux density for '+phaseCalName) != -1 and re.search('in SpW=[0-9]+(?: \(.*?\))? is: [0-9]+\.[0-9]+', fc[i], re.DOTALL|re.IGNORECASE) != None:
        line = (re.search('in SpW=[0-9]+(?: \(.*?\))? is: [0-9]+\.[0-9]+', fc[i], re.DOTALL|re.IGNORECASE)).group(0)
        spwId = (line.split('='))[1].split()[0]
        flux = float((line.split(':'))[1].split()[0])
        setjy(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
          field = phaseCalName.replace(';','*;').split(';')[0],
          spw = spwId,
          standard = 'manual',
          fluxdensity = [flux,0,0,0])
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.phase_int') 
  gaincal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.phase_int',
    field = '0,2,3', # J1733-1304,Titan,J1922+1530
    solint = 'int',
    refant = 'DA55',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.split.phase_int', msName='uid___A002_Xa8df68_Xb4b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.flux_inf') 
  gaincal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.flux_inf',
    field = '0,2,3', # J1733-1304,Titan,J1922+1530
    solint = 'inf',
    refant = 'DA55',
    gaintype = 'T',
    calmode = 'a',
    gaintable = ['uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch', 'uid___A002_Xa8df68_Xb4b.ms.split.phase_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.split.flux_inf', msName='uid___A002_Xa8df68_Xb4b.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.phase_inf') 
  gaincal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    caltable = 'uid___A002_Xa8df68_Xb4b.ms.split.phase_inf',
    field = '0,2,3', # J1733-1304,Titan,J1922+1530
    solint = 'inf',
    refant = 'DA55',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa8df68_Xb4b.ms.split.phase_inf', msName='uid___A002_Xa8df68_Xb4b.ms.split', interactive=False) 
  

# Save flags before applycal
mystep = 16
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    mode = 'save',
    versionname = 'BeforeApplycal')
  
  

# Application of the bandpass and gain cal tables
mystep = 17
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i in ['0', '2']: # J1733-1304,Titan
    applycal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
      field = str(i),
      gaintable = ['uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch', 'uid___A002_Xa8df68_Xb4b.ms.split.phase_int', 'uid___A002_Xa8df68_Xb4b.ms.split.flux_inf'],
      gainfield = ['', i, i],
      interp = 'linear,linear',
      calwt = T,
      flagbackup = F)
  
  applycal(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    field = '3,4~40', # w51
    gaintable = ['uid___A002_Xa8df68_Xb4b.ms.split.bandpass_smooth20ch', 'uid___A002_Xa8df68_Xb4b.ms.split.phase_inf', 'uid___A002_Xa8df68_Xb4b.ms.split.flux_inf'],
    gainfield = ['', '3', '3'], # J1922+1530
    interp = 'linear,linear',
    calwt = T,
    flagbackup = F)
  

# Split out corrected column
mystep = 18
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.cal') 
  os.system('rm -rf uid___A002_Xa8df68_Xb4b.ms.split.cal.flagversions') 
  split(vis = 'uid___A002_Xa8df68_Xb4b.ms.split',
    outputvis = 'uid___A002_Xa8df68_Xb4b.ms.split.cal',
    datacolumn = 'corrected',
    keepflags = T)
  
  

# Save flags after applycal
mystep = 19
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa8df68_Xb4b.ms.split.cal',
    mode = 'save',
    versionname = 'AfterApplycal')
  
  

