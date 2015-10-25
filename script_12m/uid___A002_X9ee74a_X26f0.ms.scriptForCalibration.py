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



# Using reference antenna = DA57 #DV18 contained some flags

# Import of the ASDM
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  if os.path.exists('uid___A002_X9ee74a_X26f0.ms') == False:
    importasdm('uid___A002_X9ee74a_X26f0', asis='Antenna Station Receiver Source CalAtmosphere CalWVR CorrelatorMode SBSummary', bdfflags=True, lazy=False, process_caldevice=False)
  if applyonly != True: es.fixForCSV2555('uid___A002_X9ee74a_X26f0.ms')

# Fix of SYSCAL table times
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  from recipes.almahelpers import fixsyscaltimes
  fixsyscaltimes(vis = 'uid___A002_X9ee74a_X26f0.ms')

print "# A priori calibration"

# listobs
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.listobs')
  listobs(vis = 'uid___A002_X9ee74a_X26f0.ms',
    listfile = 'uid___A002_X9ee74a_X26f0.ms.listobs')
  
  

# A priori flagging
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]



  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms',
    mode = 'manual',
    spw = '1~24',
    autocorr = T,
    flagbackup = F)
  
  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms',
    mode = 'manual',
    intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*',
    flagbackup = F)
  
  flagcmd(vis = 'uid___A002_X9ee74a_X26f0.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'plot',
    plotfile = 'uid___A002_X9ee74a_X26f0.ms.flagcmd.png')
  
  flagcmd(vis = 'uid___A002_X9ee74a_X26f0.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'apply')

# flagmanager(vis = 'uid___A002_Xa8df68_Xb4b.ms', mode = 'save', versionname = 'Apriori')


# Generation and time averaging of the WVR cal table
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.wvr') 
  
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.wvrgcal') 
  
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_X9ee74a_X26f0.ms.wvrgcal')
  
  wvrgcal(vis = 'uid___A002_X9ee74a_X26f0.ms',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.wvr',
    spw = [17, 19, 21, 23],
    smooth = '6.048s',
    toffset = 0,
    tie = ['w51,J1922+1530'],
    statsource = 'w51')
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: aU.plotWVRSolutions(caltable='uid___A002_X9ee74a_X26f0.ms.wvr', spw='17', antenna='DV18', #next time this should be plotted to DA57
    yrange=[-199,199],subplot=22, interactive=False,
    figfile='uid___A002_X9ee74a_X26f0.ms.wvr.plots/uid___A002_X9ee74a_X26f0.ms.wvr') 
  
  #Note: If you see wraps in these plots, try changing yrange or unwrap=True 
  #Note: If all plots look strange, it may be a bad WVR on the reference antenna.
  #      To check, you can set antenna='' to show all baselines.
  

# Generation of the Tsys cal table
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.tsys') 
  gencal(vis = 'uid___A002_X9ee74a_X26f0.ms',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.tsys',
    caltype = 'tsys')
  
  # Flagging edge channels
  
  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.tsys',
    mode = 'manual',
    spw = '9:0~3;124~127,11:0~3;124~127,13:0~3;124~127,15:0~3;124~127',
    flagbackup = F)

  if applyonly != True: aU.plotbandpass(caltable='uid___A002_X9ee74a_X26f0.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='92.1875%',showfdm=True, showBasebandNumber=True, showimage=False, 
    field='', figfile='uid___A002_X9ee74a_X26f0.ms.tsys.plots.overlayTimeBeforeTsysFlag/uid___A002_X9ee74a_X26f0.ms.tsys') 
  
  
 # Flagging ms.tsys with bad tsys calibration

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.tsys',
    mode = 'manual',
    spw = '15:10~20;26~36;41~51',
    antenna='DA45',
    flagbackup = F)

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.tsys',
    mode = 'manual',
    spw = '15:94~107',
    antenna='DA54',
    flagbackup = F)

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.tsys',
    mode = 'manual',
    spw = '9,11,13,15',
    antenna='DV15',
    flagbackup = F)

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.tsys',
    mode = 'manual',
    spw = '11:22~30,13:40~51',
    antenna='DV18',
    flagbackup = F)

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.tsys',
    mode = 'manual',
    spw = '13:4~11',
    antenna='',
    scan='10,15,20',
    flagbackup = F)



  if applyonly != True: aU.plotbandpass(caltable='uid___A002_X9ee74a_X26f0.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='92.1875%',showfdm=True, showBasebandNumber=True, showimage=False, 
    field='', figfile='uid___A002_X9ee74a_X26f0.ms.tsys.plots.overlayTime/uid___A002_X9ee74a_X26f0.ms.tsys') 
  
  
  if applyonly != True: es.checkCalTable('uid___A002_X9ee74a_X26f0.ms.tsys', msName='uid___A002_X9ee74a_X26f0.ms', interactive=False) 
  

# Generation of the antenna position cal table
mystep = 6
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Position for antenna DV12 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV18 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DA64 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Note: the correction for antenna DA49 is larger than 2mm.
  
  # Position for antenna DA49 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DA62 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DA61 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Note: the correction for antenna DA60 is larger than 2mm.
  
  # Position for antenna DA60 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DA45 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV10 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DA47 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DA46 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Note: the correction for antenna DA41 is larger than 2mm.
  
  # Position for antenna DA41 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DA43 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV16 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Note: the correction for antenna DA63 is larger than 2mm.
  
  # Position for antenna DA63 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV15 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV25 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DV08 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV11 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV19 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DA58 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV20 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Note: the correction for antenna DA53 is larger than 2mm.
  
  # Position for antenna DA53 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DA50 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DA51 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV24 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DA57 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DA54 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DV07 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DV04 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV05 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV02 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Position for antenna DV17 is derived from baseline run made on 2015-01-13 02:43:15.
  
  # Position for antenna DV13 is derived from baseline run made on 2014-12-21 05:53:17.
  
  # Note: no baseline run found for antenna DV21.
  
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.antpos') 
  gencal(vis = 'uid___A002_X9ee74a_X26f0.ms',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.antpos',
    caltype = 'antpos',
    antenna = 'DV19,DV18,DA64,DA49,DA62,DA61,DA60,DA45,DV10,DA47,DA46,DA41,DA43,DV16,DA51,DA63,DA57,DA54,DV11,DV07,DV04,DV20,DA53,DA50,DV12,DV24,DV25,DV08,DV15,DA58,DV05,DV02,DV17,DV13',
    parameter = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  #  parameter = [0.000763224595131,-0.000872836862611,-0.00100151686389,0.000568790779515,-0.000960048260564,-0.00102838046779,0.000376327428967,-0.00146782677621,-0.00112532451749,0.000575381331146,-0.00161136779934,-0.00117099285126,0.000670030074481,-0.00150872308554,-0.000943140154195,0.000350542975176,-0.00136158594062,-0.00102771511075,0.00102918455377,-0.00232664868236,-0.00126710766926,0.000471355393529,-0.00119987875223,-0.000884058885276,0.000297389547474,-0.00103880338778,-0.000898549931055,0.000490408856422,-0.00135660730302,-0.00102397520095,0.000566912652314,-0.000942276669298,-0.00089417600233,0.000858597457409,-0.00189473666251,-0.00125291012228,0.00067268848631,-0.00138315141508,-0.000907520447897,0.000569481262685,-0.00104891929093,-0.000820085003252,0.000748454127461,-0.00149133522063,-0.000856994651258,0.00189472222701,-0.00622595380992,-0.00312222354114,0.000514037961512,-0.00116958478532,-0.000836663256079,0.000382543263633,-0.00130245732194,-0.00111832579152,0.00037343815848,-0.0013356736041,-0.000845442467062,0.000829665135942,-0.00124130951181,-0.00117025657789,0.000513039702663,-0.00117243510926,-0.00120666779938,0.000449148240922,-0.00121965529975,-0.000785080128059,0.000661999452859,-0.00169800035655,-0.0011000004597,-0.000144578050822,-0.000582538545132,-0.000986936502159,0.000500225694647,-0.000744895727704,-0.000412322215341,0.00032130877439,-0.00116064453922,-0.000723131333211,0.000717123296325,-0.00148645733893,-0.000767197339289,0.000355394796746,-0.000718360463543,-0.00069032730126,0.000800814772724,-0.00122360248343,-0.00107297561125,0.000467471672268,-0.0012010330936,-0.00114111118815,0.000594357991325,-0.00133694808869,-0.000945342999283,0.000233996569865,-0.00104491981857,-0.000759456501754,0.000469437917569,-0.000701300406205,-0.000656311497542,0.000550707336515,-0.000622102059424,-0.000776093453169])
  

# Application of the WVR, Tsys and antpos cal tables
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  
  from recipes.almahelpers import tsysspwmap
  tsysmap = tsysspwmap(vis = 'uid___A002_X9ee74a_X26f0.ms', tsystable = 'uid___A002_X9ee74a_X26f0.ms.tsys', tsysChanTol = 1)
  
  
  
  applycal(vis = 'uid___A002_X9ee74a_X26f0.ms',
    field = '0',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X9ee74a_X26f0.ms.tsys', 'uid___A002_X9ee74a_X26f0.ms.wvr', 'uid___A002_X9ee74a_X26f0.ms.antpos'],
    gainfield = ['0', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1517-2422 didn't have any Tsys measurement, and I couldn't find any close measurement. But this is not a science target, so this is probably Ok.
  
  applycal(vis = 'uid___A002_X9ee74a_X26f0.ms',
    field = '2',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X9ee74a_X26f0.ms.tsys', 'uid___A002_X9ee74a_X26f0.ms.wvr', 'uid___A002_X9ee74a_X26f0.ms.antpos'],
    gainfield = ['2', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1922+1530 didn't have any Tsys measurement, so I used the one made on w51. This is probably Ok.
  
  applycal(vis = 'uid___A002_X9ee74a_X26f0.ms',
    field = '3',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X9ee74a_X26f0.ms.tsys', 'uid___A002_X9ee74a_X26f0.ms.wvr', 'uid___A002_X9ee74a_X26f0.ms.antpos'],
    gainfield = ['4', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_X9ee74a_X26f0.ms',
    field = '4~40',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X9ee74a_X26f0.ms.tsys', 'uid___A002_X9ee74a_X26f0.ms.wvr', 'uid___A002_X9ee74a_X26f0.ms.antpos'],
    gainfield = ['4', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  if applyonly != True: es.getCalWeightStats('uid___A002_X9ee74a_X26f0.ms') 
  

# Split out science SPWs and time average
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split') 
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.flagversions') 
  split(vis = 'uid___A002_X9ee74a_X26f0.ms',
    outputvis = 'uid___A002_X9ee74a_X26f0.ms.split',
    datacolumn = 'corrected',
    spw = '17,19,21,23',
    keepflags = T)
  
  

print "# Calibration"

# Listobs, clear pointing table, and save original flags
mystep = 9
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.listobs')
  listobs(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    listfile = 'uid___A002_X9ee74a_X26f0.ms.split.listobs')
  
  tb.open('uid___A002_X9ee74a_X26f0.ms.split/POINTING', nomodify = False)
  a = tb.rownumbers()
  tb.removerows(a)
  tb.close()
  
  if not os.path.exists('uid___A002_X9ee74a_X26f0.ms.split.flagversions/Original.flags'):
    flagmanager(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
      mode = 'save',
      versionname = 'Original')
  
  

# Initial flagging
mystep = 10
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Flagging shadowed data
  
  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    mode = 'shadow',
    flagbackup = F)
  
  # Flagging atmospheric line(s)
  
  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    mode = 'manual',
    spw = '2:0~2611',
    field = '2',
    flagbackup = F)
  
# flag the remaining tsys orginated problem, based on checks on the bandpass source, field 0

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms',
    mode = 'manual',
    spw = '3:380~387;892~897;1405~1410;1918~1921;2943~2945;3454~3457',
    antenna='DA45',
    flagbackup = F)

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms',
    mode = 'manual',
    spw = '3:3041~3230',
    antenna='DA54',
    flagbackup = F)

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms',
    mode = 'manual',
    spw = '1:569~913,2:1273~1483',
    antenna='DV18',
    flagbackup = F)
 
 #FLAGS on flux calibrator, lines which are absent in model

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
     mode = 'manual',
     spw = '0:2165~2185,1:2260~2280;2909~2915',
     field = '2',
     flagbackup = F)

# FLAGS on all sources due to atmospheric lines, based on checks on the phase calibrator, field 3 

# The lines 230.525-230.566 GHz (a few very narrow lines) cannot be flagged as they fall on the CO(2-1) line, affected channels: 2:188~278
# The lines near 231.28 GHz cannot be flagged as they are on top of 13CS(5-4), affetcted channels: 2:1569~1884

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
     mode = 'manual',
     spw = '3:1744~1937',
     field = '',
     flagbackup = F)

# Flag phase calibrator (field 3) for a part of scan 12, due to large phase scatter

  flagdata(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
      mode = 'manual',
      spw = '',
      field = '3',
      antenna = 'DA43',
      timerange='2015/04/23/09:54:00~2015/04/23/09:55:00',
      flagbackup = F)


# Putting a model for the flux calibrator(s)
mystep = 11
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  setjy(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    field = '2', # Titan
    spw = '0,1,2,3',
    standard = 'Butler-JPL-Horizons 2012')


 
  if applyonly != True:
    os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.setjy.field*.png') 
    for i in ['2']:
      plotms(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
        xaxis = 'uvdist',
        yaxis = 'amp',
        ydatacolumn = 'model',
        field = str(i),
        spw = '0,1,2,3',
        avgchannel = '9999',
        coloraxis = 'spw',
        plotfile = 'uid___A002_X9ee74a_X26f0.ms.split.setjy.field'+i+'.png')
  

# Save flags before bandpass cal
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]


  
  flagmanager(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    mode = 'save',
    versionname = 'BeforeBandpassCalibration',
    merge='replace')
  
 


# Bandpass calibration
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.ap_pre_bandpass') 
  
  gaincal(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.split.ap_pre_bandpass',
    field = '0', # J1751+0939
    spw = '0:1536~2304,1:1536~2304,2:1236~2304,3:1536~2304', #modified for spw 2, due to large abs bands near the center of the band
    scan = '1,2,4',
    solint = 'int',
    refant = 'DA57', #script suggested DV18
    calmode = 'p')
  
  if applyonly != True: es.checkCalTable('uid___A002_X9ee74a_X26f0.ms.split.ap_pre_bandpass', msName='uid___A002_X9ee74a_X26f0.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.bandpass') 
  bandpass(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.split.bandpass',
    field = '0', # J1751+0939
    scan = '1,2,4',
    solint = 'inf',
    combine = 'scan',
    refant = 'DA57', #script suggested DV18
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_X9ee74a_X26f0.ms.split.ap_pre_bandpass')
  
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.bandpass_smooth20ch') 
  
  bandpass(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.split.bandpass_smooth20ch',
    field = '0', # J1751+0939
    scan = '1,2,4',
    solint = 'inf,20ch',
    combine = 'scan',
    refant = 'DA57', #script suggested DV18
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_X9ee74a_X26f0.ms.split.ap_pre_bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_X9ee74a_X26f0.ms.split.bandpass_smooth20ch', msName='uid___A002_X9ee74a_X26f0.ms.split', interactive=False) 
  
  if applyonly != True: es.checkCalTable('uid___A002_X9ee74a_X26f0.ms.split.bandpass', msName='uid___A002_X9ee74a_X26f0.ms.split', interactive=False) 
  

# Save flags before gain cal
mystep = 14
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    mode = 'save',
    versionname = 'BeforeGainCalibration')
  
  

# Gain calibration
mystep = 15
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.phase_int') 
  
  gaincal(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.split.phase_int',
    field = '0~0,2~3', # J1751+0939,Titan,J1922+1530
    solint = 'int',
    refant = 'DA57', #script suggested DV18
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X9ee74a_X26f0.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_X9ee74a_X26f0.ms.split.phase_int', msName='uid___A002_X9ee74a_X26f0.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.ampli_inf') 
  gaincal(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.split.ampli_inf',
    field = '0~0,2~3', # J1751+0939,Titan,J1922+1530
    solint = 'inf',
    refant = 'DA57', #script suggested DV18
    gaintype = 'G', #was 'T'
    calmode = 'ap', # was 'a' 
    gaintable = ['uid___A002_X9ee74a_X26f0.ms.split.bandpass_smooth20ch', 'uid___A002_X9ee74a_X26f0.ms.split.phase_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_X9ee74a_X26f0.ms.split.ampli_inf', msName='uid___A002_X9ee74a_X26f0.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.flux_inf') 
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.fluxscale') 
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_X9ee74a_X26f0.ms.split.fluxscale')
  
  fluxscaleDict = fluxscale(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.split.ampli_inf',
    fluxtable = 'uid___A002_X9ee74a_X26f0.ms.split.flux_inf',
    reference = '2') # Titan
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: es.fluxscale2(caltable = 'uid___A002_X9ee74a_X26f0.ms.split.ampli_inf', removeOutliers=True, msName='uid___A002_X9ee74a_X26f0.ms', writeToFile=True, preavg=10000)
  
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.phase_inf') 
  gaincal(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    caltable = 'uid___A002_X9ee74a_X26f0.ms.split.phase_inf',
    field = '0~0,2~3', # J1751+0939,Titan,J1922+1530
    solint = 'inf',
    refant = 'DA57', #script suggested DV18
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X9ee74a_X26f0.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_X9ee74a_X26f0.ms.split.phase_inf', msName='uid___A002_X9ee74a_X26f0.ms.split', interactive=False) 
  

# Save flags before applycal
mystep = 16
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    mode = 'save',
    versionname = 'BeforeApplycal')
  
  

# Application of the bandpass and gain cal tables
mystep = 17
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i in ['0', '2']: # J1751+0939,Titan
    applycal(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
      field = str(i),
      gaintable = ['uid___A002_X9ee74a_X26f0.ms.split.bandpass_smooth20ch', 'uid___A002_X9ee74a_X26f0.ms.split.phase_int', 'uid___A002_X9ee74a_X26f0.ms.split.flux_inf'],
      gainfield = ['', i, i],
      interp = 'linear,linear',
      calwt = T,
      flagbackup = F)
  
  applycal(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    field = '3,4~40', # w51
    gaintable = ['uid___A002_X9ee74a_X26f0.ms.split.bandpass_smooth20ch', 'uid___A002_X9ee74a_X26f0.ms.split.phase_inf', 'uid___A002_X9ee74a_X26f0.ms.split.flux_inf'],
    gainfield = ['', '3', '3'], # J1922+1530
    interp = 'linear,linear',
    calwt = T,
    flagbackup = F)
  

# Split out corrected column
mystep = 18
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.cal') 
  os.system('rm -rf uid___A002_X9ee74a_X26f0.ms.split.cal.flagversions') 
  split(vis = 'uid___A002_X9ee74a_X26f0.ms.split',
    outputvis = 'uid___A002_X9ee74a_X26f0.ms.split.cal',
    datacolumn = 'corrected',
    keepflags = T)
  
  

# Save flags after applycal
mystep = 19
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X9ee74a_X26f0.ms.split.cal',
    mode = 'save',
    versionname = 'AfterApplycal')
  
  

