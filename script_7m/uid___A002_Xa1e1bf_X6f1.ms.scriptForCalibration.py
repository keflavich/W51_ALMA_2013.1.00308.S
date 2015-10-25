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


if re.search('^4.3.1', casadef.casa_version) == None:
 sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.3.1')


# CALIBRATE_AMPLI: Titan
# CALIBRATE_ATMOSPHERE: J1751+0939,Titan,w51
# CALIBRATE_BANDPASS: J1751+0939
# CALIBRATE_FLUX: Titan
# CALIBRATE_FOCUS: 
# CALIBRATE_PHASE: J1922+1530
# CALIBRATE_POINTING: J1517-2422,J1751+0939,J1922+1530
# OBSERVE_TARGET: w51

# Using reference antenna = 

# Import of the ASDM
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  if os.path.exists('uid___A002_Xa1e1bf_X6f1.ms') == False:
    importasdm('uid___A002_Xa1e1bf_X6f1', asis='Antenna Station Receiver Source CalAtmosphere CalWVR CorrelatorMode', bdfflags=False, lazy=False)
    os.system(os.environ['CASAPATH'].split()[0]+'/bin/bdflags2MS -f "COR DELA INT MIS SIG SYN TFB WVR ZER" uid___A002_Xa1e1bf_X6f1 uid___A002_Xa1e1bf_X6f1.ms')
  if applyonly != True: es.fixForCSV2555('uid___A002_Xa1e1bf_X6f1.ms')

# Fix of SYSCAL table times
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  from recipes.almahelpers import fixsyscaltimes
  fixsyscaltimes(vis = 'uid___A002_Xa1e1bf_X6f1.ms')

print "# A priori calibration"

# listobs
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.listobs')
  listobs(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    listfile = 'uid___A002_Xa1e1bf_X6f1.ms.listobs')
  
  

# A priori flagging
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  flagdata(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    mode = 'manual',
    spw = '0~23',
    autocorr = T,
    flagbackup = F)
  
  flagdata(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    mode = 'manual',
    intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*',
    flagbackup = F)
  
  flagcmd(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'plot',
    plotfile = 'uid___A002_Xa1e1bf_X6f1.ms.flagcmd.png')
  
  flagcmd(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'apply')
  

# Generation and time averaging of the WVR cal table
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  

# Generation of the Tsys cal table
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.tsys') 
  gencal(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.tsys',
    caltype = 'tsys')
  
  # Flagging edge channels
  
  flagdata(vis = 'uid___A002_Xa1e1bf_X6f1.ms.tsys',
    mode = 'manual',
    spw = '8:0~3;124~127,10:0~3;124~127,12:0~3;124~127,14:0~3;124~127',
    flagbackup = F)
  
  if applyonly != True: aU.plotbandpass(caltable='uid___A002_Xa1e1bf_X6f1.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='92.1875%',showfdm=True, showBasebandNumber=True, 
    field='', figfile='uid___A002_Xa1e1bf_X6f1.ms.tsys.plots.overlayTime/uid___A002_Xa1e1bf_X6f1.ms.tsys') 
  
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa1e1bf_X6f1.ms.tsys', msName='uid___A002_Xa1e1bf_X6f1.ms', interactive=False) 
  

# Generation of the antenna position cal table
mystep = 6
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Position for antenna CM02 is derived from baseline run made on 2014-06-02 04:15:36.
  
  # Position for antenna CM03 is derived from baseline run made on 2014-06-02 04:15:36.
  
  # Position for antenna CM12 is derived from baseline run made on 2015-04-05 02:00:43.
  
  # Position for antenna CM05 is derived from baseline run made on 2015-04-05 02:00:43.
  
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.antpos') 
  gencal(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.antpos',
    caltype = 'antpos',
    antenna = 'CM12,CM02,CM03,CM05',
    parameter = [0,0,0,0,0,0,0,0,0,0,0,0])
  #  parameter = [-1.26590310937e-05,0.000133968636039,-5.30862834409e-05,9.85723924856e-05,-0.000690061029368,-0.00037334201369,0.000443109311163,-0.000620041042566,-0.000481208786368,-8.22705465483e-05,0.000106050712869,1.76862508909e-05])
  

# Application of the WVR, Tsys and antpos cal tables
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  
  from recipes.almahelpers import tsysspwmap
  tsysmap = tsysspwmap(vis = 'uid___A002_Xa1e1bf_X6f1.ms', tsystable = 'uid___A002_Xa1e1bf_X6f1.ms.tsys', tsysChanTol = 1)
  
  
  
  applycal(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    field = '0',
    spw = '16,18,20,22',
    gaintable = ['uid___A002_Xa1e1bf_X6f1.ms.tsys', 'uid___A002_Xa1e1bf_X6f1.ms.antpos'],
    gainfield = ['0', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1517-2422 didn't have any Tsys measurement, and I couldn't find any close measurement. But this is not a science target, so this is probably Ok.
  
  applycal(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    field = '2',
    spw = '16,18,20,22',
    gaintable = ['uid___A002_Xa1e1bf_X6f1.ms.tsys', 'uid___A002_Xa1e1bf_X6f1.ms.antpos'],
    gainfield = ['2', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1922+1530 didn't have any Tsys measurement, so I used the one made on w51. This is probably Ok.
  
  applycal(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    field = '3',
    spw = '16,18,20,22',
    gaintable = ['uid___A002_Xa1e1bf_X6f1.ms.tsys', 'uid___A002_Xa1e1bf_X6f1.ms.antpos'],
    gainfield = ['4', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    field = '4~16',
    spw = '16,18,20,22',
    gaintable = ['uid___A002_Xa1e1bf_X6f1.ms.tsys', 'uid___A002_Xa1e1bf_X6f1.ms.antpos'],
    gainfield = ['4', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[]],
    calwt = T,
    flagbackup = F)
  
  if applyonly != True: es.getCalWeightStats('uid___A002_Xa1e1bf_X6f1.ms') 
  

# Split out science SPWs and time average
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split') 
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.flagversions') 
  split(vis = 'uid___A002_Xa1e1bf_X6f1.ms',
    outputvis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    datacolumn = 'corrected',
    spw = '16,18,20,22',
    keepflags = T)
  
  

print "# Calibration"

# Listobs, clear pointing table, and save original flags
mystep = 9
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.listobs')
  listobs(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    listfile = 'uid___A002_Xa1e1bf_X6f1.ms.split.listobs')
  
  tb.open('uid___A002_Xa1e1bf_X6f1.ms.split/POINTING', nomodify = False)
  a = tb.rownumbers()
  tb.removerows(a)
  tb.close()
  
  if not os.path.exists('uid___A002_Xa1e1bf_X6f1.ms.split.flagversions/Original.flags'):
    flagmanager(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
      mode = 'save',
      versionname = 'Original')
  
  

# Initial flagging
mystep = 10
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Flagging shadowed data
  
  flagdata(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    mode = 'shadow',
    flagbackup = F)
  
  # Flagging atmospheric line(s)
  
  flagdata(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    mode = 'manual',
    spw = '2:0~2740',
    field = '2',
    flagbackup = F)
  
  # Flagging atmospheric line(s)
  
  flagdata(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    mode = 'manual',
    spw = '1:0~104',
    field = '2',
    flagbackup = F)
  
  # Flagging atmospheric line(s)
  
  flagdata(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    mode = 'manual',
    spw = '2:0~2740',
    field = '2',
    flagbackup = F)
  
  # Flagging atmospheric line(s)
  
  flagdata(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    mode = 'manual',
    spw = '1:0~104',
    field = '2',
    flagbackup = F)
  
  

# Putting a model for the flux calibrator(s)
mystep = 11
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  setjy(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    field = '2', # Titan
    spw = '0,1,2,3',
    standard = 'Butler-JPL-Horizons 2012')
  
  if applyonly != True:
    os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.setjy.field*.png') 
    for i in ['2']:
      plotms(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
        xaxis = 'uvdist',
        yaxis = 'amp',
        ydatacolumn = 'model',
        field = str(i),
        spw = '0,1,2,3',
        avgchannel = '9999',
        coloraxis = 'spw',
        plotfile = 'uid___A002_Xa1e1bf_X6f1.ms.split.setjy.field'+i+'.png')
  

# Save flags before bandpass cal
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    mode = 'save',
    versionname = 'BeforeBandpassCalibration')
  
  

# Bandpass calibration
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.ap_pre_bandpass') 
  
  gaincal(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.ap_pre_bandpass',
    field = '0', # J1751+0939
    spw = '0:1638~2457,1:1632~2448,2:1632~2448,3:1632~2448',
    scan = '1,2,4',
    solint = 'int',
    refant = 'CM03',
    calmode = 'p')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa1e1bf_X6f1.ms.split.ap_pre_bandpass', msName='uid___A002_Xa1e1bf_X6f1.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.bandpass') 
  bandpass(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.bandpass',
    field = '0', # J1751+0939
    scan = '1,2,4',
    solint = 'inf',
    combine = 'scan',
    refant = 'CM03',
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_Xa1e1bf_X6f1.ms.split.ap_pre_bandpass')
  
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.bandpass_smooth20ch') 
  
  bandpass(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.bandpass_smooth20ch',
    field = '0', # J1751+0939
    scan = '1,2,4',
    solint = 'inf,20ch',
    combine = 'scan',
    refant = 'CM03',
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_Xa1e1bf_X6f1.ms.split.ap_pre_bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa1e1bf_X6f1.ms.split.bandpass_smooth20ch', msName='uid___A002_Xa1e1bf_X6f1.ms.split', interactive=False) 
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa1e1bf_X6f1.ms.split.bandpass', msName='uid___A002_Xa1e1bf_X6f1.ms.split', interactive=False) 
  

# Save flags before gain cal
mystep = 14
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    mode = 'save',
    versionname = 'BeforeGainCalibration')
  
  

# Gain calibration
mystep = 15
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.phase_int') 
  
  gaincal(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.phase_int',
    field = '0~0,2~3', # J1751+0939,Titan,J1922+1530
    solint = 'int',
    refant = 'CM03',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_Xa1e1bf_X6f1.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa1e1bf_X6f1.ms.split.phase_int', msName='uid___A002_Xa1e1bf_X6f1.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.ampli_inf') 
  gaincal(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.ampli_inf',
    field = '0~0,2~3', # J1751+0939,Titan,J1922+1530
    solint = 'inf',
    refant = 'CM03',
    gaintype = 'T',
    calmode = 'a',
    gaintable = ['uid___A002_Xa1e1bf_X6f1.ms.split.bandpass_smooth20ch', 'uid___A002_Xa1e1bf_X6f1.ms.split.phase_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa1e1bf_X6f1.ms.split.ampli_inf', msName='uid___A002_Xa1e1bf_X6f1.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.flux_inf') 
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.fluxscale') 
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_Xa1e1bf_X6f1.ms.split.fluxscale')
  
  fluxscaleDict = fluxscale(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.ampli_inf',
    fluxtable = 'uid___A002_Xa1e1bf_X6f1.ms.split.flux_inf',
    reference = '2') # Titan
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: es.fluxscale2(caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.ampli_inf', removeOutliers=True, msName='uid___A002_Xa1e1bf_X6f1.ms', writeToFile=True, preavg=10000)
  
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.phase_inf') 
  gaincal(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.phase_inf',
    field = '0~0,2~3', # J1751+0939,Titan,J1922+1530
    solint = 'inf',
    refant = 'CM03',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_Xa1e1bf_X6f1.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa1e1bf_X6f1.ms.split.phase_inf', msName='uid___A002_Xa1e1bf_X6f1.ms.split', interactive=False) 
  

# Save flags before applycal
mystep = 16
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    mode = 'save',
    versionname = 'BeforeApplycal')
  
  

# Application of the bandpass and gain cal tables
mystep = 17
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i in ['0', '2']: # J1751+0939,Titan
    applycal(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
      field = str(i),
      gaintable = ['uid___A002_Xa1e1bf_X6f1.ms.split.bandpass_smooth20ch', 'uid___A002_Xa1e1bf_X6f1.ms.split.phase_int', 'uid___A002_Xa1e1bf_X6f1.ms.split.flux_inf'],
      gainfield = ['', i, i],
      interp = 'linear,linear',
      calwt = T,
      flagbackup = F)
  
  applycal(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    field = '3,4~16', # w51
    gaintable = ['uid___A002_Xa1e1bf_X6f1.ms.split.bandpass_smooth20ch', 'uid___A002_Xa1e1bf_X6f1.ms.split.phase_inf', 'uid___A002_Xa1e1bf_X6f1.ms.split.flux_inf'],
    gainfield = ['', '3', '3'], # J1922+1530
    interp = 'linear,linear',
    calwt = T,
    flagbackup = F)
  

# Split out corrected column
mystep = 18
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.cal') 
  os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.cal.flagversions') 
  split(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split',
    outputvis = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal',
    datacolumn = 'corrected',
    keepflags = T)
  
  

# Save flags after applycal
mystep = 19
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal',
    mode = 'save',
    versionname = 'AfterApplycal')
  
  

