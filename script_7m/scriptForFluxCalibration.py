import re

#if re.search('^4.3.1', casadef.casa_version) == None:
# sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.3.1')


print "# Flux calibration of the data."

#                                       "J1922+1530"     0     218.37      0.14836      0.15330       2015-04-30T09:05:22     uid___A002_X9f852b_X134a.ms.split.cal
#                                       "J1922+1530"     0     218.37      0.16204      0.15330       2015-05-05T09:08:08     uid___A002_X9ff365_X2a83.ms.split.cal
#                                       "J1922+1530"     0     218.37      0.14980      0.15330       2015-05-30T05:47:27      uid___A002_Xa1e1bf_X290.ms.split.cal
#                                       "J1922+1530"     0     218.37      0.14901      0.15330       2015-05-30T07:19:41      uid___A002_Xa1e1bf_X6f1.ms.split.cal
#                                       "J1922+1530"     1     219.36      0.14302      0.14666       2015-04-30T09:05:22     uid___A002_X9f852b_X134a.ms.split.cal
#                                       "J1922+1530"     1     219.36      0.15837      0.14666       2015-05-05T09:08:08     uid___A002_X9ff365_X2a83.ms.split.cal
#                                       "J1922+1530"     1     219.36      0.14369      0.14666       2015-05-30T05:47:27      uid___A002_Xa1e1bf_X290.ms.split.cal
#                                       "J1922+1530"     1     219.36      0.14086      0.14666       2015-05-30T07:19:41      uid___A002_Xa1e1bf_X6f1.ms.split.cal
#                                       "J1922+1530"     2     231.37      0.13636      0.14115       2015-04-30T09:05:22     uid___A002_X9f852b_X134a.ms.split.cal
#                                       "J1922+1530"     2     231.37      0.14903      0.14115       2015-05-05T09:08:08     uid___A002_X9ff365_X2a83.ms.split.cal
#                                       "J1922+1530"     2     231.37      0.13088      0.14115       2015-05-30T05:47:27      uid___A002_Xa1e1bf_X290.ms.split.cal
#                                       "J1922+1530"     2     231.37      0.13424      0.14115       2015-05-30T07:19:41      uid___A002_Xa1e1bf_X6f1.ms.split.cal
#                                       "J1922+1530"     3     233.98      0.13645      0.14058       2015-04-30T09:05:22     uid___A002_X9f852b_X134a.ms.split.cal
#                                       "J1922+1530"     3     233.98      0.14676      0.14058       2015-05-05T09:08:08     uid___A002_X9ff365_X2a83.ms.split.cal
#                                       "J1922+1530"     3     233.98      0.12993      0.14058       2015-05-30T05:47:27      uid___A002_Xa1e1bf_X290.ms.split.cal
#                                       "J1922+1530"     3     233.98      0.13429      0.14058       2015-05-30T07:19:41      uid___A002_Xa1e1bf_X6f1.ms.split.cal

setjy(vis = 'uid___A002_X9f852b_X134a.ms.split.cal',
  field = 'J1922+1530',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.15330, 0, 0, 0])

setjy(vis = 'uid___A002_X9ff365_X2a83.ms.split.cal',
  field = 'J1922+1530',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.15330, 0, 0, 0])

setjy(vis = 'uid___A002_Xa1e1bf_X290.ms.split.cal',
  field = 'J1922+1530',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.15330, 0, 0, 0])

setjy(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal',
  field = 'J1922+1530',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.15330, 0, 0, 0])

setjy(vis = 'uid___A002_X9f852b_X134a.ms.split.cal',
  field = 'J1922+1530',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.14666, 0, 0, 0])

setjy(vis = 'uid___A002_X9ff365_X2a83.ms.split.cal',
  field = 'J1922+1530',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.14666, 0, 0, 0])

setjy(vis = 'uid___A002_Xa1e1bf_X290.ms.split.cal',
  field = 'J1922+1530',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.14666, 0, 0, 0])

setjy(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal',
  field = 'J1922+1530',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.14666, 0, 0, 0])

setjy(vis = 'uid___A002_X9f852b_X134a.ms.split.cal',
  field = 'J1922+1530',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.14115, 0, 0, 0])

setjy(vis = 'uid___A002_X9ff365_X2a83.ms.split.cal',
  field = 'J1922+1530',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.14115, 0, 0, 0])

setjy(vis = 'uid___A002_Xa1e1bf_X290.ms.split.cal',
  field = 'J1922+1530',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.14115, 0, 0, 0])

setjy(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal',
  field = 'J1922+1530',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.14115, 0, 0, 0])

setjy(vis = 'uid___A002_X9f852b_X134a.ms.split.cal',
  field = 'J1922+1530',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.14058, 0, 0, 0])

setjy(vis = 'uid___A002_X9ff365_X2a83.ms.split.cal',
  field = 'J1922+1530',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.14058, 0, 0, 0])

setjy(vis = 'uid___A002_Xa1e1bf_X290.ms.split.cal',
  field = 'J1922+1530',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.14058, 0, 0, 0])

setjy(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal',
  field = 'J1922+1530',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.14058, 0, 0, 0])

os.system('rm -rf uid___A002_X9f852b_X134a.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X9f852b_X134a.ms.split.cal',
  caltable = 'uid___A002_X9f852b_X134a.ms.split.cal.ampli_inf',
  field = 'J1922+1530',
  solint = 'inf',
  combine = 'scan',
  refant = '',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X9f852b_X134a.ms.split.cal',
  field = '3,4~16', # J1922+1530,w51
  gaintable = 'uid___A002_X9f852b_X134a.ms.split.cal.ampli_inf',
  gainfield = '3', # J1922+1530
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_X9ff365_X2a83.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X9ff365_X2a83.ms.split.cal',
  caltable = 'uid___A002_X9ff365_X2a83.ms.split.cal.ampli_inf',
  field = 'J1922+1530',
  solint = 'inf',
  combine = 'scan',
  refant = '',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X9ff365_X2a83.ms.split.cal',
  field = '3,4~16', # J1922+1530,w51
  gaintable = 'uid___A002_X9ff365_X2a83.ms.split.cal.ampli_inf',
  gainfield = '3', # J1922+1530
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_Xa1e1bf_X290.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_Xa1e1bf_X290.ms.split.cal',
  caltable = 'uid___A002_Xa1e1bf_X290.ms.split.cal.ampli_inf',
  field = 'J1922+1530',
  solint = 'inf',
  combine = 'scan',
  refant = '',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_Xa1e1bf_X290.ms.split.cal',
  field = '3,4~16', # J1922+1530,w51
  gaintable = 'uid___A002_Xa1e1bf_X290.ms.split.cal.ampli_inf',
  gainfield = '3', # J1922+1530
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_Xa1e1bf_X6f1.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal',
  caltable = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal.ampli_inf',
  field = 'J1922+1530',
  solint = 'inf',
  combine = 'scan',
  refant = '',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal',
  field = '3,4~16', # J1922+1530,w51
  gaintable = 'uid___A002_Xa1e1bf_X6f1.ms.split.cal.ampli_inf',
  gainfield = '3', # J1922+1530
  calwt = F,
  flagbackup = F)

print "# Concatenating the data."

concat(vis = ['uid___A002_X9f852b_X134a.ms.split.cal', 'uid___A002_X9ff365_X2a83.ms.split.cal', 'uid___A002_Xa1e1bf_X290.ms.split.cal', 'uid___A002_Xa1e1bf_X6f1.ms.split.cal'],
  concatvis = 'calibrated.ms')


