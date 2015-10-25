# ALMA Data Reduction Script
# $Id: scriptForPI.py,v 1.14 2015/02/06 09:32:30 dpetry Exp $

# Calibration application

import os
import sys
import glob

applyonly = True

savingslevel=0
if globals().has_key("SPACESAVING"):
    print 'SPACESAVING =', SPACESAVING
    if (type(SPACESAVING)!=int or SPACESAVING<0):
        sys.exit('ERROR: SPACESAVING value \"'+str(SPACESAVING)+'\" not permitted, must be int>0.\n'
                 +'Valid values: 0 = no saving,\n'
                 + '              1 = delete *.ms.split,\n'
                 + '              2 = delete *.ms and *.ms.split,\n'
                 + '            >=3 = delete *.ms, *.ms.split, and if possible *.ms.split.cal')
        
    savingslevel = SPACESAVING

if (os.path.basename(os.getcwd()) != 'script'):
    sys.exit('ERROR: Please start this script in directory \"script\".')

scriptnames = glob.glob('uid*.ms.scriptForCalibration.py')

pscriptnames = glob.glob('casa_piperestorescript.py')

pprnames = glob.glob('PPR*')

if ((len(scriptnames) + len(pscriptnames))  == 0):
    sys.exit('ERROR: No calibration script found.')

pprasdms = []
if (len(pprnames)>0):
    for line in open(pprnames[0]):
        if "<AsdmDiskName>" in line:
            pprasdms.append(line[line.index('uid'):line.index('</')])

try:
    os.chdir('../raw')
except:
    sys.exit('ERROR: directory \"raw\" not present.\n'
             '       Please download your raw data and unpack it to create and fill directory \"raw\".') 

# check available disk space
tmppipe = os.popen("df -P -m $PWD | awk '/[0-9]%/{print $(NF-2)}'")
avspace = int((tmppipe.readline()).rstrip('\n'))
tmppipe.close()
tmppipe = os.popen("du -sm ../../* | cut -f1")
packspace = int((tmppipe.readline()).rstrip('\n'))
tmppipe.close()

fcalpresent = 0
if os.path.exists('../script/scriptForFluxCalibration.py'):
    fcalpresent = 1

impreppresent = 0
if os.path.exists('../script/scriptForImagingPrep.py'):
    impreppresent = 1

spaceneed = packspace*(11.+fcalpresent*3.)

if (savingslevel==1):
    print 'Will delete intermediate MSs named *.ms.split to save disk space.'
    spaceneed = packspace*(7.+fcalpresent*3.)   
elif (savingslevel==2):
    print 'Will delete intermediate MSs named *.ms and *.ms.split to save disk space.'
    spaceneed = packspace*(3.+fcalpresent*3.)   
elif (savingslevel>=3):
    print 'Will delete all intermediate MSs to save disk space.'
    spaceneed = packspace*(3.+fcalpresent*3.)   


print 'Found ',avspace,' MB of available free disk space.'
print 'Expect to need up to ',spaceneed,' MB of free disk space.'
if(spaceneed>avspace):
    sys.exit('ERROR: not enough free disk space. Need at least '+str(spaceneed)+' MB.')

asdmnames = glob.glob('uid*.asdm.sdm')


if len(asdmnames) == 0:
    sys.exit('ERROR: No ASDM found in directory \"raw\".')

print 'Found the following ASDMs:', asdmnames

for i in range(len(asdmnames)):
    asdmnames[i] = asdmnames[i].replace('.asdm.sdm', '')


scriptasdms = []
for i in range(len(scriptnames)):
    scriptasdms.append(scriptnames[i].replace('.ms.scriptForCalibration.py', ''))

allasdms = []
allasdms.extend(scriptasdms)
allasdms.extend(pprasdms)

missing = []

if sorted(asdmnames) != sorted(allasdms):
    print "WARNING: Inconsistency between ASDMs and calibration scripts"
    print "         Calibration info available for: ", sorted(allasdms)
    print "         ASDMs available in directory raw: ", sorted(asdmnames)
    for myname in allasdms:
        if not (myname in asdmnames):
            missing.append(myname)
    if len(missing)==0:
        print "         Only the ASDMs for which there is calibration information will be calibrated."
    else:
        print "ERROR: the following ASDMs have calibration information but are absent from directory \"raw\":"
        print missing
        print "Will try to proceed with the rest ..."
        for myname in missing:
            if myname in scriptasdms:
                scriptasdms.remove(myname)
            if myname in pprasdms:
                pprasdms.remove(myname)
            if myname in allasdms:
                allasdms.remove(myname)
        if(len(allasdms)==0):
            sys.exit('ERROR: Nothing to process.')

ephnames = glob.glob('../calibration/*.eph')

if len(ephnames)>0:
    print "Note: this dataset uses external ephemerides."
    print "      You can find them in directory \"calibration\"."

if os.path.exists('../calibrated') and not globals().has_key("USEMS"):
    os.chdir('../calibrated')
    sys.exit('WARNING: will stop here since directory '+os.path.abspath(os.path.curdir)
             +' already exists.\nPlease delete it first and then try again.')
    
if not globals().has_key("USEMS"):
    print 'Creating destination directory for calibrated data.'
    os.mkdir('../calibrated')
else:
    print 'You have set USEMS. Will use your pre-imported MSs rather than importing them from the ASDMs.'
    for asdmname in scriptasdms:
        if not os.path.exists('../calibrated/'+asdmname+'.calibration/'+asdmname+'.ms'):
            print 'When USEMS is set, you must have created the directory \"calibrated\" and'
            print 'put the imported raw MSs \"uid*.ms\" in individual working directories'
            print 'named \"uid*.calibration\" inside \"calibrated\".'
            sys.exit('ERROR: cannot find calibrated/'+asdmname+'.calibration/'+asdmname+'.ms')
    
os.chdir('../calibrated')


for asdmname in scriptasdms:

    print 'Processing ASDM '+asdmname
    
    if not globals().has_key("USEMS"):
        os.mkdir(asdmname+'.calibration')

    os.chdir(asdmname+'.calibration')

    if not os.path.exists('../../raw/'+asdmname+'.asdm.sdm'):
        sys.exit('ERROR: cannot find raw/'+asdmname+'.asdm.sdm')

    os.system('ln -sf ../../raw/'+asdmname+'.asdm.sdm '+asdmname)

    for ephname in ephnames: 
        os.system('ln -sf ../'+ephname)

    execfile('../../script/'+asdmname+'.ms.scriptForCalibration.py')

    if not os.path.exists(asdmname+'.ms.split.cal'):
        print 'ERROR: '+asdmname+'.ms.split.cal was not created.'
    else:
        print asdmname+'.ms.split.cal was produced successfully, moving it to \"calibrated\" directory.'
        os.system('mv '+asdmname+'.ms.split.cal ..')
        if (savingslevel>=2):
            print 'Deleting intermediate MS ', asdmname+'.ms'
            os.system('rm -rf '+asdmname+'.ms')
        if (savingslevel>=1):
            print 'Deleting intermediate MS ', asdmname+'.ms.split'
            os.system('rm -rf '+asdmname+'.ms.split')

    os.chdir('..')

if (len(pprasdms)>0):

    print 'Processing the ASDMs ', pprasdms, ' using pipeline restore.'

    os.mkdir('rawdata')
    os.chdir('rawdata')
    for asdmname in pprasdms:
        if not os.path.exists('../../raw/'+asdmname+'.asdm.sdm'):
            sys.exit('ERROR: cannot find raw/'+asdmname+'.asdm.sdm')

        os.system('ln -sf ../../raw/'+asdmname+'.asdm.sdm '+asdmname)

    os.chdir('..')
        
    os.system('ln -sf ../calibration products')

    os.mkdir('working')
    os.chdir('working')

    print "now processing ", pscriptnames[0]

    execfile('../../script/'+pscriptnames[0])

    for asdmname in pprasdms:
        if not os.path.exists(asdmname+'.ms'):
            print 'ERROR: '+asdmname+'.ms was not created.'
        else:
            msmd.open(asdmname+'.ms')
            targetspws = msmd.spwsforintent('OBSERVE_TARGET*')
            sciencespws = ''
            outputspws = ''
            i = 0
            for myspw in targetspws:
                if msmd.nchan(myspw)>4:
                            sciencespws += str(myspw)+','
                            outputspws += str(i)+','
                            i += 1
            sciencespws = sciencespws.rstrip(',')
            outputspws = outputspws.rstrip(',')
            msmd.close()
            print 'Splitting out science SPWs for '+asdmname+': '+sciencespws+' -> '+outputspws
            split(vis=asdmname+'.ms', outputvis=asdmname+'.ms.split.cal', spw = sciencespws)
            if not os.path.exists(asdmname+'.ms.split.cal'):
                print 'ERROR: '+asdmname+'.ms.split.cal was not created.'
            else:
                if (savingslevel>=2):
                    print 'Deleting intermediate MS ', asdmname+'.ms'
                    os.system('rm -rf '+asdmname+'.ms')
                os.system('mv '+asdmname+'.ms.split.cal ..')

    os.chdir('..')

if fcalpresent > 0:
    print 'Executing scriptForFluxCalibration.py ...'
    execfile('../script/scriptForFluxCalibration.py')

if impreppresent > 0:
    print 'Executing scriptForImagingPrep.py ...'
    execfile('../script/scriptForImagingPrep.py')

if (savingslevel>=3) and os.path.exists('calibrated.ms'):
    for asdmname in allasdms:
        print 'Deleting intermediate MS ', asdmname+'.ms.split.cal'
        os.system('rm -rf '+asdmname+'.ms.split.cal')

print 'Done. Please find results in current directory.'
