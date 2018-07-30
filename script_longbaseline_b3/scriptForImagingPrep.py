

######################################
# Splitting off the calibrated data


import glob

vislist = glob.glob('*[!_t].ms')  # match full ms, not target.ms

for myvis in vislist:

    msmd.open(myvis)
    targetspws = msmd.spwsforintent('OBSERVE_TARGET*')  
    sciencespws = []                                      
    for myspw in targetspws:                               
        if msmd.nchan(myspw)>4:
            sciencespws.append(myspw)
    sciencespws = ','.join(map(str,sciencespws))
    msmd.close()
    
    split(vis=myvis,outputvis=myvis+'.split.cal',spw=sciencespws)


########################################
# Getting a list of ms files to image

import glob

vislist=glob.glob('*.ms.split.cal')

##################################################
# Flag Bad Data [OPTIONAL]


# Save original flags
for vis in vislist:
    flagmanager(vis=vis,
                mode='save',
                versionname='original_flags')

# Inspect the science data
#fieldlist = ['3'] # list of science data fields to inspect
#spwlist = ['1'] # list of science spws to inspect

# loop through science data fields and spws to inspect.

#for vis in vislist:
#    for field in fieldlist:
#        for spw in spwlist:
#            plotms(vis=vis,xaxis='uvwave',yaxis='amp',avgtime='3e8',
#                   field=field,spw=spw) 
#            raw_input("push enter to continue")

#            plotms(vis=vis,xaxis='chan',yaxis='amp',avgtime='3e8',
#                   field=field,spw=spw) 
#            raw_input("push enter to continue")

# Flag the offending data. See flagdata help for more info.
#flagdata(vis='',mode='manual',action='apply',flagbackup=False)

# If you need to restore original flags, use the following command.
#flagmanager(vis='',mode='restore',versionname='original_flags')

###############################################################
# Combining Measurement Sets from Multiple Executions 


# If you have multiple executions, you will want to combine the
# scheduling blocks into a single ms using concat for ease of imaging
# and self-calibration. Each execution of the scheduling block will
# generate multiple spectral windows with different sky frequencies,
# but the same rest frequency, due to the motion of the Earth. Thus,
# the resulting concatentated file will contain n spws, where n is
# (#original science spws) x (number executions).  In other words, the
# multiple spws associated with a single rest frequency will not be
# regridded to a single spectral window in the ms.

concatvis='calibrated.ms'

rmtables(concatvis)
os.system('rm -rf ' + concatvis + '.flagversions')
concat(vis=vislist,
       concatvis=concatvis)

###################################
# Splitting off science target data

concatvis='calibrated.ms'


sourcevis='calibrated_source.ms'
rmtables(sourcevis)
os.system('rm -rf ' + sourcevis + '.flagversions')
split(vis=concatvis,
      intent='*TARGET*', # split off the target sources
      outputvis=sourcevis,
      datacolumn='data')

############################################
# Rename and backup data set

os.system('mv -i ' + sourcevis + ' ' + 'calibrated_final.ms')

# At this point you should create a backup of your final data set in
# case the ms you are working with gets corrupted by clean. 

os.system('cp -ir calibrated_final.ms calibrated_final.ms.backup')


############################################
# Output a listobs file

listobs(vis='calibrated_final.ms',listfile='calibrated_final.ms.listobs.txt') 
