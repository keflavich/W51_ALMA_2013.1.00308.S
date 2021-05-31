#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail	
#SBATCH --ntasks=16
#SBATCH --mem=128gb # Job memory request PER NODE
#SBATCH --nodes=1 # exactly 1 node
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg
#SBATCH --output=w51_2017_b3_lb_e2_%j.log
#SBATCH --job-name=w51_2017_b3_lb_e2
#SBATCH --export=ALL

export FIELD_ID="W51"
export BAND_TO_IMAGE=B3
export LOGFILENAME="casa_log_w51e2lbcont_${FIELD_ID}_${BAND_TO_IMAGE}_12M_$(date +%Y-%m-%d_%H_%M_%S).log"

WORK_DIR='/orange/adamginsburg/w51/2017.1.00293.S/uvdata'
WORK_DIR='/orange/adamginsburg/w51/2017.1.00293.S/may2021_imaging'

export CASA=/orange/adamginsburg/casa/casa-release-5.6.0-60.el7/bin/casa
export CASA=/orange/adamginsburg/casa/casa-release-5.7.0-134.el7/bin/casa

module load git

which python
which git

git --version
echo $?

# not useing almaimf pipeline for this
#export ALMAIMF_ROOTDIR="/orange/adamginsburg/ALMA_IMF/reduction/reduction"
#cd ${ALMAIMF_ROOTDIR}
#python getversion.py
#export PYTHONPATH=$ALMAIMF_ROOTDIR

imaging_script=/orange/adamginsburg/w51/W51_ALMA_2013.1.00308.S/script_longbaseline_b3/script_for_imaging_selfcal.py

cd ${WORK_DIR}
echo ${WORK_DIR}

pycode="field='w51e2'; cleanmask='cleanmask_e2.crtf'"

#xvfb-run -d ${CASA} --nogui --nologger --logfile=split_${LOGFILENAME} -c "execfile('$ALMAIMF_ROOTDIR/split_windows.py')"
echo "Logfilename is ${LOGFILENAME}"
echo xvfb-run -d ${CASA} --nogui --nologger --logfile=${LOGFILENAME} -c "${pycode}; execfile('${imaging_script}')"

xvfb-run -d ${CASA} --nogui --nologger --logfile=${LOGFILENAME} -c "${pycode}; execfile('${imaging_script}')"
