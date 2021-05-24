#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail	
#SBATCH --ntasks=16
#SBATCH --mem=96gb # Job memory request PER NODE
#SBATCH --nodes=1 # exactly 1 node
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg
#SBATCH --output=w51_2015_b6_lb_%j.log
#SBATCH --job-name=w51_2015_b6_lb
#SBATCH --export=ALL

export FIELD_ID="W51"
export BAND_TO_IMAGE=B6
export LOGFILENAME="casa_log_w51lbcont_${FIELD_ID}_${BAND_TO_IMAGE}_12M_$(date +%Y-%m-%d_%H_%M_%S).log"

WORK_DIR='/orange/adamginsburg/w51/w51-alma-longbaseline'

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

imaging_script=/orange/adamginsburg/w51/W51_ALMA_2013.1.00308.S/script_longbaseline_b6/script_longbaseline_big.py

cd ${WORK_DIR}
echo ${WORK_DIR}

#xvfb-run -d ${CASA} --nogui --nologger --logfile=split_${LOGFILENAME} -c "execfile('$ALMAIMF_ROOTDIR/split_windows.py')"
echo "Logfilename is ${LOGFILENAME}"
echo xvfb-run -d ${CASA} --nogui --nologger --logfile=${LOGFILENAME} -c "execfile('${imaging_script}')"

xvfb-run -d ${CASA} --nogui --nologger --logfile=${LOGFILENAME} -c "execfile('${imaging_script}')"
