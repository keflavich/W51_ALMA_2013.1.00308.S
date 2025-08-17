#!/bin/bash

SCRIPT_ROOT='/orange/adamginsburg/w51'

# Set working directory
WORK_DIR='/red/adamginsburg/w51_imaging'
cd ${WORK_DIR}

# Create logs directory if it doesn't exist
LOGDIR="/red/adamginsburg/w51_imaging/logs"
mkdir -p ${LOGDIR}

export SAVEDIR='/orange/adamginsburg/w51/2015.1.01596.S/cubes/'

# Create a parameter file
PARAM_FILE="${LOGDIR}/line_cutout_params.txt"
echo "Creating parameter file: ${PARAM_FILE}"
> ${PARAM_FILE}  # Clear/create the file

# Get all source names from source_ids.py
SOURCES=$(python3 -c "import sys; sys.path.append('$SCRIPT_ROOT/W51_ALMA_2013.1.00308.S/script_longbaseline_b6'); from source_ids import sources_fmtd; print(' '.join(sources_fmtd.keys()))")

# The parameters we will iterate over
SUFFIX="clarkclean1e5"
NITER=100000
ROBUST=0.5
SPW_INDICES="0 1 2 3 4 5 6 7"


# Create all combinations of parameters
INDEX=0
for SRC in $SOURCES; do
    for SPW in $SPW_INDICES; do
        IMAGENAME="${SRC}_only.B6.robust${ROBUST}.spw${SPW}.${SUFFIX}.1536.bigpix.image.pbcor.fits"
        if [ -e $IMAGENAME ] || [ -e $SAVEDIR/$IMAGENAME ]; then
            echo "Image ${IMAGENAME} already exists, skipping"
        else
            # DEBUG  echo "$INDEX $SRC $SPW $ROBUST $SUFFIX $NITER because $IMAGENAME or $SAVEDIR/$IMAGENAME does not exist"
            echo "$INDEX $SRC $SPW $ROBUST $SUFFIX $NITER" >> ${PARAM_FILE}
            INDEX=$((INDEX+1))
        fi
    done
done

# Calculate the total number of jobs (zero-indexed)
TOTAL_JOBS=$((INDEX-1))
echo "Total number of jobs: $((TOTAL_JOBS+1))"
echo "Array range: 0-${TOTAL_JOBS}"

# Submit the job array with the calculated range
JOB_SCRIPT="${SCRIPT_ROOT}/W51_ALMA_2013.1.00308.S/script_longbaseline_b6/run_line_cutouts_array.sh"
echo "Submitting job array using script: ${JOB_SCRIPT}"
sbatch --array=0-${TOTAL_JOBS} ${JOB_SCRIPT}

echo "Job submission complete"