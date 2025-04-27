#!/bin/bash

# Set working directory
WORK_DIR='/orange/adamginsburg/w51'
cd ${WORK_DIR}

# Create logs directory if it doesn't exist
LOGDIR="/red/adamginsburg/w51/logs"
mkdir -p ${LOGDIR}

# Create a parameter file
PARAM_FILE="${LOGDIR}/line_cutout_params.txt"
echo "Creating parameter file: ${PARAM_FILE}"
> ${PARAM_FILE}  # Clear/create the file

# Get all source names from source_ids.py
SOURCES=$(python3 -c "import sys; sys.path.append('$WORK_DIR/W51_ALMA_2013.1.00308.S/script_longbaseline_b3'); from source_ids import sources_fmtd; print(' '.join(sources_fmtd.keys()))")

# The parameters we will iterate over
SUFFIX="clarkclean1e5"
NITER=100000
ROBUST=0.5
SPW_INDICES="0 1 2 3"

# Create all combinations of parameters
INDEX=0
for SRC in $SOURCES; do
    for SPW in $SPW_INDICES; do
        echo "$INDEX $SRC $SPW $ROBUST $SUFFIX $NITER" >> ${PARAM_FILE}
        INDEX=$((INDEX+1))
    done
done

# Calculate the total number of jobs (zero-indexed)
TOTAL_JOBS=$((INDEX-1))
echo "Total number of jobs: $((TOTAL_JOBS+1))"
echo "Array range: 0-${TOTAL_JOBS}"

# Submit the job array with the calculated range
JOB_SCRIPT="${WORK_DIR}/W51_ALMA_2013.1.00308.S/script_longbaseline_b3/run_line_cutouts_array.sh"
echo "Submitting job array using script: ${JOB_SCRIPT}"
sbatch --array=0-${TOTAL_JOBS} ${JOB_SCRIPT}

echo "Job submission complete"