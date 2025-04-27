#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --ntasks=32
#SBATCH --mem=256gb # Job memory request PER NODE
#SBATCH --nodes=1 # exactly 1 node
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=astronomy-dept-b
#SBATCH --account=astronomy-dept
#SBATCH --output=/red/adamginsburg/w51/logs/w51_line_cutouts_%A_%a.log
#SBATCH --job-name=w51_line_cutouts
#SBATCH --export=ALL

# Path to CASA
export CASA=/orange/adamginsburg/casa/casa-6.6.0-2-py3.8.el8/bin/casa

# Define log directory
export LOGDIR="/red/adamginsburg/w51/logs"

# Set working directory
WORK_DIR='/orange/adamginsburg/w51'
cd ${WORK_DIR}

# Parameter file containing all job configurations
PARAM_FILE="${LOGDIR}/line_cutout_params.txt"

# Check if parameter file exists
if [ ! -f "${PARAM_FILE}" ]; then
    echo "Error: Parameter file ${PARAM_FILE} not found!"
    echo "Please run submit_line_cutouts_job.sh first to generate the parameter file."
    exit 1
fi

# Read the parameters for this job
PARAMS=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $PARAM_FILE)
echo "Parameters for this job: $PARAMS"

# Extract the parameters
read -r _ SRC SPW ROBUST SUFFIX NITER <<< "$PARAMS"

# Create a logfile name with the job parameters
LOGFILENAME="${LOGDIR}/casa_log_w51_line_cutout_${SRC}_spw${SPW}_robust${ROBUST}_${SLURM_ARRAY_TASK_ID}.log"
echo "Logfilename is ${LOGFILENAME}"

# Path to our Python script
PYTHON_SCRIPT="$WORK_DIR/W51_ALMA_2013.1.00308.S/script_longbaseline_b3/scriptForImaging_line_cutouts_single.py"

# Command to pass to CASA
CASA_CMD="sourcename='$SRC'; spw=$SPW; robust=$ROBUST; suffix='$SUFFIX'; niter=$NITER; execfile('${PYTHON_SCRIPT}')"

# Print the command we're about to run
echo "Running: xvfb-run -d ${CASA} --nogui --nologger --logfile=${LOGFILENAME} -c \"${CASA_CMD}\""

# Run the command
xvfb-run -d ${CASA} --nogui --nologger --logfile=${LOGFILENAME} -c "${CASA_CMD}"

# Print job completion
echo "Job ${SLURM_ARRAY_TASK_ID} completed at $(date)"