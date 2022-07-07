#!/bin/bash
## Lines starting with #SBATCH are job options
## Lines starting with ## Are ignored by the job scheduler
#SBATCH -p all  
#SBATCH --job-name=fits-partdata
## Where job output is written to. CHANGE THIS TO YOUR OWN HOME DIR
#SBATCH --output=/home/clilje/summer_project/"Statistics"/"Job Results"/job_%j
## Resources requested
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
## Runtime of job. Job will be cancelled if it exceeds this time
#SBATCH --time=24:00:00
## Types of email notifications (NONE, BEGIN, END, FAIL, REQUEUE, ALL)
#SBATCH --mail-type=END
## CHANGE THIS TO YOUR OWN EMAIL ADDRESS
#SBATCH --mail-user=c.lilje@sms.ed.ac.uk
#################

cd ${SLURM_SUBMIT_DIR}
python -u graph.py
echo 'My job ran!!!'

