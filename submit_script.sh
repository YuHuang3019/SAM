#!/bin/sh
##
#SBATCH --account=glab             # The account name for the job.
#SBATCH --job-name=LBA_fft         # The job name.
#SBATCH -c 24                      # The number of cpu cores to use.
#SBATCH --time=23:59:00            # The time the job will take to run.
#SBATCH --mem=128gb
#SBATCH --mail-type=ALL            # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=scs2229@columbia.edu

module load intel-parallel-studio/2017

mpirun -np 24 ./SAM_RAD_CAM_MICRO_M2005

date

# End of script







