#!/bin/bash
# matlab batch submission script

#SBATCH --partition=campus-new     # more likely to get 12 cores in the full queue
#SBATCH --time=3-00:00:00        # 6 hours
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12       # request 16 cores
# #SBATCH --mail-user=youremail@fhcrc.org
#SBATCH --mail-type=END
#SBATCH --job-name=Matlab
#SBATCH --output=matlab-%j.out
#SBATCH --error=matlab-%j.err
#SBATCH --mem-per-cpu=10000
#SBATCH --qos=matlab

# grab environment variables for debugging
#env | grep SLURM

# set the path for the current version of matlab
. /app/lmod/lmod/init/profile
ml MATLAB/R2020a

matlab -nodisplay -nosplash -nodesktop -r "main_simulateManyWells;exit" 