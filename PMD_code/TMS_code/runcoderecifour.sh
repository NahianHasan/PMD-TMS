#!/bin/bash -l
#SBATCH --array=1-941%300
#SBATCH --job-name=coilopt
#SBATCH -p wmglab
#SBATCH --cpus-per-task=4
#SBATCH --mem=63000
#SBATCH --time=5:00:00
#SBATCH --mail-user=ljg24@duke.edu
#SBATCH --export=ALL
echo "I ran on:" "${SLURM_ARRAY_TASK_ID}"
#
cd /work/ljg24/FEM_reciprocity/matlab
#
i=$((${SLURM_ARRAY_TASK_ID}))
export OMP_STACKSIZE=4096M     # omp default stack size is too small
ulimit -c 0                    # no core dumps
ulimit -s unlimited            # unlimited stack

/usr/bin/time -v /hpchome/apps/rhel7/matlabR2017b/bin/matlab -nodesktop -nodisplay -r "cd /work/ljg24/FEM_reciprocity/matlab;modegenerator($i,1);exit;"
/usr/bin/time -v /hpchome/apps/rhel7/matlabR2017b/bin/matlab -nodesktop -nodisplay -r "cd /work/ljg24/FEM_reciprocity/matlab;modegenerator($i,2);exit;"
/usr/bin/time -v /hpchome/apps/rhel7/matlabR2017b/bin/matlab -nodesktop -nodisplay -r "cd /work/ljg24/FEM_reciprocity/matlab;modegenerator($i,3);exit;"
