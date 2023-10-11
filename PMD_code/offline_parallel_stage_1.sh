#!/bin/bash -l

#SBATCH --job-name=stage_1_offline_parallel
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --output='slurm_output/slurm-%A_%a.out'


module load ${15}
cd $1

/usr/bin/time -v matlab -nodesktop -nodisplay -r "offline_parallel_stage_1('$1',$2,'$3','$4','$5','$6',$7,'$8','$9',${10},'${11}',${12},${13},'${14}');exit;"
wait
echo "All Processess are Complete"
