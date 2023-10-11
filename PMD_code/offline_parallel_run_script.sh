#!/bin/sh -l

#change the max wall time according to the specific cluster
stage_1_max_walltime=${20}
stage_2_max_walltime=${21}
stage_3_max_walltime=${22}
stage_4_max_walltime=${23}

#submit parallel jobs
while true; do
    DP=$(sbatch -A ${19} --cpus-per-task=${15} --time=$stage_1_max_walltime $1/offline_parallel_stage_1.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${24})
    if [ "$?" = "0" ]; then
		break
	else
		sleep 600
	fi
done
while true; do
    CR=$(sbatch -A ${19} --dependency=afterany:${DP##* } --array=[1-$2] --cpus-per-task=${16} --time=$stage_2_max_walltime $1/offline_parallel_stage_2.sh $2 $7 $8 $1 ${24})
    if [ "$?" = "0" ]; then
		break
	else
		sleep 600
	fi
done
while true; do
    CC=$(sbatch -A ${19} --dependency=afterany:${CR##* } --cpus-per-task=${17} --time=$stage_3_max_walltime $1/offline_parallel_stage_3.sh $2 $7 $8 $1 ${24})
    if [ "$?" = "0" ]; then
		break
	else
		sleep 600
	fi
done
while true; do
    RC=$(sbatch -A ${19} --dependency=afterany:${CC##* } --cpus-per-task=${18} --time=$stage_4_max_walltime --array=[1-$2] $1/offline_parallel_stage_4.sh $2 $7 $8 $1 ${24})
    if [ "$?" = "0" ]; then
		break
	else
		sleep 600
	fi
done
echo "Parallel jobs submitted"
