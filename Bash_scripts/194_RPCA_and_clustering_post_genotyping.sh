#!/bin/bash

eval "$(conda shell.bash hook)"
 
  
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
Pythonscripts_path=$(echo "/home/manuel.tardaguila/Scripts/Python/")

MASTER_ROUTE=$1
analysis=$2


##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")
 
mkdir -p $Log_files
 
### POST_G_RPCA_and_recluster_at_

res_param=$(echo '0.5')

conda activate multiome_QC_DEF

type=$(echo "POST_G_RPCA_and_recluster_at_""$res_param")
outfile_POST_G_RPCA_and_recluster=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_POST_G_RPCA_and_recluster
echo -n "" > $outfile_POST_G_RPCA_and_recluster
name_POST_G_RPCA_and_recluster=$(echo "$type""_job")


Rscript_POST_G_RPCA_and_recluster=$(echo "$Rscripts_path""496_rpca_and_clustering_post_genotyping.R")
object_post_genotyping=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_QCed_cell_annotated_rpca_integrate_clustered_subcluster_majority_vote_only_genotyped.rds")
mem=$(echo "8192")
processors=$(echo "16")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"

# 15 8192
 
myjobid_POST_G_RPCA_and_recluster=$(sbatch --job-name $name_POST_G_RPCA_and_recluster --output=$outfile_POST_G_RPCA_and_recluster --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_POST_G_RPCA_and_recluster --object_post_genotyping $object_post_genotyping --res_param $res_param --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_POST_G_RPCA_and_recluster=$(sbatch --dependency=afterany:$myjobid_POST_G_RPCA_and_recluster --open-mode=append --output=$outfile_POST_G_RPCA_and_recluster --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_POST_G_RPCA_and_recluster >> $outfile_POST_G_RPCA_and_recluster")

conda deactivate

