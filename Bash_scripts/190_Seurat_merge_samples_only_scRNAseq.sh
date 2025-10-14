#!/bin/bash 
 
eval "$(conda shell.bash hook)"
  
   
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF


sample_array=$(echo 'MCO_01373_3GEX,MCO_01374_3GEX')
 
### Merge_pre_merged_per_sample

type=$(echo "$sample_array_sel""_""Merge_pre_merged_per_sample")
outfile_Merge_pre_merged_per_sample=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_Merge_pre_merged_per_sample
echo -n "" > $outfile_Merge_pre_merged_per_sample
name_Merge_pre_merged_per_sample=$(echo "$type""_job")


Rscript_Merge_pre_merged_per_sample=$(echo "$Rscripts_path""492_Merge_samples_only_scRNAseq.R")


mem=$(echo "8192")
processors=$(echo "12")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"


myjobid_Merge_pre_merged_per_sample=$(sbatch --job-name $name_Merge_pre_merged_per_sample --output=$outfile_Merge_pre_merged_per_sample --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_Merge_pre_merged_per_sample --sample_array $sample_array --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_Merge_pre_merged_per_sample=$(sbatch --dependency=afterany:$myjobid_Merge_pre_merged_per_sample --open-mode=append --output=$outfile_Merge_pre_merged_per_sample --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Merge_pre_merged_per_sample >> $outfile_Merge_pre_merged_per_sampley")

##########################################################################################################################

### cluster_merged_object

type=$(echo "cluster_merged_object")
outfile_cluster_merged_object=$(echo "$Log_files""outfile_4_""$type"".log")
touch $outfile_cluster_merged_object
echo -n "" > $outfile_cluster_merged_object
name_cluster_merged_object=$(echo "$type""_job")


Rscript_cluster_merged_object=$(echo "$Rscripts_path""493_Clustering_of_merged_samples_only_scRNAseq.R")

filtered_db_object=$(echo "$output_dir""merged_unprocessed_db_filt.rds")

mem=$(echo "8192")
processors=$(echo "12")
total_memory=$(( mem * processors ))

echo "$processors"
echo "$total_memory"


# --dependency=afterany:$myjobid_Merge_pre_merged_per_sample
 
myjobid_cluster_merged_object=$(sbatch --dependency=afterany:$myjobid_Merge_pre_merged_per_sample --job-name $name_cluster_merged_object --output=$outfile_cluster_merged_object --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_cluster_merged_object --filtered_db_object $filtered_db_object --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_cluster_merged_object=$(sbatch --dependency=afterany:$myjobid_cluster_merged_object --open-mode=append --output=$outfile_cluster_merged_object --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_cluster_merged_object >> $outfile_cluster_merged_object")


conda deactivate
