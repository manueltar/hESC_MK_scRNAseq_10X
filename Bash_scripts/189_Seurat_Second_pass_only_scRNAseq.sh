#!/bin/bash

eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

conda activate multiome_QC_DEF

sample_array=$(echo 'MCO_01373_3GEX,MCO_01374_3GEX')


a=($(echo "$sample_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
do

    sample_array_sel=$i
    echo "$sample_array_sel"

     preliminary_filtered=$(echo "/scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/""processing_outputs""/""$sample_array_sel""/""intermediate""/""preliminary_filtered.rds")
     path_processing_outputs=$(echo "/scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/""processing_outputs""/""$sample_array_sel""/")
     intermediate_dir=$(echo "/scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/""processing_outputs""/""$sample_array_sel""/""intermediate""/")
     premerge_dir=$(echo "/scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/""processing_outputs""/""$sample_array_sel""/""pre_merge""/")


     output_dir=$(echo "/scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/""processing_outputs""/")
     Log_files=$(echo "$output_dir""/""Log_files/")

   echo "$preliminary_filtered"
   echo "$path_processing_outputs"
   echo "$intermediate_dir"
   echo "$premerge_dir"
   

   mkdir -p $premerge_dir

    ### Seurat_second_pass

    type=$(echo "$sample_array_sel""_""Seurat_second_pass")
    outfile_Seurat_second_pass=$(echo "$Log_files""outfile_2_""$type"".log")
    touch $outfile_Seurat_second_pass
    echo -n "" > $outfile_Seurat_second_pass
    name_Seurat_second_pass=$(echo "$type""_job")

 
    Rscript_Seurat_second_pass=$(echo "$Rscripts_path""491_Seurat_second_pass_only_scRNAseq.R")

    sample_name=$sample_array_sel
    mem=$(echo "8192")
    processors=$(echo "12")
    total_memory=$(( mem * processors ))

    echo "$processors"
    echo "$total_memory"

 

    myjobid_Seurat_second_pass=$(sbatch --job-name $name_Seurat_second_pass --output=$outfile_Seurat_second_pass --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_Seurat_second_pass --sample_name $sample_name --preliminary_filtered $preliminary_filtered --path_processing_outputs $path_processing_outputs --intermediate_dir $intermediate_dir --premerge_dir $premerge_dir --processors $processors --total_memory $total_memory --type $type --out $output_dir")
    myjobid_seff_Seurat_second_pass=$(sbatch --dependency=afterany:$myjobid_Seurat_second_pass --open-mode=append --output=$outfile_Seurat_second_pass --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Seurat_second_pass >> $outfile_Seurat_second_pass")


done

conda deactivate
