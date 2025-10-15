--------------> Path backed up: /group/soranzo/manuel.tardaguila/2025_hESC_MK_SCRNAseq_10X/no_competition/
-------------- Path in scratch: /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/

# 1. cellranger mapping and counts


$ sbatch ~/Scripts/sbatch/13_Cell_ranger_count_for_GEX_libraries.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ MCO_01373_3GEX /group/so\
ranzo/paola.benaglio/MK_diff_scRNA10x/250528_A01481_0302_AHWNTLDSXC/fastq_raw/


$ sbatch ~/Scripts/sbatch/13_Cell_ranger_count_for_GEX_libraries.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ MCO_01374_3GEX /group/so\
ranzo/paola.benaglio/MK_diff_scRNA10x/250528_A01481_0302_AHWNTLDSXC/fastq_raw/

# 2. cellbender correction

$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ MCO_01373_3GEX

$ sbatch ~/Scripts/sbatch/7_CellBender_v_scratch.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ MCO_01374_3GEX

# 3. Align unmapped reads to reference with barcodes

$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v3.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ /group/soranzo/manuel.tardaguila/Multiom\
e/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01373_3GEX

$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v3.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ /group/soranzo/manuel.tardaguila/Multiom\
e/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa MCO_01374_3GEX

# 4. Filter and keep cells with a concordant barcode asignation from 3 different UMIs

$ bash ~/Scripts/Wraper_scripts/119_Filter_Larry_and_graphs_v3_Paolas_filtering.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_SCRNAseq_10X/no_competition\
/ count_and_filter /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/deconvolute_LARRY/

# 5. Seurat First pass

$ bash ~/Scripts/Wraper_scripts/188_Seurat_First_pass_only_scRNAseq.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ processing_outputs

# 6. Seurat Second pass

$ bash ~/Scripts/Wraper_scripts/189_Seurat_Second_pass_only_scRNAseq.sh

# 7. Merge samples in 1 object

$ bash ~/Scripts/Wraper_scripts/190_Seurat_merge_samples_only_scRNAseq.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ processing_outputs

# 8. QC

----> Jupyter notebook: Final_QC_in_the_merged_object.ipynb

# 9. Recluster and export h5ad for rpca

$ bash ~/Scripts/Wraper_scripts/191_Recluster_at_low_res_and_export_h5ad.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ processing_outpu\
ts

# 11. Cell typist prediction

----> Jupyter notebook: Cell_Typist_triple_prediction_cell_identity.ipynb
----> Jupyter notebook: mapping_cell_types.ipynb

# 12. Rpca pre genotyping

$ bash ~/Scripts/Wraper_scripts/192_RPCA_and_clustering.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ processing_outputs

----> Jupyter notebook:RPCA_graphs_ALL_SAMPLES.ipynb

# 13. targeted amplification of GEX 1: processing reads

$ bash ~/Scripts/Wraper_scripts/193_cellranger_alignment_of_targeted_amp_GEX.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ targeted_amp\
licon_GEX /group/soranzo/paola.benaglio/amplicon_seq_lympho_diff_0725/250718_A02059_0175_BHHHTWDSXF/fastq_trimmed/

# 14. targeted amplification of GEX 2: cellranger count of reads

$ sbatch ~/Scripts/sbatch/14_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/t\
argeted_amplicon_GEX/ MCO_01373_3GEX /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/targeted_amplicon_GEX/cellranger/

$ sbatch ~/Scripts/sbatch/14_Cell_ranger_count_for_CHEK2_targeted_amplicon_GEX_libraries.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/t\
argeted_amplicon_GEX/ MCO_01374_3GEX /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/targeted_amplicon_GEX/cellranger/

# 15. targeted amplification of GEX 3: cell bender correction

$ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competitio\
n/targeted_amplicon_GEX/ MCO_01373_3GEX

$ sbatch /home/manuel.tardaguila/Scripts/sbatch/15_cellbender_on_targeted_GEX_amplification.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competitio\
n/targeted_amplicon_GEX/ MCO_01374_3GEX

# 16. Assign barcodes

----> Jupyter notebook: notebook_to_assign_barcodes.ipynb

# 17. Rpca post_genotyping


$ bash ~/Scripts/Wraper_scripts/194_RPCA_and_clustering_post_genotyping.sh /scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/ processing_output\
s


# 18. Final cell annotation

----> Jupyter notebook: Post_G_final_cell_annotation.ipynb # Add UMAP marker genes and UMAP Sample ID

# 19. DE results accounting for time variation


"Conclusion: Use the Full Model
The simple model's results for both cell types are unreliable because they are based on an invalid statistical unit that pools two different biological stat\
es.

Your Full Model (∼time_point+Genotype) is the correct way to handle these samples because it uses the time_point factor to statistically control for the mas\
sive biological/technical differences caused by the cell count and maturation shift.

Trust the full model's result: A true Genotype effect is only revealed when the time-dependent noise is removed."


$ bash ~/Scripts/Wraper_scripts/196_SCRNA_10X_DE_per_identity_with_time_point.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_SCRNAseq_10X/no_competition/ \
DE_per_identity_with_time_point


