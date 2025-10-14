
.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/lib/R/library"))
.libPaths()
# sessionInfo()

Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
library(reticulate)
reticulate::use_python("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
reticulate::use_condaenv("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')
suppressMessages(library("optparse"))
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(scDblFinder))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble"))
suppressMessages(library("biovizBase"))
suppressMessages(library("patchwork"))
suppressMessages(library(future))




opt = NULL

options(warn = -1)


log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}

create_the_preliminary_Seurat_object = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform processors ----
  
  processors = as.numeric(opt$processors)
  
  cat("processors\n")
  cat(sprintf(as.character(processors)))
  cat("\n")
  
  #### READ and transform memory ----
  
  #### READ and transform total_memory (memory in MB) ----
  total_memory = as.numeric(opt$total_memory) # Corrected variable name to match bash script
  cat("Total Memory (MB) for global objects:", as.character(total_memory), "\n") # Improved log message
  
  #### Assign resources -------------
  
  log_info_simple("plan stage")
  
  # Set up parallel processing: 'multiprocess' works on a single machine across cores.
  # 'total_memory' is expected in MB from the bash script, convert to bytes for future.globals.maxSize.
  plan("multicore", workers = processors)
  options(future.globals.maxSize = total_memory * 1024^2) # Corrected: Convert MB to bytes
  
  
  #### READ and transform sample_name ----
  
  sample_name = opt$sample_name
  
  cat("sample_name_\n")
  cat(sprintf(as.character(sample_name)))
  cat("\n")
  
  #### READ and transform master_path ----
  
  master_path = opt$master_path
  
  cat("master_path_\n")
  cat(sprintf(as.character(master_path)))
  cat("\n")
  
  sample_dir<-paste(master_path,sample_name,'/','outs','/',sep='')
  
  cat("sample_dir_\n")
  cat(sprintf(as.character(sample_dir)))
  cat("\n")
  
  #### READ and transform rna_min_features ----
  
  rna_min_features = opt$rna_min_features
  
  cat("rna_min_features_\n")
  cat(sprintf(as.character(rna_min_features)))
  cat("\n")
  
  # #### READ and transform atac_min_fragments ----
  # 
  # atac_min_fragments = opt$atac_min_fragments
  # 
  # cat("atac_min_fragments_\n")
  # cat(sprintf(as.character(atac_min_fragments)))
  # cat("\n")
  
  #### READ and transform MITO_max ----
  
  MITO_max = opt$MITO_max
  
  cat("MITO_max_\n")
  cat(sprintf(as.character(MITO_max)))
  cat("\n")
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  path_processing_outputs = paste(out,sample_name,'/',sep='')
  
  if (file.exists(path_processing_outputs)){
    
    
  }else{
    
    dir.create(file.path(path_processing_outputs))
    
  }#path_processing_outputs
  
  
  output_dir = paste(path_processing_outputs,'intermediate','/',sep='')
  
  if (file.exists(output_dir)){
    
    
  }else{
    
    dir.create(file.path(output_dir))
    
  }#output_dir
  
  
  
  
  ###### Read raw data -------
  
  inputdata.10x         <- Read10X_h5(file.path(sample_dir, 'raw_feature_bc_matrix.h5'))
  
  cat("inputdata.10x_0\n")
  cat(str(inputdata.10x))
  cat("\n")
  
  rna_counts            <- inputdata.10x

  cat("rna_counts_0\n")
  cat(str(rna_counts))
  cat("\n")
  
  ## Create object with rna_counts
  
  adata                 <- CreateSeuratObject(counts=rna_counts)
  
  cat("adata_0\n")
  cat(str(adata))
  cat("\n")
  
  metadata_adata<-adata[[]]
  
  cat("metadata_adata_0\n")
  cat(str(metadata_adata))
  cat("\n")
  
  ### Add the percent.mt for genes that start with MT (mithochondrial genes)
  
  adata[['percent.mt']] <- PercentageFeatureSet(adata, pattern = '^MT-')
  
  
  ##### Create the filter for cells that have equal or more than rna_min_features genes -----
  
  stored_filters = c() # Define the object to store the CB that pass the subsquent filters
  
  stored_filters['total_bc'] = length(colnames(adata[["RNA"]])) # All the cells of the adata RNA layer
  
  
  adata_sub <- subset(x = adata,subset = nFeature_RNA >= as.numeric(rna_min_features)) # subset for cells with more or equal to 500 genes
  
  cat("adata_sub_0\n")
  cat(str(adata_sub))
  cat("\n")
  
  stored_filters['after_rna_minfeat'] = length(colnames(adata_sub[["RNA"]])) # Only the cells with more or equal to 500 genes
  
  cat("stored_filters_0\n")
  cat(str(stored_filters))
  cat("\n")
  
  adata_sub[['percent.mt']] <- PercentageFeatureSet(adata_sub, pattern = '^MT-')
  
 
  
  metadata_adata_sub<-adata_sub[[]]
  
  cat("metadata_adata_sub_0\n")
  cat(str(metadata_adata_sub))
  cat("\n")
  

  
  ####### First doublet detection RNA with scDblFinder to be performed on lightly filtered data on the already filtered adata_sub object ------
  
  DefaultAssay(adata_sub) <- 'RNA'
  sce <- scDblFinder(GetAssayData(object = adata_sub, slot = "counts"))
  sce_results = data.frame(SummarizedExperiment::colData(sce))
  adata_sub@meta.data = cbind(adata_sub@meta.data,sce_results)
  
  
  metadata_adata_sub<-adata_sub[[]]
  
  cat("metadata_adata_sub_0\n")
  cat(str(metadata_adata_sub))
  cat("\n")
  
  

  
  
  ###### Remove high mito content cells (more or equal to 10% mt genes ------------------
  
  
  
  adata_sub_mito = subset(adata_sub, subset = percent.mt < as.numeric(MITO_max))
  
  stored_filters['after_mito'] = length(colnames(adata_sub_mito[["RNA"]]))
  
  cat("stored_filters_3\n")
  cat(str(stored_filters))
  cat("\n")
  
  ### Redefine adata as the object with all the previous filters implemented ---------------
  
  adata <- adata_sub_mito
  
  ### Compute intermediate bulk metrics of filtered object ----------------
  
 
  stored_filters['median_genesperCells_RNA'] = median(adata[[]][,'nFeature_RNA'])
  invisible(gc())
  
  
  #### Add sample name as the orig.ident field of the metadata ------
  
  adata@meta.data$orig.ident = sample_name
  
  
  metadata_adata<-adata[[]]
  
  cat("metadata_adata_REDEFINED\n")
  cat(str(metadata_adata))
  cat("\n")
  
 
  ########### PCA RNA analysis -------------------------
  
  cat("PCA RNA analysis\n")
  
  DefaultAssay(adata) <- 'RNA'
  suppressMessages(adata <- SCTransform(adata, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims=1:50, 
                                                                                         reduction.name='umap.rna', reduction.key='rnaUMAP_'))
  
  
  
  
  ################ Multimodal analysis using both RNA and ATAC layers for clusters ------------------------------------
  
  cat("Multimodal analysis using both RNA and ATAC layers for clusters\n")
  
  adata <- FindNeighbors(adata, reduction = 'pca', dims = 1:50)
  adata <- RunUMAP(adata, reduction = 'pca', dims = 1:50, reduction.name='umap', reduction.key='UMAP_')
  adata <- FindClusters(adata, algorithm=4, resolution = .5, verbose=FALSE) 
  
  

  
  #### SAVE object and selected barcodes --------------
  
  saveRDS(adata, file = file.path(output_dir,'preliminary_filtered.rds'))
  filtered_bcs <- colnames(adata[["RNA"]])
  write(filtered_bcs, file=(file.path(output_dir,'keep_barcodes_step1.txt')),sep='\n')
  
  
  #### output metrics ----------------------------------------------------------------------
  
  stored_filters["scDBL_RNA"] <- sum(adata@meta.data$scDblFinder.class == 'doublet')
  

  write.table(data.frame(stored_filters), (file.path(output_dir,'barcodes_stats.tsv')),
              sep="\t", quote=F, col.names=FALSE)
  
}



printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}



#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--sample_name"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--master_path"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--rna_min_features"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--atac_min_fragments"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MITO_max"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--processors"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--total_memory"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  create_the_preliminary_Seurat_object(opt)
  
 

}


###########################################################################

system.time( main() )
  