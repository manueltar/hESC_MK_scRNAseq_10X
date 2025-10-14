
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
suppressMessages(library(glmGamPoi))
suppressMessages(library(future))




opt = NULL

options(warn = -1)


log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}

create_the_premerged_Seurat_object = function(option_list)
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
  
  path_processing_outputs = opt$path_processing_outputs
  
  cat("path_processing_outputs_0\n")
  cat(sprintf(as.character(path_processing_outputs)))
  cat("\n")
  
  intermediate_dir = opt$intermediate_dir
  
  cat("intermediate_dir_0\n")
  cat(sprintf(as.character(intermediate_dir)))
  cat("\n")
  
  
  premerge_dir = opt$premerge_dir
  
  cat("premerge_dir_0\n")
  cat(sprintf(as.character(premerge_dir)))
  cat("\n")
  
  
  #### load preliminary_filtered seurat object ------------------------
  
  
  adata <- readRDS(file = opt$preliminary_filtered)
  
  # cat("adata_0\n")
  # cat(str(adata))
  # cat("\n")
  
 
  
  
  #### Read in cell bender counts -- filter them for previous bcs -- uses it as the main RNA modality ------------------
  
  cb = Read10X_h5(file.path(path_processing_outputs, 'cellbender_gex_seurat.h5'))
  
  if (is.list(cb) && length(cb) > 0) {
    cb_counts <- cb[[1]]
  } else {
    cb_counts <- cb
  }

  cb_counts   <- cb_counts[,colnames(adata)]
  
  cat("cb_counts_0\n")
  cat(str(cb_counts))
  cat("\n")
  
  adata2 = CreateSeuratObject(counts = cb_counts)
  
  adata2@meta.data$orig.ident = sample_name
  
  #add in previous raw RNA data as another assay (RNA_raw) for comparison
  
  DefaultAssay(adata) <- 'RNA'
  raw_rna <-  GetAssayData(object = adata, slot = "counts")
  raw_rna_assay <- CreateAssayObject(counts = raw_rna)
  adata2[['RNA_raw']] <- raw_rna_assay
  
  
  cat("adata2_1\n")
  # cat(str(adata2))
  cat("\n")
  
  adata2[['percent.mt']] <- PercentageFeatureSet(adata2, pattern = '^MT-')
  
  cat("adata2_0\n")
  # cat(str(adata2))
  cat("\n")
  
 
 
  ############ add previous GEX-only metadata ------------------------------------------
  
  old_meta = adata@meta.data
  
  # Only keep GEX-related doublet metrics
  colkeep = c('scDblFinder.class','scDblFinder.score','scDblFinder.weighted','scDblFinder.cxds_score')
  
  adata2@meta.data = cbind(adata2@meta.data,old_meta[,colkeep])
  
 
  # Filter out cells that somehow ended up with zero counts after CellBender
  adata2 <- adata2[, unname(which( colSums(GetAssayData(adata2, slot = "counts", assay = "RNA"))!=0))]
  
  ###### SAVE -----
  
  metadata_check<-adata2[[]]
  
  cat("metadata_check_0\n")
  cat(str(metadata_check))
  cat("\n")
  
  
  saveRDS(adata2, file = file.path(premerge_dir,'pre_merged.rds'))
  
  ###### SAVE -----
  
  
  
  saveRDS(adata2, file = file.path(premerge_dir,'pre_merged.rds'))
  
  
  
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
    make_option(c("--fragfile"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--path_processing_outputs"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--intermediate_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--snATAC_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--crange_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--premerge_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--preliminary_filtered"), type="character", default=NULL, 
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
  
  create_the_premerged_Seurat_object(opt)

}


###########################################################################

system.time( main() )
