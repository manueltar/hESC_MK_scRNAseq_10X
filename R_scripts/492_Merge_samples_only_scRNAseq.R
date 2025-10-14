
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


merge_and_filter_doublets = function(option_list)
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
  
  #### READ and transform sample_array ----
  
  sample_array = unlist(strsplit(opt$sample_array, split=","))
  
  cat("sample_array_\n")
  cat(sprintf(as.character(sample_array)))
  cat("\n")
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### LOOP TO READ in the pre merge object per sample and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  adatas <- list()
  

  
  for(i in 1:length(sample_array)){
    
    sample_array_sel<-sample_array[i]
    
    cat("---------------------------------------->\t") 
    cat(sprintf(as.character(sample_array_sel)))
    cat("\n")
    
    if(sample_array_sel == 'MCO_01373_3GEX' | sample_array_sel == 'MCO_01374_3GEX'){
      
      path_processing_outputs = paste('/scratch/manuel.tardaguila/hESC_MK_SCRNAseq_10X/no_competition/processing_outputs/',sample_array_sel,'/',sep='')
      
      Diff<-'Diff_MK_non_competition'
      
      
    }else{
      
     # Do nothing
      
    }#sample_array_sel == 'MCO_01373_3GEX'
    
    
    if (file.exists(path_processing_outputs)){
      
      
    }else{
      
      dir.create(file.path(path_processing_outputs))
      
    }#path_processing_outputs
    
    sample_dir <- path_processing_outputs
   
    
    pre_merge_dir = paste(path_processing_outputs,'pre_merge','/',sep='')
    
    if (file.exists(pre_merge_dir)){
      
      
    }else{
      
      dir.create(file.path(pre_merge_dir))
      
    }#pre_merge_dir
    
    setwd(pre_merge_dir)
    
    adata <- readRDS("pre_merged.rds")
    
    adata@meta.data$orig.ident = sample_array_sel
    adata@meta.data$Diff <- Diff
    adatas[[sample_array_sel]] <- adata
    DefaultAssay(adatas[[sample_array_sel]]) <- "RNA"

    
    
   
    
  }#i in 1:length(sample_array)
  
  
  merged = merge(x =adatas[[1]], y=adatas[2], add.cell.ids = sample_array )
  
  cat("merged_0\n")
  cat(str(merged))
  cat("\n")
  
  merged[["RNA"]] <-JoinLayers(merged[["RNA"]])
  
  cat("merged_1\n")
  cat(str(merged))
  cat("\n")
  
  rm(adatas)
  
  gc()
  
 
  
  setwd(out)
  
  saveRDS(merged, file = 'merged_unprocessed.rds')
  
 
  
  ##### Make a list of the genes corrected by CellBender (QC purposes) --------------
  

  celb_mat    = GetAssayData(merged, layer = "counts", assay = "RNA")
  raw_mat     = GetAssayData(merged, layer = "counts", assay = "RNA_raw")
  fract       = celb_mat/raw_mat
  means       = rowMeans(fract, na.rm=T)
  sorted_list = 1-sort(means)
  
  cat("fract\n")
  cat(str(fract))
  cat("\n")
  
  cat("means\n")
  cat(str(means))
  cat("\n")
  
  
  cat("sorted_list_CELL_BENDER_list\n")
  cat(str(sorted_list))
  cat("\n")
  
  Genes_corrected_more_than_0.25<-sorted_list[which(sorted_list >= 0.25)]
  
  cat("Genes_corrected_more_than_0.25\n")
  cat(str(Genes_corrected_more_than_0.25))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(Genes_corrected_more_than_0.25)))))
  cat("\n")
  cat(sprintf(as.character(summary(Genes_corrected_more_than_0.25))))
  cat("\n")
  
  ############ RMV cells with 0 counts after CellBender correction -----------------------------------------
  
  merged <- merged[, unname(which( colSums(GetAssayData(merged, slot = "counts", assay = "RNA"))!=0))]
  
  
  merged2<-subset(merged, scDblFinder.class != "doublet")
  

  

  
  ###### SAVE -----
  
  setwd(out)
  
  #### saveRDS(merged, file = 'merged_unprocessed.rds')
  
  saveRDS(merged2, file="merged_unprocessed_db_filt.rds")
  
  
  
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
    make_option(c("--sample_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--frag_file"), type="character", default=NULL, 
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
  
  merge_and_filter_doublets(opt)
 

}


###########################################################################

system.time( main() )
  