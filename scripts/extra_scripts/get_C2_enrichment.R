args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

library(BiocManager)
library(msigdbr)
library(ggplot2)
library(stats)
library(sets)
library(biomaRt)
library(clusterProfiler)
library(data.table)
library(stringi)
library(pheatmap)
library(tidyr)
library(ggpubr)
# Background set of genes
background_set <- fread("background_set.txt")

pathways <- readRDS("pathways.rds")

files <-list.files(args[1])
files <- files[grep("^GSE*", files)]
print(length(files))
compute_enrichment <- function(foreground_genes,all_genes,background_genes=NULL,pathways=NULL) {
    if (is.null(pathways)) { 
        pathways <- load_pathways()
        pathways <- lapply( pathways, function(pathway_genes) {return(pathway_genes[pathway_genes %in% all_genes])})
    }
    fisher_enrichment_dt <- data.table(pathway=names(pathways),p_value=-1,odds_ratio=-1)
    if (is.null(background_genes)) {
        background_genes <- setdiff(all_genes,foreground_genes)
        #background_genes <- unique(unlist(pathways))
    }

    for (pathway in names(pathways)) {
        pathway_genes <- pathways[[pathway]]
        non_pathway_genes <- setdiff(all_genes,pathway_genes)

        num_in_pathway_and_foreground <- intersect(pathway_genes,foreground_genes) %>% length
        num_in_pathway_and_not_foreground <- intersect(pathway_genes,background_genes) %>% length
        num_not_in_pathway_and_foreground <- intersect(non_pathway_genes,foreground_genes) %>% length
        num_not_in_pathway_and_not_foreground <- intersect(non_pathway_genes,background_genes) %>% length
        fisher_mat <- matrix(c(num_in_pathway_and_foreground,num_in_pathway_and_not_foreground,
                              num_not_in_pathway_and_foreground,num_not_in_pathway_and_not_foreground),
                             nrow=2,ncol=2,byrow=T)
        
        test_res <- fisher.test(fisher_mat,alternative="g")
        pathway_ <- pathway
        fisher_enrichment_dt[pathway==pathway_,`:=`(p_value=test_res$p.value, odds_ratio=test_res$estimate,
        num_p_fg=num_in_pathway_and_foreground,num_p_bg=num_in_pathway_and_not_foreground,
        num_not_p_fg=num_not_in_pathway_and_foreground,num_not_p_bg=num_not_in_pathway_and_not_foreground)]
    }
    fisher_enrichment_dt[,q_value:=p.adjust(p_value)]

    return(fisher_enrichment_dt)
}
get_enrichment_data <- function(current_files, current_pathways, current_dir){
    total_enrichment_pathways <- names(pathways[[current_pathways]])
    #pathwaysDF <- data.frame("pathways" = total_enrichment_pathways)
    for(i in 1:length(current_files)) {
        current_file_path <- paste(current_dir, current_files[i], sep="")
        new_file_name <- paste(current_dir,current_pathways,"_enrichment_",current_files[i], sep="")
        
        total_genes <- fread(current_file_path, header=FALSE)

        colnames(total_genes) <- c("genes")
        enrichment_test <- compute_enrichment(foreground_genes = total_genes$genes,
                                               all_genes = background_set$gene,
                                               pathways = pathways[[current_pathways]])
        subset_enrichment_test  <- enrichment_test[,c("q_value")]
        rownames(subset_enrichment_test) <-  enrichment_test$pathway

        #pathwaysDF <- cbind(pathwaysDF, subset_enrichment_test)
        saveRDS(subset_enrichment_test, file = new_file_name)
    }
    
    #all_qvals <- pathwaysDF[,2:length(pathwaysDF)]
    #rownames(all_qvals) <- total_enrichment_pathways
    #colnames(all_qvals) <- current_files[1:length(all_qvals)]
    
    #enrichment_heatmap <- pheatmap(as.matrix(log(all_qvals)), show_rownames = F)
    #return(list("matrix" = all_qvals, "heatmap" = enrichment_heatmap))
}

print(files[seq(args[2], args[3], 1)])
get_enrichment_data(files[seq(args[2], args[3], 1)], "C2", args[1])
