if (!require("pacman")) install.packages("pacman")
pacman::p_load(biomaRt, BiocManager, stats, sets, stringi,
               ggplot2, tidyr, data.table, clusterProfiler,
               stringr, msigdbr, ggpubr, pheatmap) 

# go
go_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
unique_go_genes <- unique(go_gene_sets$gene_symbol)

# kegg
kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
unique_kegg_genes <- unique(kegg_gene_sets$gene_symbol)

# oncogenic

## hallmark
H_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
unique_H_genes <- unique(H_gene_sets$gene_symbol)

## C2
c2_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
unique_c2_genes <- unique(c2_gene_sets$gene_symbol)


## C4 CM
c4_gene_sets <- msigdbr(species = "Homo sapiens", category = "C4", subcategory = "CM")
unique_c4_genes <- unique(c4_gene_sets$gene_symbol)

## cosmic
cosmic_set <- fread("../data/gene_sets/cancer_gene_census.csv")
unique_cosmic_genes <- unique(cosmic_set$`Gene Symbol`)

## total
total_msigdbr <- msigdbr(species = "Homo sapiens")

# inflammation signature
print('###################################')
inflam_signature <- total_msigdbr$gene_symbol[grep("inflam",total_msigdbr$gs_name, ignore.case=T)]
print(paste("Total inflam Signature:",length(inflam_signature),sep=" "))
inflam_signature_union <- unique(inflam_signature)
print(paste("inflam Signature Union:",length(inflam_signature_union),sep=" "))
inflam_signature_recurrent <- sort(table(inflam_signature), decreasing = TRUE)
inflam_signature_recurrent <- inflam_signature_recurrent[inflam_signature_recurrent >1]
print(paste("inflam Signature Recurrent:",length(inflam_signature_recurrent),sep=" "))

# stress signature
print('###################################')
stress_signature <- total_msigdbr$gene_symbol[grep("stress",total_msigdbr$gs_name, ignore.case=T)]
print(paste("Total Stress Signature:",length(stress_signature),sep=" "))
stress_signature_union <- unique(stress_signature)
print(paste("Stress Signature Union:",length(stress_signature_union),sep=" "))
stress_signature_recurrent <- sort(table(stress_signature), decreasing = TRUE)
stress_signature_recurrent <- stress_signature_recurrent[stress_signature_recurrent >1]
print(paste("Stress Signature Recurrent:",length(stress_signature_recurrent),sep=" "))
print(paste("Stress Signature Unique:",length(unique(total_msigdbr$gs_name[grep("stress",total_msigdbr$gs_name, ignore.case=T)])), sep=" "))

# wound signature
print('###################################')
wound_signature <- total_msigdbr$gene_symbol[grep("wound",total_msigdbr$gs_name, ignore.case=T)]
print(paste("Total Wound Signature:",length(wound_signature),sep=" "))
wound_signature_union <- unique(wound_signature)
print(paste("Wound Signature Union:",length(wound_signature_union),sep=" "))
wound_signature_recurrent <- sort(table(wound_signature), decreasing = TRUE)
wound_signature_recurrent <- wound_signature_recurrent[wound_signature_recurrent >1]
print(paste("Wound Signature Recurrent:",length(wound_signature_recurrent),sep=" "))
print(paste("Wound Healing Signature Unique:",length(unique(total_msigdbr$gs_name[grep("wound",total_msigdbr$gs_name, ignore.case=T)])), sep=" "))

# regen signature
print('###################################')
regen_signature <- total_msigdbr$gene_symbol[grep("regen",total_msigdbr$gs_name, ignore.case=T)]
print(paste("Total Regen Signature:",length(regen_signature),sep=" "))
regen_signature_union <- unique(regen_signature)
print(paste("Regen Signature Union:",length(regen_signature_union),sep=" "))
regen_signature_recurrent <- sort(table(regen_signature), decreasing = TRUE)
regen_signature_recurrent <- regen_signature_recurrent[regen_signature_recurrent >1]
print(paste("Regen Signature Recurrent:",length(regen_signature_recurrent),sep=" "))
print(paste("Regeneration Signature Unique:",length(unique(total_msigdbr$gs_name[grep("regen",total_msigdbr$gs_name, ignore.case=T)])), sep=" "))

print('###################################')

# unions
# experiment type downregulated union -> exp_wrs_list[["downregulated"]]
# experiment type upregulated union -> exp_wrs_list[["upregulated"]]
combined_unions <- unlist(exp_wrs_list, recursive=F)
combined_unions[["msigdb.regen"]] <- regen_signature_union
combined_unions[["msigdb.wound"]] <- wound_signature_union
combined_unions[["msigdb.stress"]] <- stress_signature_union

# immediate unions
# experiment type downregulated union -> exp_wrs_list[["downregulated"]]
# experiment type upregulated union -> exp_wrs_list[["upregulated"]]
immediate_combined_unions <- unlist(exp_immediate_list, recursive=F)
immediate_combined_unions[["msigdb.regen"]] <- regen_signature_union
immediate_combined_unions[["msigdb.wound"]] <- wound_signature_union
immediate_combined_unions[["msigdb.stress"]] <- stress_signature_union

# experiment type upregulated reccurent -> recurrent_exp_wrs_list[["upregulated"]]
combined_recurrent <- unlist(recurrent_exp_wrs_list, recursive=F)

# reccurent sets
# experiment type downregulated recurrent -> recurrent_exp_wrs_list[["downregulated"]]
# experiment type upregulated reccurent -> recurrent_exp_wrs_list[["upregulated"]]
combined_recurrent <- unlist(recurrent_exp_wrs_list, recursive=F)
combined_recurrent[["msigdb.regen"]] <- names(regen_signature_recurrent)
combined_recurrent[["msigdb.wound"]] <- names(wound_signature_recurrent)
combined_recurrent[["msigdb.stress"]] <- names(stress_signature_recurrent)
#combined_recurrent[["msigdb.inflam"]] <- names(inflam_signature_recurrent)

# immediate reccurent sets
# experiment type downregulated recurrent -> recurrent_exp_wrs_list[["downregulated"]]
# experiment type upregulated reccurent -> recurrent_exp_wrs_list[["upregulated"]]
immediate_combined_recurrent <- unlist(recurrent_exp_immediate_list, recursive=F)
immediate_combined_recurrent[["msigdb.regen"]] <- names(regen_signature_recurrent)
immediate_combined_recurrent[["msigdb.wound"]] <- names(wound_signature_recurrent)
immediate_combined_recurrent[["msigdb.stress"]] <- names(stress_signature_recurrent)
#combined_recurrent[["msigdb.inflam"]] <- names(inflam_signature_recurrent)

# +
human <- readRDS("../data/biomart_orthologs/human.RDS")
 total_msigdbr <- msigdbr(species = "Homo sapiens")
 unique_total_msigdbr <- unique(total_msigdbr$gene_symbol)
 print(paste("Total H Gene Set Genes:", length(unique_total_msigdbr), sep = " "))
 total_msigdbr_sets_ensembl <- getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol",
                           values=unique_total_msigdbr,
                           mart=human,attributesL=c("hgnc_symbol"),
                           martL=human)
 print(paste("Total Available H Gene Set Genes:", nrow(total_msigdbr_sets_ensembl), sep = " "))

# background_set <- data.frame("gene"=total_msigdbr_sets_ensembl$HGNC.symbol)

# write.table(background_set,
#             file = "./background_set.txt", quote = FALSE, sep = "\t",
#             row.names = FALSE, col.names = TRUE)

# Background set of genes
#background_set <- fread("/data/timonaj/cancer_as_wound/background_set.txt")

# +
# pathway_list <- list("GO" = go_gene_sets[,c("gs_name", "gene_symbol")],
#                      "KEGG" = kegg_gene_sets[,c("gs_name", "gene_symbol")],
#                      "HALLMARK" = H_gene_sets[,c("gs_name", "gene_symbol")],
#                      "C2" = c2_gene_sets[,c("gs_name", "gene_symbol")],
#                      "C4" = c4_gene_sets[,c("gs_name", "gene_symbol")])
# pathways <- list()
# for(i in 1:length(pathway_list)) {
#     current_set <- pathway_list[[i]]
#     current_name <- names(pathway_list)[i]
#     pathway_name <- unique(current_set$gs_name)
#     pathways[[current_name]] <- list()
#     for(j in 1:length(pathway_name)) {
#         current_symbols <- current_set[current_set$gs_name == pathway_name[j],]$gene_symbol
#         final_symbols <- current_symbols[current_symbols  %in% background_set$gene]
#         pathways[[current_name]][[pathway_name[j]]] <- final_symbols
#     }
# }
# cosmic_genes <- fread("/data/timonaj/cancer_as_wound/cosmic_hallmark_genes.txt",header=F)%>% as.data.frame()
# categories <- unique(cosmic_genes$V2)
# cosmic <- list()
# for(i in 1:length(categories)) {
#     cosmic[[categories[i]]] <- cosmic_genes[cosmic_genes$V2 ==categories[i],]$V1
# }
# pathways[["Cosmic"]] <- cosmic
# pathways[["COSMIC"]] <- list("COSMIC" = cosmic_genes$V1)
# saveRDS(pathways, file = "pathways.rds")

pathways <- readRDS("../data/pathways.rds")
# -

# do enrichment analysis with 2 sets of gene
get_overlap_data <- function(current_files,overlap_files,background_gene_list=NULL){


 
    total_enrichment_pathways <- names(current_files)

    total_enrichment_overlap <- names(overlap_files)
    
    # dfs to return
    pathwaysDF <- data.frame("pathways" = total_enrichment_overlap)
    pathwaysDF_odds <- data.frame("pathways" = total_enrichment_overlap)

    for(i in 1:length(current_files)) {
        spec_exptype <- names(current_files)[i]
        overlap_vec <- numeric(length(overlap_files))
        fisher_enrichment_dt <- data.table(comparison=names(overlap_files),p_value=-1,odds_ratio=-1)
        
        for(j in 1:length(overlap_files)) {
            #set up matrix
            
            if (is.null(background_gene_list)) {
                background_genes <- background_set$gene
            }
            else {
                background_genes <- background_gene_list[[j]]
            }
            
            foreground_genes <- current_files[[i]]
            comparison_genes <- overlap_files[[j]]
            non_comparison_genes <- setdiff(background_genes,comparison_genes)
            non_foreground_genes <- setdiff(background_genes,foreground_genes)
            
            num_in_comparison_and_foreground <- intersect(comparison_genes,foreground_genes) %>% length
            num_in_comparison_and_not_foreground <- intersect(comparison_genes,non_foreground_genes) %>% length
            num_not_in_comparison_and_foreground <- intersect(non_comparison_genes,foreground_genes) %>% length
            num_not_in_comparison_and_not_foreground <- intersect(non_comparison_genes,non_foreground_genes) %>% length
            
            fisher_mat <- matrix(c(num_in_comparison_and_foreground,num_in_comparison_and_not_foreground,
                                  num_not_in_comparison_and_foreground,num_not_in_comparison_and_not_foreground),
                                 nrow=2,ncol=2,byrow=T)
            test_res <- fisher.test(fisher_mat,alternative="g")
            
            # 
            comparison_ <- names(overlap_files)[j]
            fisher_enrichment_dt[comparison==comparison_,`:=`(p_value=test_res$p.value,
                                                              odds_ratio=test_res$estimate,
                                                              num_p_fg=num_in_comparison_and_foreground,
                                                              num_p_bg=num_in_comparison_and_not_foreground,
                                                              num_not_p_fg=num_not_in_comparison_and_foreground,
                                                              num_not_p_bg=num_not_in_comparison_and_not_foreground)]

        }
        fisher_enrichment_dt[,q_value:=p.adjust(p_value)]
        pathwaysDF <- cbind(pathwaysDF, fisher_enrichment_dt$q_value)
        pathwaysDF_odds <- cbind(pathwaysDF_odds, fisher_enrichment_dt$odds_ratio)
        #saveRDS(subset_enrichment_test, file = new_file_name)
    }
    
    all_qvals <- pathwaysDF[,2:length(pathwaysDF)]
    colnames(all_qvals) <- names(current_files)
    rownames(all_qvals) <- names(overlap_files)
    
    all_odds <- pathwaysDF_odds[,2:length(pathwaysDF_odds)]
    colnames(all_odds) <- names(current_files)
    rownames(all_odds) <- names(overlap_files)
    
    #enrichment_heatmap <- pheatmap(as.matrix(all_qvals),
    #                               fontsize = 8)
    #enrichment_heatmap_unclustered <- pheatmap(as.matrix(all_qvals),
    #                               fontsize = 8,cluster_rows=FALSE,
    #                                    cluster_cols=FALSE)
    return(list("matrix" = all_qvals, "odds_matrix" = all_odds))
}

# +
# funtion to color non significant values white
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
change_matrix <- function(fdr_mat, odds_mat, threshhold) {
    new_mat <- matrix(nrow = nrow(fdr_mat), ncol = ncol(fdr_mat))
    for(i in 1:nrow(fdr_mat)) {
        for(j in 1:ncol(fdr_mat)) {
            
            if(fdr_mat[i,j] < .05) {
                new_mat[i,j] <- odds_mat[i,j]
                if(odds_mat[i,j] == "Inf") { 
                    new_mat[i,j] <- threshhold
                }
                
            }else {
                new_mat[i,j] <- threshhold
            }
        }
    }

    
    rownames(new_mat) <- rownames(fdr_mat)
    colnames(new_mat) <- colnames(fdr_mat)
    
#     cutoff.distance <-  0
#     cols <- makeColorRampPalette(c("white", "white",    # distances -1 to 0 colored from white to red
#                                    "lightblue", "navy"), # distances 1 to max(distmat) colored from green to black
#                                  cutoff.distance / max(new_mat),
#                                  500)

    new_heatmap <- pheatmap(new_mat,
                            color = colorRampPalette(c("white", "white", "white", "white", "white","white","lightblue2", "lightblue3", "lightblue4"))(10),
                            fontsize = 14,
                            breaks = c(-1, 0, 0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4))
    new_heatmap_unclustered <- pheatmap(new_mat,
                                        color = colorRampPalette(c("white", "white", "white", "white", "white", "white","lightblue2", "lightblue3", "lightblue4"))(10),
                                        fontsize = 14,
                                        breaks = c(-1, 0,  0.2, 0.4, 0.6, 0.8, 1, 2, 3, 4),
                                        cluster_rows=FALSE,
                                        cluster_cols=FALSE)
    
    return(list("matrix"=new_mat, "heatmap" = new_heatmap, "unclustered_heatmap" = new_heatmap_unclustered))                                      
# -

}

# +
change_matrix_odds <- function(fdr_mat, odds_mat, threshhold) {
    new_mat <- matrix(nrow = nrow(fdr_mat), ncol = ncol(fdr_mat))
    for(i in 1:nrow(fdr_mat)) {
        for(j in 1:ncol(fdr_mat)) {
            
            new_mat[i,j] <- odds_mat[i,j]
            if(odds_mat[i,j] == "Inf") { 
                new_mat[i,j] <- threshhold
            }
        }
    }

    
    rownames(new_mat) <- rownames(fdr_mat)
    colnames(new_mat) <- colnames(fdr_mat)
    
#     cutoff.distance <-  0
#     cols <- makeColorRampPalette(c("white", "white",    # distances -1 to 0 colored from white to red
#                                    "lightblue", "navy"), # distances 1 to max(distmat) colored from green to black
#                                  cutoff.distance / max(new_mat),
#                                  500)

    new_heatmap <- pheatmap(new_mat,
                            color = colorRampPalette(c("white", "white", "white", "white", "white",'lightblue1',"lightblue2", "lightblue3", "lightblue4"))(10),
                            fontsize = 8,
                            breaks = c(-1, 0.2, 0.4, 0.6, 0.8, 0, 1, 2, 3, 4))
    new_heatmap_unclustered <- pheatmap(new_mat,
                                        color = colorRampPalette(c("white", "white", "white", "white", "white",'lightblue1',"lightblue2", "lightblue3", "lightblue4"))(10),
                                        fontsize = 8,
                                        breaks = c(-1, 0.2, 0.4, 0.6, 0.8, 0, 1, 2, 3, 4),
                                        cluster_rows=FALSE,
                                        cluster_cols=FALSE)
    
    return(list("matrix"=new_mat, "heatmap" = new_heatmap, "unclustered_heatmap" = new_heatmap_unclustered))  
}
# -

get_reccurent_genes <- function(duplicated_list,species_exptype,total_exps,title,length) {
    print(paste(species_exptype,title, sep =" "))
    
    if(total_exps == 1) {
        print(paste(species_exptype, "total_datasets :", total_exps, sep =" "))
        print("################################################")
        return(duplicated_list)
    }
    
    duplicated_list <- sort(table(duplicated_list), decreasing = TRUE)
    
    if(length(duplicated_list) > length) {
        duplicated_list <- duplicated_list[1:length]
    }
    
    duplicated_list <- duplicated_list[duplicated_list > 1]
    
    print(paste(species_exptype, "total_datasets :", total_exps, sep =" "))
    print(paste(species_exptype, "mean reccurence :", mean(duplicated_list), sep =" "))
    print(paste(species_exptype, "min reccurence :", min(duplicated_list), sep =" "))
    print(paste(species_exptype, "max reccurence :", max(duplicated_list), sep =" "))
    print("################################################")
    return(names(duplicated_list))
}

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
    pathwaysDF <- data.frame("pathways" = total_enrichment_pathways)
    for(i in 1:length(current_files)) {
        current_file_path <- paste(current_dir, current_files[i], sep="")
        new_file_name <- paste(current_dir,current_pathways,"_enrichment_",current_files[i], sep="")
        
        if (file.exists(new_file_name)) {
            subset_enrichment_test <- readRDS(new_file_name)
            pathwaysDF <- cbind(pathwaysDF, subset_enrichment_test)
        }
        else {
            
            total_genes <- fread(current_file_path, header=FALSE)
            
            if(nrow(total_genes) < 1) {next}

            colnames(total_genes) <- c("genes")
            enrichment_test <- compute_enrichment(foreground_genes = total_genes$genes,
                                                   all_genes = background_set$gene,
                                                   pathways = pathways[[current_pathways]])
            subset_enrichment_test  <- enrichment_test[,c("q_value")]
            rownames(subset_enrichment_test) <-  enrichment_test$pathway
            
            pathwaysDF <- cbind(pathwaysDF, subset_enrichment_test)
            saveRDS(subset_enrichment_test, file = new_file_name)
        }
        
    }
    
    all_qvals <- pathwaysDF[,2:length(pathwaysDF)]
    rownames(all_qvals) <- total_enrichment_pathways
    colnames(all_qvals) <- current_files[1:length(all_qvals)]
    
    #enrichment_heatmap <- pheatmap(as.matrix(log(all_qvals)), show_rownames = F)
    return(list("matrix" = all_qvals))
}

load_enrichment <- function(directory, file_regex,current_pathways){
    current_directory_files <- list.files(directory)
    current_files <- current_directory_files[grep(file_regex, current_directory_files)]

    total_enrichment_pathways <- names(pathways[[current_pathways]])
    pathwaysDF <- data.frame("pathways" = total_enrichment_pathways)
    for(i in 1:length(current_files)){
        current_enrichment<-readRDS(paste(directory,current_files[i], sep=""))
        pathwaysDF <- cbind(pathwaysDF, current_enrichment)
    }
    all_qvals <- pathwaysDF[,2:length(pathwaysDF)]
    rownames(all_qvals) <- total_enrichment_pathways
    colnames(all_qvals) <- current_files[1:length(all_qvals)]

    # enrichment_heatmap <- pheatmap(as.matrix(log(all_qvals)), show_rownames = F)
    return(list("matrix" = all_qvals))
}

wrs_list <- list("upregulated" = list(),
                 "downregulated" = list())
recurrent_wrs_list <- list("upregulated" = list(),
                           "downregulated" = list())
recurrent_wrs_list_500 <- list("upregulated" = list(),
                           "downregulated" = list())
recurrent_wrs_list_100 <- list("upregulated" = list(),
                           "downregulated" = list())
recurrent_wrs_list_50 <- list("upregulated" = list(),
                           "downregulated" = list())


immediate_list <- list("upregulated" = list(),
                 "downregulated" = list())
recurrent_immediate_list <- list("upregulated" = list(),
                           "downregulated" = list())
recurrent_immediate_list_500 <- list("upregulated" = list(),
                           "downregulated" = list())
recurrent_immediate_list_100 <- list("upregulated" = list(),
                           "downregulated" = list())
recurrent_immediate_list_50 <- list("upregulated" = list(),
                           "downregulated" = list())

exp_wrs_list <- list("upregulated" = list(),
                 "downregulated" = list())
recurrent_exp_wrs_list <- list("upregulated" = list(),
                           "downregulated" = list())

exp_immediate_list <- list("upregulated" = list(),
                 "downregulated" = list())
recurrent_exp_immediate_list <- list("upregulated" = list(),
                           "downregulated" = list())

wrs_files <-list.files('../data/geo_degs/')
wrs_foi <- wrs_files[grep("^[a-z].*", wrs_files)]
species_exptype <- unique(sub('_[A-Z].*$', '',wrs_foi))
exptype <- unique(sub('^[a-z].*_','',species_exptype))

length(species_exptype)

library(stringr)
length(unique(sapply(wrs_foi, FUN=function(x){str_split(x, "_")[[1]][3]})))

immediate_files <-list.files('../data/immediate_genes/')
immediate_foi <- immediate_files[grep("^[a-z].*", immediate_files)]
immediate_species_exptype <- unique(sub('_[A-Z].*$', '',immediate_foi))
immediate_exptype <- unique(sub('^[a-z].*_','',immediate_species_exptype))

wrs_foi_up <- wrs_files[grep("^[a-z].*upregulated.*", wrs_files)]
wrs_foi_down <- wrs_files[grep("^[a-z].*downregulated.*", wrs_files)]

immediate_foi_up <- immediate_files[grep("^[a-z].*upregulated.*", immediate_files)]
immediate_foi_down <- immediate_files[grep("^[a-z].*downregulated.*", immediate_files)]

for(i in 1:length(species_exptype)) {
    
    current_species_up <- wrs_foi_up[grep(species_exptype[i], wrs_foi_up)]
    current_species_down <- wrs_foi_down[grep(species_exptype[i], wrs_foi_down)]
    total_spec_exp_up <- character(0)
    total_spec_exp_down <- character(0)
    
    if(length(current_species_up) == length(current_species_down)) {
        for(j in 1:length(current_species_up)) {
            current_file_path_up <- paste('../data/geo_degs/', current_species_up[j], sep="")
            current_file_path_down <- paste('../data/geo_degs/', current_species_down[j], sep="")
            
            total_genes_up <- fread(current_file_path_up, header=FALSE)$V1
            total_genes_down <- fread(current_file_path_down, header=FALSE)$V1
            
            total_spec_exp_up <- append(total_spec_exp_up,
                                        unique(total_genes_up),
                                        length(total_spec_exp_up))
            total_spec_exp_down <- append(total_spec_exp_down,
                                          unique(total_genes_down),
                                          length(total_spec_exp_down))
        }
    }
    # union of all genes
    wrs_list[["upregulated"]][[species_exptype[i]]] <- unique(total_spec_exp_up)
    wrs_list[["downregulated"]][[species_exptype[i]]] <- unique(total_spec_exp_down)
    
    # top recurrent genes
    recurrent_wrs_list[["upregulated"]][[species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_up,
                                                                                     species_exptype[i],
                                                                                     length(current_species_up),
                                                                                     "upregulated",1000)
    recurrent_wrs_list[["downregulated"]][[species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_down,
                                                                                       species_exptype[i],
                                                                                       length(current_species_down),
                                                                                     "downregulated",1000)
    recurrent_wrs_list_500[["upregulated"]][[species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_up,
                                                                                     species_exptype[i],
                                                                                     length(current_species_up),
                                                                                     "upregulated",500)
    recurrent_wrs_list_500[["downregulated"]][[species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_down,
                                                                                       species_exptype[i],
                                                                                       length(current_species_down),
                                                                                     "downregulated",500)
    
    recurrent_wrs_list_100[["upregulated"]][[species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_up,
                                                                                     species_exptype[i],
                                                                                     length(current_species_up),
                                                                                     "upregulated",100)
    recurrent_wrs_list_100[["downregulated"]][[species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_down,
                                                                                       species_exptype[i],
                                                                                       length(current_species_down),
                                                                                     "downregulated",100)
    
    recurrent_wrs_list_50[["upregulated"]][[species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_up,
                                                                                     species_exptype[i],
                                                                                     length(current_species_up),
                                                                                     "upregulated",50)
    recurrent_wrs_list_50[["downregulated"]][[species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_down,
                                                                                       species_exptype[i],
                                                                                       length(current_species_down),
                                                                                     "downregulated",50)
}

for(i in 1:length(immediate_species_exptype)) {
    
    current_species_up <- immediate_foi_up[grep(immediate_species_exptype[i], immediate_foi_up)]
    current_species_down <- immediate_foi_down[grep(immediate_species_exptype[i], immediate_foi_down)]
    total_spec_exp_up <- character(0)
    total_spec_exp_down <- character(0)
    
    if(length(current_species_up) == length(current_species_down)) {
        for(j in 1:length(current_species_up)) {
            current_file_path_up <- paste('../data/immediate_genes/', current_species_up[j], sep="")
            current_file_path_down <- paste('../data/immediate_genes/', current_species_down[j], sep="")
            
            total_genes_up <- fread(current_file_path_up, header=FALSE)$V1
            total_genes_down <- fread(current_file_path_down, header=FALSE)$V1
            
            total_spec_exp_up <- append(total_spec_exp_up,
                                        unique(total_genes_up),
                                        length(total_spec_exp_up))
            total_spec_exp_down <- append(total_spec_exp_down,
                                          unique(total_genes_down),
                                          length(total_spec_exp_down))
        }
    }
    # union of all genes
    immediate_list[["upregulated"]][[immediate_species_exptype[i]]] <- unique(total_spec_exp_up)
    immediate_list[["downregulated"]][[immediate_species_exptype[i]]] <- unique(total_spec_exp_down)
    
    # top recurrent genes
    recurrent_immediate_list[["upregulated"]][[immediate_species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_up,
                                                                                     immediate_species_exptype[i],
                                                                                     length(current_species_up),
                                                                                     "upregulated",1000)
    recurrent_immediate_list[["downregulated"]][[immediate_species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_down,
                                                                                       immediate_species_exptype[i],
                                                                                       length(current_species_down),
                                                                                     "downregulated",1000)
    recurrent_immediate_list_500[["upregulated"]][[immediate_species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_up,
                                                                                     immediate_species_exptype[i],
                                                                                     length(current_species_up),
                                                                                     "upregulated",500)
    recurrent_immediate_list_500[["downregulated"]][[immediate_species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_down,
                                                                                       immediate_species_exptype[i],
                                                                                       length(current_species_down),
                                                                                     "downregulated",500)
    
    recurrent_immediate_list_100[["upregulated"]][[immediate_species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_up,
                                                                                     immediate_species_exptype[i],
                                                                                     length(current_species_up),
                                                                                     "upregulated",100)
    recurrent_immediate_list_100[["downregulated"]][[immediate_species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_down,
                                                                                       immediate_species_exptype[i],
                                                                                       length(current_species_down),
                                                                                     "downregulated",100)
    
    recurrent_immediate_list_50[["upregulated"]][[immediate_species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_up,
                                                                                     immediate_species_exptype[i],
                                                                                     length(current_species_up),
                                                                                     "upregulated",50)
    recurrent_immediate_list_50[["downregulated"]][[immediate_species_exptype[i]]] <- get_reccurent_genes(total_spec_exp_down,
                                                                                       immediate_species_exptype[i],
                                                                                       length(current_species_down),
                                                                                     "downregulated",50)
}


for(i in 1:length(exptype)) {
    
    current_exp_up <- wrs_foi_up[grep(exptype[i], wrs_foi_up)]
    current_exp_down <- wrs_foi_down[grep(exptype[i], wrs_foi_down)]
    total_exp_up <- character(0)
    total_exp_down <- character(0)
    
    if(length(current_exp_up) == length(current_exp_down)) {
        for(j in 1:length(current_exp_up)) {
            current_file_path_up <- paste('../data/geo_degs/', current_exp_up[j], sep="")
            current_file_path_down <- paste('../data/geo_degs/', current_exp_down[j], sep="")
            
            total_genes_up <- fread(current_file_path_up, header=FALSE)$V1
            total_genes_down <- fread(current_file_path_down, header=FALSE)$V1
            
            total_exp_up <- append(total_exp_up,
                                        unique(total_genes_up),
                                        length(total_exp_up))
            total_exp_down <- append(total_exp_down,
                                          unique(total_genes_down),
                                          length(total_exp_down))
        }
    }
    # union of all genes
    exp_wrs_list[["upregulated"]][[exptype[i]]] <- unique(total_exp_up)
    exp_wrs_list[["downregulated"]][[exptype[i]]] <- unique(total_exp_down)
    
    # top recurrent genes
    recurrent_exp_wrs_list[["upregulated"]][[exptype[i]]] <- get_reccurent_genes(total_exp_up,
                                                                                     exptype[i],
                                                                                     length(current_exp_up),
                                                                                     "upregulated",1000)
    recurrent_exp_wrs_list[["downregulated"]][[exptype[i]]] <- get_reccurent_genes(total_exp_down,
                                                                                       exptype[i],
                                                                                       length(current_exp_down),
                                                                                     "downregulated",1000)
}

for(i in 1:length(immediate_exptype)) {
    
    current_exp_up <- immediate_foi_up[grep(immediate_exptype[i], immediate_foi_up)]
    current_exp_down <- immediate_foi_down[grep(immediate_exptype[i], immediate_foi_down)]
    total_exp_up <- character(0)
    total_exp_down <- character(0)
    
    if(length(current_exp_up) == length(current_exp_down)) {
        for(j in 1:length(current_exp_up)) {
            current_file_path_up <- paste('../data/immediate_genes/', current_exp_up[j], sep="")
            current_file_path_down <- paste('../data/immediate_genes/', current_exp_down[j], sep="")
            
            total_genes_up <- fread(current_file_path_up, header=FALSE)$V1
            total_genes_down <- fread(current_file_path_down, header=FALSE)$V1
            
            total_exp_up <- append(total_exp_up,
                                        unique(total_genes_up),
                                        length(total_exp_up))
            total_exp_down <- append(total_exp_down,
                                          unique(total_genes_down),
                                          length(total_exp_down))
        }
    }
    # union of all genes
    exp_immediate_list[["upregulated"]][[immediate_exptype[i]]] <- unique(total_exp_up)
    exp_immediate_list[["downregulated"]][[immediate_exptype[i]]] <- unique(total_exp_down)
    
    # top recurrent genes
    recurrent_exp_immediate_list[["upregulated"]][[immediate_exptype[i]]] <- get_reccurent_genes(total_exp_up,
                                                                                     immediate_exptype[i],
                                                                                     length(current_exp_up),
                                                                                     "upregulated",1000)
    recurrent_exp_immediate_list[["downregulated"]][[immediate_exptype[i]]] <- get_reccurent_genes(total_exp_down,
                                                                                       immediate_exptype[i],
                                                                                       length(current_exp_down),
                                                                                     "downregulated",1000)
}


