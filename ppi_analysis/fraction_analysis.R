#!/usr/bin/env Rscript



args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] <- "out.txt"
}

suppressMessages(suppressWarnings(library(data.table)))

euclidean <- function(x1,y1,x2,y2) {
    return(sqrt((x2-x1)^2+(y2-y1)^2))
}
closeness_measures <-function(df) {
    pairwise_euclidean_distances <- numeric()
    start <- as.numeric(args[3])
    end <-  as.numeric(args[4])
    for(i in start:end) {
        for(j in i:nrow(df)) {
            if(j != i) {
                pairwise_euclidean_distances <- append(pairwise_euclidean_distances,
                                                       euclidean(df[i,1],df[i,2],df[j,1],df[j,2]),
                                                       length(pairwise_euclidean_distances))
            }
        }
    }
    
    return(pairwise_euclidean_distances)
}

current_file <- readRDS(args[1])

output <- closeness_measures(current_file)

saveRDS(output, paste("/data/timonaj/cancer_as_wound/ppi_analysis/pairwise_euclidean/",args[2],sep=""))
