#!/usr/bin/env Rscript
library(Seurat)
library(stringr)
library(tidyverse)

arg_list <- commandArgs(trailingOnly=TRUE)
input_file <- arg_list[1]

gene_selection_method <- arg_list[2]
print(paste('Using gene selection method:',gene_selection_method,sep=""))

if (gene_selection_method %in% c("vst","mean.var.plot","dispersion")) {
    gene_selection_n <- as.numeric(arg_list[3])
} else {
    gene_name_file <- read.csv(arg_list[3])
}

read_input_file <- readRDS(input_file)

input_file_path <- unlist(str_split(input_file,'/'))
input_file_path <- paste(input_file_path[1:length(input_file_path)-1],collapse='/')


if (gene_selection_method %in% c("vst","mean.var.plot","dispersion")) {
    variable_genes <- HVFInfo(object = read_input_file, method = gene_selection_method)

    # Sorting from highest to lowest
    #top_n_genes <- variable_genes |>
    #                base::sort_by(~ list(-variable_genes$residual_variance)) |>
    #                head(gene_selection_n)

    top_n_genes <- variable_genes %>%
        arrange(desc(variable_genes$residual_variance)) %>%
        slice(1:gene_selection_n)

    selected_names <- rownames(top_n_genes)

} else if (gene_selection_method %in% c("moransi","markvariogram")) {
    print('Spatial methods not implemented for this assay')
} else {
    # Checking first column
    print(gene_name_file[1:5,])
    selected_names <- as.character(gene_name_file[ ,1])

}

extract_data <- FetchData(read_input_file,vars=selected_names)

write.csv(extract_data,paste(input_file_path,"selected_genes.csv",sep="/"))
