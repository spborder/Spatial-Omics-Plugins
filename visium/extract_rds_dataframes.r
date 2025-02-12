#!/usr/bin/env Rscript
library(Seurat)
library(stringr)

key_list <- commandArgs(trailingOnly=TRUE)
input_file <- key_list[1]
key_list <- key_list[2:length(key_list)]

read_input_file <- readRDS(input_file)

# Extracting dataframes from RDS file and saving as csvs in the temporary directory
for (k in key_list){
    print(k)
    if (k %in% names(read_input_file@assays)){
        write.csv(read_input_file[[k]]@data,paste(gsub("\\.","_",k),".csv",sep=""))
    }
}

# Writing spot coordinates
write.csv(read_input_file$images[["slice1"]]@coordinates@data,"./spot_coordinates.csv")
