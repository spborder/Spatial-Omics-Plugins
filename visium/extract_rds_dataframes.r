#!/usr/bin/env Rscript
library(Seurat)
library(stringr)

key_list <- commandArgs(trailingOnly=TRUE)
input_file <- key_list[1]
key_list <- key_list[2:length(key_list)]

read_input_file <- readRDS(input_file)

input_file_path <- unlist(str_split(input_file,'/'))
input_file_path <- paste(input_file_path[1:length(input_file_path)-1],collapse='/')
print(paste("input_file_path:",input_file_path,sep=""))

# Extracting dataframes from RDS file and saving as csvs in the temporary directory
for (k in key_list){
    print(paste("key:",k,sep=""))
    if (k %in% names(read_input_file@assays)){
        save_path <- paste(input_file_path,paste(gsub("\\.","_",k),".csv",sep=""),sep="/")
        print(paste("save_path:",save_path,sep=""))
        write.csv(read_input_file[[k]]@data,save_path)
    }
}

# Writing spot coordinates
spot_save_path <- paste(input_file_path,"spot_coordinates.csv",sep='/')
print(paste("spot_save_path: ",spot_save_path,sep=''))

if ("coordinates" %in% names(read_input_file@images$slice1)){
    # This is VisiumV1 format
    write.csv(read_input_file@images[["slice1"]]@coordinates,spot_save_path)
} else if ("centroids" %in% names(read_input_file@images$slice1)){
    # This is VisiumV2 format
    centroids <- as.data.frame(read_input_file@images$slice1$centroids@coords)
    rownames(centroids) <- read_input_file@images$slice1$centroids@cells

    write.csv(centroids,spot_save_path)
}
