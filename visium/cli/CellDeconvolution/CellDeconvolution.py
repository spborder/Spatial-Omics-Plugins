"""Implementing cell composition deconvolution
"""
import os
import sys

from ctk_cli import CLIArgumentParser
import girder_client

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri


ORGAN_REF_KEY = {
    "Azimuth Adipose Reference": "adiposeref",
    "Azimuth Bone Marrow Reference": "bonemarrowref",
    "Azimuth Fetus Reference": "fetusref",
    "Azimuth Heart Reference": "heartref",
    "Azimuth Human Cortex Reference": "humancortexref",
    "Azimuth Kidney Reference": "kidneyref",
    "KPMP Atlas Kidney": "kidneykpmp",
    "Azimuth Lung Reference": "lungref",
    "Azimuth Pancreas Reference": "pancreasref",
    "Azimuth Mouse Cortex Reference": "mousecortexref",
    "Azimuth PBMC Reference": "pbmcref",
    "Azimuth Tonsil Reference": "tonsilref"
}


INTEGRATION_DATA_KEYS = {
    'adiposeref': ['celltype.l1','celltype.l2'],
    'bonemarrowref': ['celltype.l1','celltype.l2'],
    'fetusref': ['annotation.l2','annotation.l1'],
    'heartref': ['celltype.l1','celltype.l2'],
    'humancortexref': ['class','subclass','cluster','cross-species cluster'],
    'kidneyref': ['annotation.l1','annotation.l2','annotation.l3'],
    'kidneykpmp': ['subclass.l2','subclass.l1'],
    'lungref': ['annotation.l1','annotation.l2'],
    'pancreasref': ['celltype.l1','celltype.l2'],
    'mousecortexref': ['class','subclass','cluster','cross-species cluster'],
    'pbmcref': ['celltype.l1','celltype.l2','celltype.l3'],
    'tonsilref': ['celltype.l1','celltype.l2']
}

# pancreas might be in there also ['annotation.l1']
# liver might be in there also ['celltype.l1','celltype.l2']
# tonsil has a v2, ['celltype.l1','celltype.l2']
# mouse pansci = ['Main_cell_type']
# lung has a v1 and v2, v2 = ['ann_level_1','ann_level_2','ann_level_3','ann_level_4','ann_level_5','ann_finest_level']

robjects.r('library(Seurat)')
robjects.r('library(stringr)')
robjects.r('library(SeuratDisk)')
robjects.r('library(Azimuth)')
robjects.r('library(SeuratData)')
robjects.r('library(patchwork)')
robjects.r('library(dplyr)')
robjects.r('library(tools)')
robjects.r('options(timeout=300)')

robjects.r('''
            # Function for reading in different types of files
            read_data_formats <- function(input_file_path){
                print(input_file_path)
                file_extension <- file_ext(input_file_path)
                print(file_extension)
                if (tolower(file_extension) == "rds"){
                    read_file <- readRDS(input_file_path)
                } else if (tolower(file_extension) == "h5"){
                    read_file <- Read10X_h5(input_file_path)
                } else if (tolower(file_extension) == "h5ad"){
                    read_file <- LoadH5Seurat(input_file_path)
                }
           
                return(read_file)
            }
           ''')

robjects.r('''
            # General function for integration using reference object
           integrate_spatial <- function(input_file, organ_key){

                # Reading input file
                read_input_file <- read_data_formats(input_file)
                file_extension <- file_ext(input_file)
           
                if (organ_key == "kidneykpmp"){
                    print("Using KPMP Reference")
                    integrated_spatial_data <- integrate_kpmp_atlas(read_input_file)
                } else {
                    integrated_spatial_data <- RunAzimuth(read_input_file,organ_key)
                }

                output_path <- str_replace(input_file,paste(".",file_extension,sep=""),'_integrated.rds')
                saveRDS(integrated_spatial_data, output_path)
           }
           ''')

robjects.r('''
            # Function for integration using KPMP atlas
            integrate_kpmp_atlas <- function(spatial){
                DefaultAssay(spatial) <- "SCT"
                
                atlas_path <- "../KidneyAtlas_snCV3_20percent.h5Seurat"
                kpmp_atlas <- LoadH5Seurat(atlas_path, assays = c("counts","scale.data"),tools = TRUE,images=False)
           
                Idents(kpmp_atlas) <- kpmp_atlas@meta.data$subclass.l2
           
                kpmp_atlas <- subset(kpmp_atlas, idents = "NA", invert = T)
                kpmp_atlas <- UpdateSeuratObject(kpmp_atlas)
                kpmp_atlas[["RNA"]] <- as(object = kpmp_atlas[["RNA"]],Class="SCTAssay")

                DefaultAssay(kpmp_atlas) <- "RNA"
                Idents(kpmp_atlas) <- kpmp_atlas@meta.data[["subclass.l2"]]
           
                anchors <- FindTransferAnchors(
                    reference = kpmp_atlas, query = spatial, normalization.method = "SCT",
                    query.assay = "SCT", recompute.residuals = FALSE
                )

                predictions.assay <- TransferData(
                    anchorset = anchors, refdata = kpmp_atlas@meta.data[["subclass.l2"]],
                    prediction.assay = TRUE,
                    weight.reduction = spatial[["pca"]], dims = 1:30
                )
                spatial[["pred_subclass_l2"]] <- predictions.assay
           
                df_pred <- predictions.assay@data
                max_pred <- apply(df_pred, 2, function(x) max.col(t(x),"first"))
                max_pred_val <- apply(df_pred, 2, function(x) max(t(x)))
                max_pred <- as.data.frame(max_pred)
                max_pred$Seurat_subset <- rownames(df_pred)[max_pred$max_pred]
                max_pred$score <- max_pred_val
                max_pred$Barcode <- rownames(max_pred)
           
                spatial@meta.data.subclass.l2 <- max_pred$Seurat_subset
                spatial@meta.data$subclass.l2_score <- max_pred$score
           
                Idents(kpmp_atlas) <- kpmp_atlas@meta.data[["subclass.l1"]]
           
                anchors <- FindTransferAnchors(
                    reference = kpmp_atlas, query = spatial, normalization.method = "SCT",
                    query.assay = "SCT", recompute.residuals = FALSE
                )
                predictions.assay <- TransferData(
                    anchorset = anchors, refdata = kpmp_atlas@meta.data[["subclass.l1"]],
                    prediction.assay = TRUE,
                    weight.reduction = spatial[["pca"]],dims = 1:30
                )
           
                spatial[["pred_subclass_l1"]] <- predictions.assay
           
                df_pred <- predictions.assay@data
                max_pred <- apply(df_pred, 2, function(x) max.col(t(x), "first"))
                max_pred_val <- apply(df_pred,2, function(x) max(t(x)))
           
                max_pred <- as.data.frame(max_pred)
                max_pred$Seurat_subset <- rownames(df_pred)[max_pred$max_pred]
                max_pred$score <- max_pred_val
                max_pred$Barcode <- rownames(max_pred)
           
                spatial@meta.data$subclass.l1 <- max_pred$Seurat_subset
                spatial@meta.data$subclass.l1_score <- max_pred$score
           
                return(spatial)
            }
           ''')


def main(args):

    gc = girder_client.GirderClient(
        apiUrl = args.girderApiUrl
    )
    gc.setToken(args.girderToken)

    print('Input arguments:')
    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

    if not args.organ == 'Not Listed':
        # print contents of current working directory, see if files were copied over
        print('Contents of working directory')
        print(os.listdir(os.getcwd()+'/'))
        file_info = gc.get(f'/file/{args.counts_file}')

        # Downloading counts file to cwd
        gc.downloadFile(
            args.counts_file,
            path = f'{os.getcwd()}/{file_info["name"]}'
        )
        print('Updated contents of directory')
        print(os.listdir(os.getcwd()+'/'))

        print(f'Running cell deconvolution for: {args.organ}')
        integrator = robjects.globalenv['integrate_spatial']
        integrator(
            file_info['name'],
            ORGAN_REF_KEY[args.organ]
        )

        print(os.listdir(os.getcwd()+'/'))

        print(f'Uploading file to {file_info["itemId"]}')
        # Posting integration results to item
        file_ext = file_info['name'].split('.')[-1]
        uploaded_file = gc.uploadFileToItem(
            itemId = file_info['itemId'],
            filepath = f'./{file_info["name"].replace("."+file_ext,"_integrated.rds")}'
        )
        
if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

