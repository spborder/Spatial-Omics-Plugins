## Running cell deconvolution methods
library(Seurat)
library(STdeconvolve)
library(stringr)
library(SeuratDisk)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(dplyr)
library(tools)
options(timeout=300)

# Function for reading in different types of counts files
read_data_formats <- function(input_file_path){
    file_extension <- file_ext(input_file_path)

    if (tolower(file_extension) == "rds"){
        read_file <- readRDS(input_file_path)
    } else if (tolower(file_extension) == "h5"){
        read_file <- Read10X_h5(input_file_path)
        read_file <- CreateSeuratObject(counts = read_file)
    } else if (tolower(file_extension) == "h5ad"){
        read_file <- LoadH5Seurat(input_file_path)
    }

    return(read_file)
}

# Function for running STdeconvolve
RunSTDeconvolve <- function(read_input_file){

    counts <- read_input_file@assays[[read_input_file@active.assay]]$counts

    # Using default parameters from their GitHub
    counts <- cleanCounts(counts,min.lib.size=100)
    ## feature select for genes
    corpus <- restrictCorpus(counts,removeAbove=1.0,removeBelow=0.05)
    ## choose optimal number of cell-types
    ldsas <- fitLDA(t(as.matrix(corpus)),Ks=seq(2,9,by=1))
    ## getting best model results
    optLDA <- optimalModel(models=ldas,opt="min")
    ## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
    results <- getBetaTheta(optLDA,perc.filt = 0.05, betaScale = 1000)
    deconProp <- results$theta
    deconGexp <- results$beta

    # Modifying column names in deconProp
    colnames(deconProp) <- lapply(colnames(deconProp),function(i){paste("ST Topic",i,sep=" ")})

    read_input_file@assays[["stdeconvolve_results"]] <- CreateAssayObject(data=deconProp)

    return(read_input_file)
}

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


# General function for getting deconvolution results
get_cell_deconvolution <- function(input_file, organ_key){

    # Reading input file
    read_input_file <- read_data_formats(input_file)
    file_extension <- file_ext(input_file)

    if (organ_key == "kidneykpmp"){
        print("Using KPMP Reference")
        integrated_spatial_data <- integrate_kpmp_atlas(read_input_file)
    } else if (organ_key =="st_deconvolve"){
        print("Using STdeconvolve")
        integrated_spatial_data <- RunSTDeconvolve(read_input_file)
    } else {
        integrated_spatial_data <- RunAzimuth(read_input_file,organ_key)
    }

    output_path <- str_replace(input_file,paste(".",file_extension,sep=""),'_integrated.rds')
    saveRDS(integrated_spatial_data,output_path)
}


arg_list <- commandArgs(trailingOnly=TRUE)
input_file <- arg_list[1]
organ_key <- arg_list[2]

get_cell_deconvolution(input_file,organ_key)

