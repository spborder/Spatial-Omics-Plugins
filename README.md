# *Spatial --Omics Plugins*
Plugins to support integrative analysis of structures in large histology images and their associated spatial --omics data.

This set of plugins uses *fusion-tools*, a package created for the analysis of spatial --omics data and the generation of modular interactive visualization dashboards. For more details on *fusion-tools* see the [repo](https://github.com/spborder/fusion-tools/)

## Assay types supported:

1. "general"
    - This includes multi-use plugins which can be applied to any assay as well as with non-spatial data

    1. `CreateAnnotation`
        - This plugin includes some utilities for manipulating annotations stored in DSA. Operations include combining/subtracting one annotation from another (+/-) as well as "within" and combinations of operations when passing a JSON-style list in the advanced parameters.
    
    2. `ExportAnnotations`
        - This plugin lets users export annotations to either GeoJSON or XML formats. Translated annotations will appear in the items "files".

    3. `SpatialAggregation`
        - This plugin enables the spatial aggregation of properties within other structures. Aggregated properties are derived from intersecting structures' "user" field and will be transferred to the overlapping structure's (or structures') "user" field.

2. "codex"
    - This includes plugins for processing *CODEX* (Co-Detection by Indexing) or *MxIF* (Multiplexed Immunofluorescence) slides.

    1. `CellSegmentation`
        - This plugin includes CellPose for cell segmentation. Updating this plugin for other methods of cell segmentation should be straightforward as the patching->prediction->recombination procedure should be the same regardless of the model.

    2. `RegisterImages`
        - This plugin allows for registration of same-section/adjacent-section histology based on segmented nuclei. Should also allow for uploading affine transform matrices for other methods of alignment.

3. "visium"
    - This includes plugins for *10x Visium* Spatial Transcriptomics.

    1. `CellDeconvolution`
        - This plugin implements cell composition deconvolution using *Azimuth* (*Seurat*) as well as alignment using the snRNA-seq reference created by the Kidney Precision Medicine Project (KPMP)

    2. `SpotAnnotation`
        - This plugin creates the "spot" annotations using known properties of *Visium* spots (55um diameter, 100um apart). Slide microns-per-pixel (MPP) values are automatically derived using this relationship. Users can define what properties are automatically included for quick reference in each spot.

4. "visium_hd"
    - This includes plugins for *10x VisiumHD*.

5. "xenium"
    - This includes plugins for *10x Xenium*.





    




