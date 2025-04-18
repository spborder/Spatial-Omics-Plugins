*Visium Plugins*
=================

This set of CLIs is for processing data from *10x Visium* datasets.

This includes creation of spot annotations and cell deconvolution using a variety of methods.

CellDeconvolution
------------------

This CLI uses spot-level transcriptomics data to derive cell composition per-spot. Several different methods 
for cell deconvolution are implemented for *10x Visium* sections including *LabelTransfer* (*Seurat v4.0*), *Azimuth*, 
and *STDeconvolve*. Each of these deconvolution techniques are executed as R scripts with a set of command-line inputs.

The resulting deconvolved data file is then posted to the original item (image) to which the raw data is also attached.


SpotAnnotation
---------------

Spot annotations are created from the outputs of the CellDeconvolution CLI. Dataframes created as a result of cell deconvolution 
are aligned with the index of another CSV file containing the centroids of each spot. As some properties of the 
spots are already known (55um diameter and 100um apart), the pixel radius that each spot should be buffered by is 
readily calculable.

The *shapely* library is used to create circular ROIs from *Point* geometries provided in the centroid coordinate file. Data
that is extracted from the cell deconvolution output file are then added to the "properties" field of each Spot following GeoJSON 
convention before conversion to *large-image* format and POST-ing to the original item (image).


