*PhenoCycler Plugins*
======================

These plugins are related to the processing of *PhenoCycler* datasets but can also 
be applied to any multiplexed immunofluorescence dataset.

CellSegmentation
-----------------

This CLI uses *CellPose* for segmentation of nuclei from the DAPI/nuclei stain channel of an incoming IF image.

Cell segmentations generated per-patch are then converted into GeoJSON formatted annotations and merged 
with those from previous patches for increased performance.


RegisterImages
---------------

Registration of an IF image with a histology image from the same section can be facilitated using this CLI. 

Nuclei are first segmented from both images at low-resolution to create a basis for alignment. For IF images 
this process involves simple thresholding and some morphological operations on the DAPI/nuclei/first channel (frame). 
For histology images, color deconvolution is used as a basis for grayscale conversion of H&E/PAS stained images prior to 
thresholding and simple morphological operations.

Once nuclei are available for both types of images, the library *pystackreg* is used to determine the optimal affine 
transform to use in order to align the two images. Conversion of the IF image to align with the histology image is facilitated 
using the *large-image-converter* library which is a subsidiary of the *large-image* library.
