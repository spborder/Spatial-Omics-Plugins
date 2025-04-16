"""Plugin implementation of cell segmentation methods
"""

import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
import json

#from fusion_tools.dataset import SegmentationDataset
from fusion_tools.utils.shapes import export_annotations
import large_image
from tqdm import tqdm

from cellpose import models
#import deepcell
#import tensorflow as tf

from pathlib import Path

import rasterio
import rasterio.features

import girder_client
from ctk_cli import CLIArgumentParser

# SegmentationDataset-specific imports
import random
from math import floor

import large_image.constants
import numpy as np
import pandas as pd
import geopandas as gpd
import large_image
from tqdm import tqdm

from shapely.geometry import shape, box, Polygon
from shapely.validation import make_valid
from shapely.ops import unary_union
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_holes
from skimage.measure import label, find_contours
from skimage.draw import polygon2mask
from PIL import Image
from io import BytesIO

from typing import Callable
from typing_extensions import Union

from joblib import Parallel, delayed
from fusion_tools.dataset import SegmentationDataset
from fusion_tools.utils.shapes import geojson_to_histomics

class CellModel:
    def __init__(self,
                 model_type: str,
                 model_args: dict):
        
        self.model_type = model_type.lower().strip()
        self.model_args = model_args

        if self.model_type=='cellpose':
            import torch

            if torch.cuda.is_available() and self.model_args['use_gpu']:
                use_gpu = True
            else:
                use_gpu = False

            self.model = models.CellposeModel(pretrained_model = '../.cellpose/models/cyto3',gpu=use_gpu)
        elif self.model_type=='deepcell':
            raise NotImplementedError
            """
            self.model = deepcell.applications.NuclearSegmentation(
                model = tf.keras.models.load_model(
                    Path("../.deepcell/models") / 'NuclearSegmentation'
                )
            )
            """
        else:
            raise NotImplementedError
        
        print(f'Model Loaded: {self.model_type}')

    def predict(self, input_image, bbox):

        masks = None
        if self.model_type=='cellpose':
            input_image = input_image.astype(np.uint8)
            masks,flows,styles = self.model.eval([input_image],channel_axis=2)
        
        elif self.model_type=='deepcell':
            raise NotImplementedError
            #masks = self.model.predict(input_image)
        
        if not masks is None:
            mask_shapes = []
            for m in masks:
                mask_shapes.extend(self.mask_to_shape(m,bbox))
        else:
            print('Masks is None')
            mask_shapes = []

        return mask_shapes

    def mask_to_shape(self, mask: np.ndarray, bbox:list)->list:
        mask_features = []
        
        for geo,val in rasterio.features.shapes(mask,mask=mask>0):
            mask_features.append({
                'type': 'Feature',
                'geometry': {
                    'type': geo['type'],
                    'coordinates': [[
                        [float(i[0]+bbox[0]),float(i[1]+bbox[1])]
                        for i in geo['coordinates'][0]
                    ]]
                },
                'properties': {
                    'name': 'Segmented Cell'
                }
            })
        
        return mask_features


def main(args):

    sys.stdout.flush()

    # Initialize girder client
    gc = girder_client.GirderClient(
        apiUrl = args.girderApiUrl
    )
    gc.setToken(args.girderToken)

    print('Input arguments: ')
    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

    file_info = gc.get(f'/file/{args.input_image}')
    image_id = file_info['itemId']
    image_name = file_info['name']

    region_annotation = None
    if args.input_region==[-1,-1,-1,-1]:
        patch_region = 'all'
    else:
        patch_region = {
            'type': 'FeatureCollection',
            'features': [
                {
                    'type': 'Feature',
                    'geometry': {
                        'type': 'Polygon',
                        'coordinates': [[
                            [args.input_region[0],args.input_region[1]],
                            [args.input_region[2],args.input_region[1]],
                            [args.input_region[2],args.input_region[3]],
                            [args.input_region[0],args.input_region[3]],
                            [args.input_region[0],args.input_region[1]]
                        ]]
                    },
                    'properties': {
                        'name': 'Segmentation Region'
                    }
                }
            ],
            'properties': {
                'name': 'Segmentation Region'
            }
        }

        if args.return_segmentation_region:
            region_annotation = geojson_to_histomics(patch_region)



    if args.use_frame_index:
        seg_transform = lambda img: img[:,:,args.segmentation_frame][:,:,None]
    else:
        seg_transform = None

    # Downloading the slide 
    gc.downloadFile(
        fileId = args.input_image,
        path = os.getcwd()+f'/{image_name}'
    )
    
    cell_seg_dataset = SegmentationDataset(
        slides = [
            os.getcwd()+f'/{image_name}'
        ],
        annotations = None,
        use_parallel = False,
        verbose = True,
        patch_mode = 'all',
        patch_region = patch_region,
        patch_size = [args.patch_size,args.patch_size],
        transforms=seg_transform,
        shuffle = False
    )

    # This can be added to the xml some other time
    model_args = {
        'diameter': None,
        'cellprob_threshold': 0.0,
        'channels': [0,0],
        'use_gpu': args.use_gpu
    }
    cell_model = CellModel(
        model_type = args.method,
        model_args = model_args
    )

    all_cells_gdf = gpd.GeoDataFrame()

    cell_count = 0
    bbox_list = []
    with tqdm(cell_seg_dataset, total = len(cell_seg_dataset)) as pbar:
        for idx, (patch,_) in enumerate(cell_seg_dataset):
            
            bbox = cell_seg_dataset.data[idx]['bbox']
            bbox_list.append({
                'type': 'Feature',
                'geometry': {
                    'type': 'Polygon',
                    'coordinates': [list(box(*bbox).exterior.coords)],
                },
                'properties': {
                    'bbox_idx': idx
                }
            })
            mask_features = cell_model.predict(patch,bbox)
            cell_count=len(mask_features)
            if len(mask_features)>0:
                # Converting masks to annotations
                if all_cells_gdf.empty:
                    all_cells_gdf = gpd.GeoDataFrame.from_features(mask_features)

                else:
                    new_cells = gpd.GeoDataFrame.from_features(mask_features)
                    all_cells_gdf = pd.concat([all_cells_gdf,new_cells],axis=0,ignore_index=True)
                    merged_geoms = all_cells_gdf.union_all().geoms
                    all_cells_gdf = gpd.GeoDataFrame({'geometry': merged_geoms, 'name': ["Segmented Cells"]*len(merged_geoms)})

            pbar.set_description(f'Predicting on patches in SegmentationDataset: {cell_count} in last mask, {all_cells_gdf.shape[0]} post-merge')
            pbar.update(1)

    pbar.close()   

    if not all_cells_gdf.empty:
        annotations_geo = all_cells_gdf.to_geo_dict(show_bbox=True)
        annotations_geo['properties'] = {
            'name': 'Segmented Cells'
        }

        print(f'Found: {len(annotations_geo["features"])} Cells!')
        new_annotations = geojson_to_histomics(annotations_geo)

        gc.post(
            f'/annotation/item/{image_id}?token={args.girderToken}',
            data = json.dumps(new_annotations),
            headers = {
                'X-HTTP-Method': 'POST',
                'Content-Type': 'application/json'
            }
        )

        if args.return_segmentation_region and not region_annotation is None:
            bboxes = geojson_to_histomics(
                {
                    'type': 'FeatureCollection',
                    'properties': {
                        'name': 'Segmentation Patches'
                    },
                    'features': bbox_list
                }
            )

            gc.post(
                f'/annotation/item/{image_id}?token={args.girderToken}',
                data = json.dumps(bboxes),
                headers = {
                    'X-HTTP-Method': 'POST',
                    'Content-Type': 'application/json'
                }
            )

            gc.post(
                f'/annotation/item/{image_id}?token={args.girderToken}',
                data = json.dumps(region_annotation),
                headers = {
                    'X-HTTP-Method': 'POST',
                    'Content-Type': 'application/json'
                }
            )


    # Putting job metadata
    job_submitter = gc.get('/user/me')
    job_meta = {
        'patch_size': args.patch_size,
        'use_gpu': args.use_gpu,
        'segmentation_frame': args.segmentation_frame,
        'use_frame_index': args.use_frame_index,
        'input_region': args.input_region,
        'user': job_submitter['login']
    }

    gc.put(f'/item/{args.input_image}/metadata',parameters={'metadata': job_meta})


if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

