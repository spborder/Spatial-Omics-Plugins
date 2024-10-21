"""Plugin implementation of cell segmentation methods
"""

import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
import json

from fusion_tools.dataset import SegmentationDataset
from fusion_tools.utils.shapes import export_annotations
import large_image
from tqdm import tqdm

import cellpose
import deepcell
import tensorflow as tf

from pathlib import Path

import rasterio
import rasterio.features

import girder_client
from ctk_cli import CLIArgumentParser


class CellModel:
    def __init__(self,
                 model_type: str,
                 model_args: dict):
        
        self.model_type = model_type
        self.model_args = model_args

        if self.model_type=='cellpose':
            self.model = cellpose.models.CellposeModel(pretrained_model = '../.cellpose/models/tissuenet_cp3',gpu=True)
        elif self.model_type=='deepcell':
            self.model = deepcell.applications.NuclearSegmentation(
                model = tf.keras.models.load_model(
                    Path("../.deepcell/models") / 'NuclearSegmentation'
                )
            )

    def predict(self, input_image, bbox):

        if self.model_type=='cellpose':
            masks, _, _ = self.model.eval(input_image,
                                          diameter = self.model_args['diameter'],
                                          cellprop_threshold = self.model_args['cellprob_threshold'],
                                          channels = self.model_args['channels'])
        
        elif self.model_type=='deepcell':
            masks = self.model.predict(input_image)

        mask_shapes = self.mask_to_shape(masks,bbox)

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
                        for i in geo['coordinates']
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
    
    cell_seg_dataset = SegmentationDataset(
        slides = [],
        annotations = None,
        use_parallel = False,
        verbose = True,
        patch_mode = 'all',
        patch_region = 'tissue',
        patch_size = [224,224]
    )

    cell_model = CellModel(
        model_type = args.model_type,
        model_args = args.model_args
    )

    all_cells_gdf = gpd.GeoDataFrame()

    with tqdm(cell_seg_dataset, total = len(cell_seg_dataset)) as pbar:
        pbar.set_description('Predicting on patches in SegmentationDataset')
        for idx, (patch,_) in enumerate(cell_seg_dataset):
            
            bbox = cell_seg_dataset.data[idx]['bbox']
            mask_features = cell_model.predict(patch,bbox)

            if len(mask_features)>0:
                # Converting masks to annotations
                if all_cells_gdf.empty:
                    all_cells_gdf = gpd.GeoDataFrame.from_features(mask_features)

                else:
                    new_cells = gpd.GeoDataFrame.from_features(mask_features)
                    all_cells_gdf = pd.concat([all_cells_gdf,new_cells],axis=0,ignore_index=True)
                    merged_geoms = all_cells_gdf.union_all().geoms
                    all_cells_gdf = gpd.GeoDataFrame({'geometry': merged_geoms, 'name': ["Segmented Cells"]*len(merged_geoms)})

            pbar.update(1)

    pbar.close()   

    annotations_geo = all_cells_gdf.to_geo_dict(show_bbox=True)
    export_annotations(annotations_geo,format='histomics',save_path='./new_annotations.json')

    with open('./new_annotations.json','r') as f:
        new_annotations = json.load(f)

        f.close()



    gc.post(
        f'/annotation/item/{image_id}?token={args.girderToken}',
        data = json.dumps(new_annotations),
        headers = {
            'X-HTTP-Method': 'POST',
            'Content-Type': 'application/json'
        }
    )






if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

