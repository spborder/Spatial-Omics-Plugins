"""Plugin for programmatic feature extraction
"""
import os
import sys
import numpy as np
import json

import girder_client
from ctk_cli import CLIArgumentParser

from fusion_tools.feature_extraction import (
    ParallelFeatureExtractor, distance_transform_features, color_features,
    morphological_features, texture_features
)

import geopandas as gpd
from shapely.geometry import shape
from fusion_tools.tileserver import DSATileServer
from fusion_tools.handler import DSAHandler

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
    
    image_id = gc.get(f'/file/{args.input_image}')['itemId']

    tile_server = DSATileServer(
        api_url = args.girderApiUrl,
        item_id = image_id
    )

    annotations = DSAHandler(girderApiUrl=args.girderApiUrl).get_annotations(
        item = image_id
    )

    feature_extractor = ParallelFeatureExtractor(
        image_source = tile_server,
        feature_list = [
            lambda image, mask, coords: distance_transform_features(image, mask, coords),
            lambda image, mask, coords: color_features(image, mask, coords),
            lambda image, mask, coords: texture_features(image, mask, coords),
            lambda image, mask, coords: morphological_features(image, mask, coords)
        ],
        preprocess = None,
        n_jobs = 4,
        verbose = True
    )

    # Add these dataframes as files? Or store in the "user" field?
    histomics_anns = []
    for ann_idx,ann in enumerate(annotations):
        print(f'On annotation: {ann_idx+1}/{len(annotations)}')
        feature_df = feature_extractor.start(ann['features'])

        histomics_ann = {
            'annotation': {
                'name': ann['properties']['name'],
                'elements': []
            }
        }

        # Feature_df will have bounding box columns (min_x, min_y, max_x, max_y) which can be used to align with structures
        if args.save_to_files:
            save_path = os.getcwd()+f'\{ann["properties"]["name"]}.csv'
            feature_df.to_csv(save_path)
            gc.uploadFileToItem(image_id, save_path, reference=None, mimeType=None, filename=None, progressCallback=None)
        
        if args.save_to_elements:
            # They should just match by index but this will double check to ensure that bounding boxes align between Features (in the geojson) and feature rows in the dataframe
            bbox_coords = feature_df[['min_x','min_y','max_x','max_y']].values.tolist()

            for f in ann['features']:
                ann_bbox = list(shape(f['geometry']).bounds)
                feature_df_idx = bbox_coords.index(ann_bbox)

                f['properties'] = f['properties'] | feature_df.iloc[feature_df_idx,:].to_dict()

                histomics_el = {
                    'type': 'polyline',
                    'points': [i+[0] for i in f['geometry']['coordinates'][0]],
                    'user': f['properties']
                }

                histomics_ann['annotation']['elements'].append(histomics_el)

        histomics_anns.append(histomics_ann)
            
    if args.save_to_elements:
        # Just posting a copy for now
        gc.post(
            f'/annotation/item/{image_id}?token={args.girderToken}',
            data = json.dumps(histomics_anns),
            headers = {
                'X-HTTP-Method': 'POST',
                'Content-Type': 'application/json'
            }
        )




if __name__=='__main__':
    main(CLIArgumentParser().parse_args())

